############################################################
##                                                        ##
##               equispecseq.gi                      ##
## HAP subpackage for GAP (Groups Algorithms Programming) ##
##         under the GNU GPL license (v. 3),  2012        ##
##      by Alexander D. Rahm & Bui Anh Tuan               ##
############################################################



InstallGlobalFunction("EquivariantSpectralSequencePage", function( C, n)
#########################################################################

#########################################################################
local reducedTorsionCells,celldata,j,N,l,w,P,Q,RP,RQ,g,Pt,k,pos,E1PageRec,p,q,
stabgrp,cell,i, E1page, T, EnPage, Differential, CohomologyOfGroup, stabres,stabcohom, inclusionMaps, groupname,name,sb,se,
maps,map,eqmap,tmp,BI,SGN,LstEl,s,r,t,multiple,E2page,CH1,Mat1, MatRank
;

if IsString(C) then
groupname:=Filtered(C,x->not(x='(' or x=')' or x=',' or x='[' or x=']'));

   Read(Concatenation( 	DirectoriesPackageLibrary("HAP")[1]![1],
			"Perturbations/Gcomplexes/",groupname));
   celldata := StructuralCopy(HAP_GCOMPLEX_LIST);
   name:=StructuralCopy(groupname);
   RemoveCharacters(name,"torsion");
   sb:=Position(name,'_');
   se:=Length(name);
   l:=Int(name{[sb+1..se]});

   reducedTorsionCells:=[];
   for i in [1..Size(celldata)] do
      reducedTorsionCells[i]:=[];
      for j in [1..Size(celldata[i])] do
          reducedTorsionCells[i][j]:=[i-1,j];
      od;
   od;
   N:=Size(reducedTorsionCells);
#   if N>2 then return fail;fi;

else


#if not IsHapTorsionSubcomplex(C) then
#    return fail;
#else
    N:=Size(C!.reducedTorsionCells);

    ## We only consider the case when length of the subcomplex is
    ## less than 2 in order to get rid of the differential d2
#    if N>2 then return fail;fi;

    reducedTorsionCells:=C!.torsionCells;
    celldata:=C!.celldata;
    l:=C!.torsion;
    groupname:=C!.groupname;
fi;
############### LIST THE STABILIZERS#######################
stabgrp:=[];
#stabres:=[];
#stabcohom:=[];
for j in [1..N] do
    stabgrp[j]:=[];

    for i in [1..Size(reducedTorsionCells[j])] do
        cell:=reducedTorsionCells[j][i];
        Add(stabgrp[j],celldata[cell[1]+1][cell[2]]!.TheMatrixStab);
    od;
od;
###################### Rank of a matrix ###################
MatRank:=function(g)
if not (IsBound(g[1]) and IsBound(g[1][1])) then return 0;
else return RankMat(g);fi;
end;
###################### E_1 Page of Cohomology##############
E1page:=function(p,q)
local w,i;
    if p=0 then return Differential(1,0,q)[1];
    else return Differential(1,p-1,q)[2];
    fi;
end;
######################End of E_1 Page######################
E2page:=function(p,q)
local t,M,lnth,N;
lnth:=Length(reducedTorsionCells);

if (p<0) or (p>lnth-1) then return 0;fi;

if p=0 then
    M:=Differential(1,0,q)[3];
    if IsEmpty(M) then return 0;
    else return Size(M[1])-RankMat(M);
    fi;
fi;
if p=lnth-1 then
    M:=Differential(1,p-1,q)[3];
    if IsEmpty(M) then return E1page(p,q);
    else return E1page(p,q)-RankMat(M);
    fi;
fi;

M:=Differential(1,p-1,q)[3];
N:=Differential(1,p,q)[3];

return E1page(p,q)-MatRank(M)-MatRank(N);

end;
##############End of E_2 Page######################
###################### Cohomology of Group ################
CohomologyOfGroup:=function(k)
local p,w,lnth;
Print("\n Users please self-aware that this attribute only works correctly in those cases whose d2 differentials are trivial \n");
    w:=0;
    lnth:=Length(reducedTorsionCells);
    for p in [0..Minimum(k,lnth-1)] do
        w:=w+E2page(p,k-p);
    od;
    return w;
return fail;
end;
######################End of Cohomology####################

            inclusionMaps:=[];
            maps:=[];
            multiple:=[];
#            N:=2;
      for k in [1..N-1] do
            maps[k]:=[];
            inclusionMaps[k]:=[];
            multiple[k]:=[];
            for i in [1..Size(reducedTorsionCells[k+1])] do
                cell:= reducedTorsionCells[k+1][i];
                maps[k][i]:=[];
                inclusionMaps[k][i]:=[];
                multiple[k][i]:=[];

                tmp:=celldata[k+1][cell[2]].BoundaryImage;
                BI:=tmp.ListIFace;
                SGN:=tmp.ListSign;
                LstEl:=tmp.ListElt;
                P:=StructuralCopy(stabgrp[k+1][i]);
                if IsPNormal(P,l) then
                    P:=Normalizer(P,Center(SylowSubgroup(P,l)));
                fi;
                RP:=ResolutionFiniteGroup(P,n);
                for r in [1..Size(reducedTorsionCells[k])] do
                    s:=reducedTorsionCells[k][r][2];

                    pos:=Positions(BI,s);
                    multiple[k][i][r]:=Sum(SGN{pos}) mod l;
#                    Print([k,i,r,s],BI,pos,multiple[k][i][r],"\n");
                    if not multiple[k][i][r]=0 then


                        Q:=StructuralCopy(stabgrp[k][r]);
                        if IsPNormal(Q,l) then
                            Q:=Normalizer(Q,Center(SylowSubgroup(Q,l)));
                        fi;
                        RQ:=ResolutionFiniteGroup(Q,n);
                        t:=Position(BI,s);
                        Pt:=ConjugateGroup(P,LstEl[t]);
                        for g in  stabgrp[k][r] do
                            if IsSubgroup(Q,ConjugateGroup(Pt,g)) then
                                break;
                            fi;
                        od;
                        map:=GroupHomomorphismByFunction(P,
                             Q,x->(LstEl[t]*g)^-1*x*(LstEl[t]*g));
                        inclusionMaps[k][i][r]:=LstEl[t];
                        eqmap:=EquivariantChainMap(RP,RQ,map);

                        T:=HomToIntegersModP(eqmap,l);
                        maps[k][i][r]:=T;

                    else

                        maps[k][i][r]:=0;
                    fi;
                od;
            od;
      od;

            stabres:=[];
      for k in [1..N-1] do

            stabres[k]:=[];
      od;
      for k in [1..N-1] do
            for j in [1..Size(maps[k][1])] do
              if not IsBound(stabres[k][j]) then
                i:=1;
                while i<=Size(maps[k]) do
                    if maps[k][i][j]=0 then
                        i:=i+1;
                    else break;
                    fi;
                od;
                if i>Size(maps[k]) then
                    P:=StructuralCopy(stabgrp[k][j]);
                    if IsPNormal(P,l) then
                    fi;
                    RP:=ResolutionFiniteGroup(P,n);
                    stabres[k][j]:=HomToIntegersModP(RP,l);
                fi;
              fi;
            od;



            for i in [1..Size(maps[k])] do
              if not IsBound(stabres[k][i]) then
                j:=1;
                while j<=Size(maps[k][i]) do
                    if maps[k][i][j]=0 then
                        j:=j+1;
                    else break;
                    fi;
                od;
                if j>Size(maps[k][i]) then
                    P:=StructuralCopy(stabgrp[k][i]);
                    if IsPNormal(P,l) then
                    fi;
                    RP:=ResolutionFiniteGroup(P,n);

                    stabres[k+1][i]:=HomToIntegersModP(RP,l);
                fi;
              fi;
            od;

       od;


###################### d1 differential#### ################
Differential:=function(k,p,q)
local w,i,j,A,B,CH,temp,x,M,Mat,BMat,t;
    if k=1 then
        N:=Length(reducedTorsionCells);
        if (p < 0) or (p > N-2) then return [];fi;

            CH:=[];
            Mat:=[];

            for i in [1..Size(reducedTorsionCells[p+2])] do
                CH[i]:=[];
                Mat[i]:=[];
                for j in [1..Size(reducedTorsionCells[p+1])] do
                    if maps[p+1][i][j]=0 then CH[i][j]:=0;
                    else
                        CH[i][j]:=Cohomology(maps[p+1][i][j],q);
                    fi;
                od;
            od;
            A:=[];

            for j in [1..Size(CH[1])] do
                i:=1;
                while i<=Size(CH) do

                    if CH[i][j]=0 then

                        i:=i+1;
                    else

                        A[j]:=Size(AbelianInvariants(Source(CH[i][j]))); break;
                    fi;
                od;
                if i>Size(CH) then
                    A[j]:=Cohomology(stabres[p+1][j],q);
                fi;
            od;

            B:=[];

            for i in [1..Size(CH)] do
                j:=1;

                while j<=Size(CH[1]) do
                    if CH[i][j]=0 then j:=j+1;
                    else B[i]:=Size(AbelianInvariants(Target(CH[i][j]))); break;
                    fi;
                od;
                if j>Size(CH[1]) then
                    B[i]:=Cohomology(stabres[p+2][i],q);
                fi;
            od;

            for i in [1..Size(reducedTorsionCells[p+2])] do
                Mat[i]:=[];
                for j in [1..Size(reducedTorsionCells[p+1])] do
                    if maps[p+1][i][j]=0 then
                        if A[j]*B[i]=0 then Mat[i][j]:=[];
                        else Mat[i][j]:=0*RandomMat(B[i],A[j]);fi;
                    else

                        Mat[i][j]:=multiple[p+1][i][j]*GroupHomomorphismToMatrix(CH[i][j],l);

                    fi;

                od;
            od;

            BMat:=[];
            for i in [1..Size(Mat)] do
                M:=[];
                for t in [1..Size(Mat[1])] do
                    t:=Size(Mat[i][1]);
                    if not t=0 then break;fi;
                od;
                for r in [1..t] do
                    M[r]:=[];
                    for j in [1..Size(Mat[1])] do
                        if not IsEmpty(Mat[i][j]) then
                            Append(M[r],Mat[i][j][r]);
                        fi;
                    od;
                od;
                Append(BMat,M);
            od;
#Print(1/0);
            return [Sum(A),Sum(B),BMat];
    else
    return fail;
    fi;
end;
######################End of d1 differential###############
######################### E_n Page ########################
EnPage:=function(k,p,q)
if k=1 then return E1page(p,q);
elif k=2 then return E2page(p,q);
else return fail;
fi;
end;
######################End of d1 differential###############


###########################################################
    return Objectify(HapEquivariantSpectralSequencePage,
                rec(
                    page:=EnPage,
                    differential:=Differential,
                    groupname:=groupname,
                    torsion:=l,
                    cohomology:=CohomologyOfGroup,
                    maps:=maps,
                    multiple:=multiple,
                    inclusionMaps:= inclusionMaps
                    ));

end);

###########################################################






InstallGlobalFunction("GroupHomomorphismToMatrix", function( phi,p)
local i,x,vectors,fgensA,fgensB,A,B,M,FreeGenerators, ElementToWord;


    ElementToWord:=function( G,L, g)
    local i,x,vectors,a;
        vectors:=CombinationDisjointSets(List([1..Size(L)],w->p));
        for x in vectors do
            a:=One(G);
            for i in [1..Size(x)] do
                a:=a*L[i]^x[i];
            od;
            if g=a then return x;fi;

        od;
    return false;
    end;
    A:=Source(phi);
    B:=Target(phi);

    FreeGenerators:=function(G)
    local gens,i,fgens;
        gens:=GeneratorsOfGroup(G);
        fgens:=[];
        for x in gens do
            if not Order(x)=1 then
                 if ElementToWord(G,fgens,x)=false then
                     Add(fgens,x);
                 fi;
            fi;
        od;
    return fgens;
    end;



    fgensA:=FreeGenerators(A);
    fgensB:=FreeGenerators(B);
    M:=[];
    for i in [1..Size(fgensA)] do
        Add(M,ElementToWord(B,fgensB,Image(phi,fgensA[i])));
    od;
    return TransposedMat(M);
end);

######################################################
DeclareGlobalFunction("MatrixSize");
InstallGlobalFunction(MatrixSize,
function(M)
return [Length(M),Length(M[1])];
end);
