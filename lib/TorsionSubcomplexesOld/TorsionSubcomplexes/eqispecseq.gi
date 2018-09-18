############################################################
##                                                        ##
##               equispecseq.gi                      ##
## HAP subpackage for GAP (Groups Algorithms Programming) ##
##         under the GNU GPL license (v. 3),  2012        ##
##      by Alexander D. Rahm & Bui Anh Tuan               ##
############################################################   	



InstallGlobalFunction("EquivariantSpectralSequencePage", function( C, m,n)
#########################################################################
                                             
#########################################################################
local reducedTorsionCells,celldata,j,N,l,P,Q,RP,RQ,g,Pt,
stabgrp,cell,i, E1page, T, EnPage, Differential, CohomologyOfGroup, stabres,stabcohom, inclusionMaps, groupname,name,sb,se,
maps,map,eqmap,tmp,BI,SGN,LstEl,s,r,t,sign,E2page,CH1,Mat1
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
   if N>2 then return fail;fi; 

else


#if not IsHapTorsionSubcomplex(C) then
    return fail;
#else
    N:=Size(C!.reducedTorsionCells);

    ## We only consider the case when length of the subcomplex is 
    ## less than 2 in order to get rid of the differential d2
    if N>2 then return fail;fi; 

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

###################### E_1 Page of Cohomology##############
E1page:=function(p,q)
local w,i;
    if p=0 then return Differential(1,0,q)[1];
    elif p=1 then return Differential(1,0,q)[2];
    else return fail;
    fi;
end;
######################End of E_1 Page######################
E2page:=function(p,q)
local t,M;
M:=Differential(1,0,q)[3];
if IsEmpty(M) then 
   if p=0 then 
       return 0;
   else return E1page(1,q);
   fi;
fi;

if p=0 then return Size(M[1])-RankMat(M);
elif p=1 then 
return E1page(1,q)-RankMat(M);
else return fail;
fi;

end;##############End of E_1 Page######################
###################### Cohomology of Group ################
CohomologyOfGroup:=function(k)
local p,w;
    w:=0;
    for p in [0..Minimum(k,1)] do
        w:=w+E2page(p,k-p);
    od;
    return w;
return fail;
end;
######################End of Cohomology####################

            inclusionMaps:=[];
            maps:=[];
            sign:=[];
            N:=2;
            for i in [1..Size(reducedTorsionCells[N])] do
                cell:= reducedTorsionCells[N][i];
                maps[i]:=[];
                inclusionMaps[i]:=[];
                sign[i]:=[];

                tmp:=celldata[N][cell[2]].BoundaryImage;
                BI:=tmp.ListIFace;
                SGN:=tmp.ListSign;
                LstEl:=tmp.ListElt;
                P:=StructuralCopy(stabgrp[N][i]);
                if IsPNormal(P,l) then 
                    P:=Normalizer(P,Center(SylowSubgroup(P,l)));
                fi;
                RP:=ResolutionFiniteGroup(P,n);
                for r in [1..Size(reducedTorsionCells[N-1])] do
                    s:=reducedTorsionCells[N-1][r][2];
                    sign[i][r]:=0;
                    if s in BI then
                        

                        Q:=StructuralCopy(stabgrp[N-1][r]);
                        if IsPNormal(Q,l) then 
                            Q:=Normalizer(Q,Center(SylowSubgroup(Q,l)));
                        fi;
                        RQ:=ResolutionFiniteGroup(Q,n);
                        t:=Position(BI,s);
                        Pt:=ConjugateGroup(P,LstEl[t]);
                        for g in  stabgrp[N-1][r] do                                   
                            if IsSubgroup(Q,ConjugateGroup(Pt,g)) then
                                break;
                            fi;
                        od;
                        map:=GroupHomomorphismByFunction(P,
                             Q,x->(LstEl[t]*g)^-1*x*(LstEl[t]*g));
                        inclusionMaps[i][r]:=LstEl[t];
                        eqmap:=EquivariantChainMap(RP,RQ,map);

                        T:=HomToIntegersModP(eqmap,l);
                        maps[i][r]:=T;
                        sign[i][r]:=SGN[t];
                    else
                        
                        maps[i][r]:=0;
                    fi;
                od;
            od;

            stabres:=[];
            for j in [1..Size(maps[1])] do
                i:=1;
                while i<=Size(maps) do
                    if maps[i][j]=0 then 
                        i:=i+1;
                    else break;
                    fi;
                od;
                if i>Size(maps) then
                    P:=StructuralCopy(stabgrp[N-1][j]);
                    if IsPNormal(P,l) then 
                        P:=Normalizer(P,Center(SylowSubgroup(P,l)));
                    fi; 
                    RP:=ResolutionFiniteGroup(P,n); 
                    stabres[j]:=HomToIntegersModP(RP,l);                 
                fi;                    
            od;

###################### d1 differential#### ################
Differential:=function(k,p,q)
local w,i,j,A,B,CH,temp,x,M,Mat,BMat,t;
    if k=1 then
        if not p=0 then return [];fi;
            CH:=[];
            Mat:=[];

            for i in [1..Size(reducedTorsionCells[N])] do
                CH[i]:=[];
                Mat[i]:=[];
                for j in [1..Size(reducedTorsionCells[N-1])] do
                    if maps[i][j]=0 then CH[i][j]:=0;
                    else
                        CH[i][j]:=Cohomology(maps[i][j],q);
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
                    A[j]:=Cohomology(stabres[j],q);
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
            od;

            for i in [1..Size(reducedTorsionCells[N])] do
                Mat[i]:=[];
                for j in [1..Size(reducedTorsionCells[N-1])] do
                    if maps[i][j]=0 then 
                        if A[j]*B[i]=0 then Mat[i][j]:=[];
                        else Mat[i][j]:=0*RandomMat(B[i],A[j]);fi;
                    else

                        Mat[i][j]:=sign[i][j]*GroupHomomorphismToMatrix(CH[i][j],l);

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
                    sign:=sign,
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