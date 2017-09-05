
#######################################################################
#0
#F  RigidFacetsSubdivision
##  Input: A pair of positive integers (m,n) 
##         
##  Output: The first n+1 terms of a free ZG-resolution  
##          where G is SL2Z(1/m)
##

InstallGlobalFunction(RigidFacetsSubdivision,
function(arg)
local W, StabRec, i, j, N, M, x, bdry, s1, s2, p, k, w, t,
      DimRec, BoundaryRec, id, dims, NotRigid, NewCell, IsAdjacent,
      Cell, Elts, Boundary, Dimension, CLeftCosetElt,
      pos, Stab, Mult,DimTemp, BoundaryTemp, Partition,
      Stabilizer, Action, IsRigidCell, ReplaceCell, SubdividingCell,
      Orbit, C;

    C:=arg[1];
    i:=0;
    while C!.dimension(i)>0 do
        i:=i+1;
    od;
    N:=i-1; # Length of the chain complex   
    M:=N;

if Length(arg)=2 then M:=Minimum(N,arg[2]);
fi;
    

    Elts:=C!.elts;
    StabRec:=[];
    DimRec:=[];

    ##################################################################
    # If g in Elts return the position of g in the list,
    # otherwise, add g to Elts and return the position.
    pos:=function(g)
    local posit;

        posit:=Position(Elts,g);
        if posit=fail then 
            Add(Elts,g);  
            return Length(Elts);
        else  
            return posit;
        fi;
    end;
    ##################################################################
    id:=pos(One(C!.group));
    ##################################################################
    # return the stabilizer of g*e,
    # 
    Stab:=function(e,g)
    return ConjugateGroup(StabRec[e[1]+1][AbsInt(e[2])],Elts[g]^-1);
    end;
    ##################################################################
    # returns  a  "canonical"  representative  of  the  right  coset 
    # Elts[g]*Stab[i+1][j]
    CLeftCosetElt:=function(i,j,g)

    return pos(CanonicalRightCountableCosetElement
                            (StabRec[AbsInt(i)+1][j],Elts[g]^-1)^-1);
    end;
    ##################################################################
    ##
    ##  Input:  A list L, degree k, position g of an element    
    ##  Output: Product of g and L.
    ##
    Mult:=function(L,k,g)
    local x,w,t,h,y,vv,LL;
        vv:=[];
        LL:=ShallowCopy(L);
        for x in [1..Length(LL)] do
            w:=Elts[g]*Elts[LL[x][2]];
            t:=CLeftCosetElt(k,AbsInt(LL[x][1]),pos(w));
            Add(vv,[LL[x][1],t]);
        od;
        return vv;
    end;
    ###################################################################
    # Store essential data: stabilizers, boundaries, dimensions


    NewCell:=[];
    for i in [1..N] do
        NewCell[i]:=[];
    od;
    
    Orbit:=[];
    for i in [1..N+1] do
        Orbit[i]:=[];
    od;

    for i in [0..N] do
        StabRec[i+1]:=[];
        DimRec[i+1]:=C!.dimension(i);
        for j in [1..C!.dimension(i)] do
            StabRec[i+1][j]:=C!.stabilizer(i,j);
        od;
    od;
    
    BoundaryRec:=[];
    for i in [1..N] do
        BoundaryRec[i]:=[];
        for j in [1..DimRec[i+1]] do
            bdry:=C!.boundary(i,j);
            BoundaryRec[i][j]:=[];
            for x in bdry do
                s1:=C!.action(i-1,AbsInt(x[1]),x[2]);
                p:=pos(CanonicalRightCountableCosetElement
                            (C!.stabilizer(i-1,AbsInt(x[1])),Elts[x[2]]^-1)^-1);
                s2:=C!.action(i-1,AbsInt(x[1]),p);
       
                Add(BoundaryRec[i][j],[s1*s2*x[1],p]);
            od;

        od;
    od;
    ##################################################################
 
    # Data type for a k-cell with stabilizer stab and boundary bdry
    Cell:=function(k,stab,bdry)
    return rec(
        dimension:=k,
        stabilizer:=stab,
        boundary:=bdry
    );
    end;
  
    ##################################################################
    # Check if the cell is whether rigid or not
 
    IsRigidCell:=function(k,m)
    local bdry, intst, L;
        bdry:=BoundaryRec[k][m];
        L:=List(bdry,w->Elements(ConjugateGroup(StabRec[k][AbsInt(w[1])],Elts[w[2]]^-1)));
        intst:=Intersection(L);
        if not Elements(StabRec[k+1][m])=Elements(intst) then
            return false;
        else return true;
        fi;

    end;
    ##################################################################
    # Check if 2 cells are adjacent
    IsAdjacent:=function(k,e,f,w)
    local bdry1,bdry2,i,l1,l2,s1,s2,s,a,s0,x;
        

        bdry1:=Mult(BoundaryRec[k][AbsInt(e[1])],k-1,e[2]);

        bdry2:=Mult(BoundaryRec[k][AbsInt(f[1])],k-1,f[2]);
#Print("\n [bdry1,bdry2] = ",[e,f,bdry1,bdry2],"\n");
        l1:=List(bdry1,w->[AbsInt(w[1]),w[2]]);
        l2:=List(bdry2,w->[AbsInt(w[1]),w[2]]);

        a:=Intersection(Elements(l1),Elements(l2));

        if IsEmpty(a) then return false;fi;
        return true; 
#        s:=Elements(Stab([k+1,w[1]],w[2]));
#        s1:=Elements(Stab([k,e[1]],e[2]));
#        s2:=Elements(Stab([k,f[1]],f[2]));
#        s0:=Elements(Stab([k-1,a[1][1]],a[1][2]));
#Print("\n [s0,s1,s2]=",[s0=s1,s1=s2],"\n");
#        if s1=s2 and s2=s0 then
#            return true;
#        else
#            return false;
#        fi;        
    end;

    ##################################################################
    # Subdividing a cell 
    SubdividingCell:=function(k,i)
    local bdry, w, x, d, y, a, b, sub, j, j1,j2, z, t, c, Flag, s, s1, ConnectToCenter,
          SearchComponent,temp, l1, bdry1, bdry2, Orb, IsSameOrbit, OrbFlag, OrbElm;

        DimRec[1]:=DimRec[1]+1;
        Add(StabRec[1],StabRec[k+1][i]);

        bdry:=ShallowCopy(BoundaryRec[k][i]);


        #############################################################
        # Components
        SearchComponent:=function(q)
        local x,z,j,w;
            x:=[];
            w:=[];
            Add(w,q);
            for j in [1..Length(Orb)] do
                if OrbFlag[j]=0 then
                    for j1 in [1..Length(Orb[j])] do
                        if Flag[j][j1]=0 and IsAdjacent(k-1,q,Orb[j][j1],[i,id]) then
                            Add(x,Orb[j][j1]);
                            Flag[j][j1]:=1;
                            OrbFlag[j]:=1;
                            break;
                        fi;
                    od;
                fi;
            od;
            for j in x do
                Append(w,SearchComponent(j));
            od;
            return w;
        end;
        #############################################################
    #################################################################
    # Connect the cell e to the barycenter of the cell f
    # e and f are in the form [k,i,g]: dimension k, obtain by sending 
    # ith-representative under the action of the element g in G 
    ConnectToCenter:=function(e)
    local bdry1, x, stab, bdrye, w, stablst, redbdry, y, a,
          LCoset, AddCell;


    ##################################################################
    # returns  a  "canonical"  representative  of  the  right  coset 
    # Elts[g]*Stab[i+1][j]
#    LCoset:=function(i,j,g)

#    return pos(CanonicalRightCountableCosetElement
#                            (StabTemp[AbsInt(i)+1][j],Elts[g]^-1)^-1);
#    end;
    ##################################################################
    # Add a k-cell with stabilizer stab and boundary bdry
    # to the cell complex
    AddCell:=function(m,stab,bdry,e)
    local j,g,s,w,u,v;

        if m<k then
            w:=e[1];
        for j in [1..DimTemp[m+1]] do
            u:=BoundaryTemp[m+1][j];
#Print("[w,u]=",[w,u],"\n");
            s:=DimRec[m+1]-DimTemp[m+1]+j;
            if AbsInt(w[1])=AbsInt(u[1]) then
#            v:=StructuralCopy(StabRec[m+1][s]);
            v:=StructuralCopy(StabRec[m][AbsInt(u[1])]);
            v:=List(v,a->Elts[w[2]]*a*Elts[u[2]]^-1);
#Print("\n v=",v,"\n");

            g:=Intersection(v,Elements(StabRec[k+1][i]));
#Print("\n [m,s,g]",[m,s,g],"\n");
            if not IsEmpty(g) then
#Print("\n [m,s,g]",[m,s,g[1]],"\n");
#Print("\n CLeftCoset", CLeftCosetElt(m,s,pos(g[1])),"\n");
                return [s, CLeftCosetElt(m,s,pos(g[1]))];
            fi;
            fi;
        od;
        
        DimTemp[m+1]:=DimTemp[m+1]+1;

        Add(BoundaryTemp[m+1],e[1]);
        DimRec[m+1]:=DimRec[m+1]+1;
        Add(StabRec[m+1],stab);
        Add(BoundaryRec[m],bdry);
        else 
            DimRec[m+1]:=DimRec[m+1]+1;
            Add(StabRec[m+1],stab);
            Add(BoundaryRec[m],bdry);
        fi;
#Print("\n [m,DimRec[m+1],id]",[m,DimRec[m+1],id],"\n");
#Print("\n CLeftCoset", CLeftCosetElt(m,DimRec[m+1],id),"\n");
        return [DimRec[m+1],CLeftCosetElt(m,DimRec[m+1],id)];    
    end;    
    ################################################################## 
 
        if e[1]=0 then 

            bdry1:=[[-SignInt(e[2][1][1])*DimRec[1],CLeftCosetElt(0,DimRec[1],id)],[e[2][1][1],e[2][1][2]]];

#            stab:=Intersection(Elements(Stab([0,e[2][1][1]],e[2][1][2])),Elements(StabRec[k+1][i]));
            stab:=Intersection(Stab([0,e[2][1][1]],e[2][1][2]),StabRec[k+1][i]);

            return AddCell(1,stab,bdry1,e[2]);
        fi;
#Print("test");        
        stablst:=[];
        bdry1:=[];


        Append(bdry1,e[2]);
        bdrye:=[];
        for y in e[2] do 
            a:=Mult(BoundaryRec[e[1]][AbsInt(y[1])],e[1]-1,y[2]);
            if y[1]<0 then a:=NegateWord(a);fi;
            Append(bdrye,a);
            Add(stablst,Stab([e[1],y[1]],y[2]));
        od;
        Add(stablst,StabRec[k+1][i]);
        bdrye:=AlgebraicReduction(bdrye);
#if e[1]=(k-1) then
#Print("\n Length of bdrye=",Length(bdrye),"\n");
#fi;
#Print("bdrye=",bdrye,"\n");
        for x in bdrye do
            w:=ConnectToCenter([e[1]-1,[[AbsInt(x[1]),x[2]]]]);
            Add(bdry1,[-SignInt(x[1])*w[1],w[2]]);
#            Add(stablst,Stab([e[1],w[1]],w[2]));
        od;

        stab:=Intersection(stablst);
#Print("\n stab=",stab,"\n");
        return AddCell(e[1]+1,stab,bdry1,e[2]);
    end;
    ##################################################################
    IsSameOrbit:=function(e,f)
    local s1,s;
        if AbsInt(e[1])=AbsInt(f[1]) then
            s:=StabRec[k+1][i];
            s1:=StabRec[k][AbsInt(e[1])];
            for x in s do
                if Elts[f[2]]^-1*x*Elts[e[2]] in s1 then
                    return SignInt(e[1])*SignInt(f[1])*pos(x);
                fi;
            od;
            return false;        
         else
            return false;
         fi;
    end;
    ###################################################################
# Group the boundary cells into orbits under the action of the stabilizer
        Orb:=[];
        Orb[1]:=[];
        OrbElm:=[];
        OrbElm[1]:=[];
        Add(Orb[1],bdry[1]);
        Add(OrbElm[1],id);
#        s:=1;

        for j in [2..Length(bdry)] do
            t:=0;
            for j1 in [1..Length(Orb)] do
                s:=IsSameOrbit(bdry[j],Orb[j1][1]);
                if not (s=false) then
                    Add(Orb[j1],bdry[j]);
                    Add(OrbElm[j1],s);
                    break;
                fi;
                t:=t+1;
            od;    
            if t=Length(Orb) then Add(Orb,[bdry[j]]);Add(OrbElm,[id]);fi;
        od;
#Print("\n bdry = ",bdry,"\n");
#Print("\n Orb=:",Orb,"\n");
     ##################################################################
# Divide the boundary into rigid parts in the big cell [k,i]
       

        Flag:=[];
        for j in [1..Length(Orb)] do
            Flag[j]:=[];
            for j1 in [1..Length(Orb[j])] do
                Flag[j][j1]:=0;
            od;
        od;
        sub:=[];
        for j in [1..Length(Orb[1])] do
            sub[j]:=[Orb[1][j]]; 
            Flag[1][j]:=1;       
        od; 
        

#        for j in [1..Length(sub)] do
           OrbFlag:=[];
           for j1 in [1..Length(Orb)] do
               OrbFlag[j1]:=0;
           od;
           OrbFlag[1]:=1;
#           sub[j]:=SearchComponent(sub[j][1]);
           sub[1]:=SearchComponent(sub[1][1]);
#        od;

     ###################################################################
        w:=[];
        BoundaryTemp:=[];
        DimTemp:=[];
        for j in [0..(k)] do
            BoundaryTemp[j+1]:=[];
            DimTemp[j+1]:=0;
        od;
#Print("\n Flag ",Flag,"\n");
#Print("sub=",sub);
        d:=ConnectToCenter([k-1,sub[1]]);
#Print("\n connecttocenter=",d,"\n");
#Print("dim1=",DimTemp[2]);
        Add(w,d);
        for j in [2..Length(OrbElm[1])] do
            Add(w,[SignInt(OrbElm[1][j])*d[1],pos(Elts[AbsInt(OrbElm[1][j])]*Elts[d[2]])]);
        od;
#        for j in [2..Length(sub)] do
#            bdry1:=ShallowCopy(sub[1]);
#            bdry2:=ShallowCopy(sub[j]);
#            for a in Elements(StabRec[k+1][i]) do
#                t:=pos(a);
#                l1:=Mult(bdry1,k-1,t);
#
#                if Set(l1)=Set(bdry2) then
#                    Add(w,[d[1],t]);
#                    break;
#                fi;
#                if Set(l1)=Set(NegateWord(bdry2)) then
#                   Add(w,[-d[1],t]);
#                   break;
#                fi;
#            od;
#
#        od;
#Print("\n subdividing cell=",w,"\n");
        return w;
    end;
    ##################################################################
    # Replacing a cell by its subdivision
    ReplaceCell:=function(k,m)
    local i, j, p, w, x, bdry, y, ww;
        w:=ShallowCopy(SubdividingCell(k,m));
        if k=N then 
            Partition:=StructuralCopy(w);
            for i in [1..Length(Partition)] do
                Partition[i][1]:=AbsInt(Partition[i][1]-SignInt(Partition[i][1]));
            od;
        fi;

#        if k<N then
        if k<=M and k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=ShallowCopy(BoundaryRec[k+1][i]);
            p:=PositionsProperty(bdry,w->AbsInt(w[1])=m);
            for j in p do
                x:=bdry[j];
                ww:=ShallowCopy(w);
                if x[1]<0 then ww:=NegateWord(ww);fi;
                ww:=Mult(ww,k,x[2]);
#Print("\n ww=",ww,"\n");
                Append(bdry,ww);
            od;
            y:=bdry{p};
            bdry:=Set(bdry);
            SubtractSet(bdry,y);
            BoundaryRec[k+1][i]:=bdry;
#Print("\n bdry=",bdry,"\n");
        od;
        fi;
        BoundaryRec[k][m]:="del";
        StabRec[k+1][m]:="del";
    end;
    ##################################################################
    # Main part: subdividing the fundamental domain
    NotRigid:=[];
    dims:=ShallowCopy(DimRec);
    i:=1;
#    Print("The cells which are not rigid: \n");
    while i<=M do
#    while i<=N do
        j:=1;
        while j<=dims[i+1] do
            if not IsRigidCell(i,j) then
#                Print([i,j]);
                Add(NotRigid,[i,j]);
            fi;
            j:=j+1;
        od;
        i:=i+1;
    od;
    for x in NotRigid do
#        Print("\n The cell ",x," is in process of subdividing \n");
        ReplaceCell(x[1],x[2]);
    od;

    
    #Delete cells which are already replaced by its subdivision
#    Print("Deleting cells which are already replaced by its subdivision... \n");
    t:=1;
    for w in [1..Length(NotRigid)] do
        k:=NotRigid[w][1];
        j:=NotRigid[w][2];
        if k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=BoundaryRec[k+1][i];

            if not IsString(bdry) then
               for x in bdry do
                    if AbsInt(x[1])>j then 
                        x[1]:=x[1]-SignInt(x[1]);
                    fi;
                od;
             fi;
             BoundaryRec[k+1][i]:=bdry;

         od;
         fi;
         dims[k+1]:=dims[k+1]-1;
         DimRec[k+1]:=DimRec[k+1]-1;
         Remove(BoundaryRec[k],j);
         Remove(StabRec[k+1],j);
         if IsBound(NotRigid[w+1]) and NotRigid[w+1][1]=NotRigid[w][1] then
             NotRigid[w+1][2]:=NotRigid[w+1][2]-t;
             t:=t+1;
         else
             t:=1;
         fi;  

    od;
#    Print("Done!","\n");
    ##################################################################
    Boundary:=function(k,m)
        if m<0 then return NegateWord(BoundaryRec[k][AbsInt(m)]);fi;
        return BoundaryRec[k][m];
    end;

    Stabilizer:=function(k,m)
        return StabRec[k+1][m];
    end;

    Dimension:=function(k)
        if k>N then return 0;fi;
        return DimRec[k+1];
    end;

    Action:=function(k,i,j)
        return 1;
    end;
    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    Partition:=Partition,
    boundary:=Boundary,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    subdividing:=SubdividingCell,
    replacecell:=ReplaceCell,
#    issameorbit:=IsSameOrbit,
    isrigid:=IsRigidCell,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################

