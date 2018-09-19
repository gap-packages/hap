
#######################################################################
#0
#F  ControlledSubdivision
##  Input: A pair of positive integers (m,n) 
##         
##  Output: The first n+1 terms of a free ZG-resolution  
##          where G is SL2Z(1/m)
##

InstallGlobalFunction(BaryCentricSubdivision,
function(C)
local W, StabRec, i, j, N, x, bdry, s1, s2, p, k, w, t,
      DimRec, BoundaryRec, id, dims, NotRigid, NewCell,
      AddCell, Cell, Elts, Boundary, Dimension, CLeftCosetElt,
      pos, IsSameOrbit, Stab, Mult, ConnectToCenter,
      Stabilizer, Action, IsRigidCell, ReplaceCell, SubdividingCell;
    
    

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
    return ConjugateGroup(StabRec[e[1]+1][e[2]],Elts[g]^-1);
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
    local x,w,t,h,y,vv;
        vv:=[];
        for x in [1..Length(L)] do
            w:=Elts[g]*Elts[L[x][2]];
            t:=CLeftCosetElt(k,AbsInt(L[x][1]),pos(w));
            Add(vv,[L[x][1],t]);
        od;
        return vv;
    end;
    ###################################################################
    # Store essential data: stabilizers, boundaries, dimensions

    i:=0;
    while C!.dimension(i)>0 do
        i:=i+1;
    od;
    N:=i-1; # Length of the chain complex
    NewCell:=[];
    for i in [1..N] do
        NewCell[i]:=[];
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
#            BoundaryRec[i][j]:=ShallowCopy(C!.boundary(i,j));
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
    # Add a k-cell with stabilizer stab and boundary bdry
    # to the cell complex
    AddCell:=function(k,stab,bdry)
    local i,g;
        if k=0 then 
            DimRec[k+1]:=DimRec[k+1]+1;
            Add(StabRec[k+1],stab);
            return [DimRec[k+1],CLeftCosetElt(0,DimRec[k+1],id)];
        fi;

        for i in [(dims[k+1]+1)..DimRec[k+1]] do
            g:=IsSameOrbit([k,StabRec[k+1][i],
                               BoundaryRec[k][i]],[k,stab,bdry]);
            if not g=false then
#Print("the cell ",[i, CLeftCosetElt(k,i,g)],"\n");
                return [i, CLeftCosetElt(k,i,g)];
            fi;
        od;
        DimRec[k+1]:=DimRec[k+1]+1;
        Add(StabRec[k+1],stab);
        Add(BoundaryRec[k],bdry);
        return [DimRec[k+1],CLeftCosetElt(k,DimRec[k+1],id)];    
    end;    
    ##################################################################
    # check if two k-cells are in the same orbit
    IsSameOrbit:=function(e,f)
    local p, bdry1, bdry2, i, a, b, x;
        if not e[1]=f[1] then
            return false;
        fi;
        bdry1:=ShallowCopy(e[3]);
        bdry2:=ShallowCopy(f[3]);
        bdry2:=List(bdry2,w->[w[1],CLeftCosetElt(e[1]-1,AbsInt(w[1]),w[2])]);
#Print("bdry1 ",bdry1,"\n");
#Print("bdry2 ",bdry2,"\n");
        p:=PositionsProperty(bdry2,w->AbsInt(w[1])=AbsInt(bdry1[1][1]));
#Print("p ",p,"\n");
        for i in p do
            for a in Elements(StabRec[e[1]][AbsInt(bdry1[1][1])]) do
                b:=Elts[bdry2[i][2]]*a*Elts[bdry1[1][2]]^-1;
                x:=List(bdry1,w->[w[1],CLeftCosetElt(e[1]-1,
                   AbsInt(w[1]),pos(b*Elts[w[2]]))]);
                if Set(x)=Set(bdry2) then 
#Print("b ",pos(b),"\n");
                    return pos(b);
                fi;
            od;    
        od;
        return false;
    end;
    ##################################################################
    # Connect the cell e to the barycenter of the cell f
    # e and f are in the form [k,i,g]: dimension k, obtain by sending 
    # ith-representative under the action of the element g in G 
    ConnectToCenter:=function(e,f)
    local bdry, x, stab, bdrye, w, stablst;
      
        if e[1]=0 then 
            stab:=Intersection(Stab([e[1],e[2]],e[3]),Stab([f[1],f[2]],f[3]));
            bdry:=[[-f[2],f[3]],[e[2],e[3]]];
#Print(e,"  ",AddCell(e[1]+1,stab,bdry),"\n");
            return AddCell(1,stab,bdry);
        fi;
        stablst:=[];
        Add(stablst,Stab([e[1],e[2]],e[3]));
#        stab:=Intersection(Stab([e[1],e[2]],e[3]),Stab([f[1],f[2]],f[3]));
        bdry:=[];
        Add(bdry,[e[2],e[3]]);
        bdrye:=Mult(BoundaryRec[e[1]][e[2]],e[1]-1,e[3]);
        for x in bdrye do
            w:=ConnectToCenter([e[1]-1,AbsInt(x[1]),x[2]],f);
            Add(bdry,[-SignInt(x[1])*w[1],w[2]]);
            Add(stablst,Stab([e[1],w[1]],w[2]));
        od;
        stab:=Intersection(stablst);
#Print(e,"  ",AddCell(e[1]+1,stab,bdry),"\n");
        return AddCell(e[1]+1,stab,bdry);
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
    # Subdividing a cell 
    SubdividingCell:=function(k,i)
    local bdry, w, x, d, y;
        y:=AddCell(0,StabRec[k+1][i],[]);
        bdry:=BoundaryRec[k][i];
        w:=[];
#Print([k,i],"  ",bdry,"\n");
        for x in bdry do
            d:=ConnectToCenter([k-1,AbsInt(x[1]),x[2]],[0,y[1],y[2]]);
            if x[1]<0 then 
                Add(w,[-d[1],d[2]]);
            else
                Add(w,d);
            fi;
        od;
        return w;
    end;
    ##################################################################
    # Replacing a cell by its subdivision
    ReplaceCell:=function(k,m)
    local i, j, p, w, x, bdry, y, ww;
        w:=ShallowCopy(SubdividingCell(k,m));
        if k<N then
        for i in [1..DimRec[k+2]] do
            bdry:=ShallowCopy(BoundaryRec[k+1][i]);
            p:=PositionsProperty(bdry,w->AbsInt(w[1])=m);
            for j in p do
                x:=bdry[j];
                ww:=ShallowCopy(w);
                if x[1]<0 then ww:=NegateWord(ww);fi;
                ww:=Mult(ww,k,x[2]);
                Append(bdry,ww);
            od;
            y:=bdry{p};
            bdry:=Set(bdry);
            SubtractSet(bdry,y);
            BoundaryRec[k+1][i]:=bdry;
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
    while i<=N do
        j:=1;
        while j<=dims[i+1] do
#            if not IsRigidCell(i,j) then
#                Print([i,j]);
                Add(NotRigid,[i,j]);
#            fi;
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
    boundary:=Boundary,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    subdividing:=SubdividingCell,
    replacecell:=ReplaceCell,
    issameorbit:=IsSameOrbit,
    isrigid:=IsRigidCell,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################

