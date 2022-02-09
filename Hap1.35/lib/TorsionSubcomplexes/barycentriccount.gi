#######################################################################
#0
#F  ControlledSubdivision
##  Input: A pair of positive integers (m,n) 
##         
##  Output: The first n+1 terms of a free ZG-resolution  
##          where G is SL2Z(1/m)
##

InstallGlobalFunction(CountingCellsOfBaryCentricSubdivision,
function(C)
local Cells,coBoundaries,N,i,Dims,pos,j,x,y,w,id,t,k,ck,c,s,a,v,g,b,d,
    Elts,Rep,mult,ListUnion, Chains, IsSameOrbit, AddReturn,
    Orbit, Dimension, StabRec, Action, Stabilizer, Boundary,
    NChains, BoundaryRec, FinalBoundary, Partition;

    
    Elts:=C!.elts;
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
    # returns  a  "canonical"  representative  of  the  right  coset 
    # Elts[g]*Stab[i+1][j]
    Rep:=function(i,j,g)

    return pos(CanonicalRightCountableCosetElement
                            (C!.stabilizer(i,j),Elts[g]^-1)^-1);
    end;
    ##################################################################
    AddReturn:=function(a,g)
    local b;
        b:=StructuralCopy(a);    
        Add(b,g);
    return b;
    end;
    ##################################################################
    mult:=function(L,g)

    return List(L,a->[a[1],pos(Elts[g]*Elts[a[2]])]);

    end;
    ##################################################################
    ListUnion:=function(x,y)
    local a;
        for a in y do
            if not a in x then
                Add(x,a);
            fi;
        od;
    end;

    ##################################################################
    Dims:=[];
    for i in [0..Length(C)] do
        if C!.dimension(i)=0 then N:=i-1; break; fi;
        Dims[i+1]:=C!.dimension(i);
    od;
    Cells:=[];
    coBoundaries:=[];
    id:=pos(One(C!.group));
    for i in [1..N+1] do 
        Cells[i]:=[];
        coBoundaries[i]:=[];
    od;
    for j in [1..Dims[N+1]] do
        Add(Cells[N+1],[j,id]);
    od;

# Construct the list of cells and the corresponding coboundary of those cells
    i:=N;
    while i>0 do
        for k in [1..Length(Cells[i+1])] do
            x:=Cells[i+1][k];
            w:=StructuralCopy(C!.boundary(i,AbsInt(x[1])));
            w:=mult(w,x[2]);
            w:=List(w,a->[AbsInt(a[1]),Rep(i-1,AbsInt(a[1]),a[2])]);
            ListUnion(Cells[i],w);
            for y in w do
                t:=Position(Cells[i],y);
                if not IsBound(coBoundaries[i][t]) then
                    coBoundaries[i][t]:=[];
                fi;
                Add(coBoundaries[i][t],k);
            od;

        od;
        i:=i-1;
    od;        

# Record k-chains as a list 
    
    Chains:=[];

# Record the 1-chains 
    Chains[1]:=[];
    for i in [1..1] do

        for j in [1..Length(Cells[i])] do
            Add(Chains[1],[[i-1,j]]);
        od;
    od;
# Construct the list of N-chains
    for k in [1..(N)] do
        Chains[k+1]:=[];

            for i in [1..Length(Chains[k])] do
                ck:=StructuralCopy(Chains[k][i]);
                c:=ck[k];
                w:=List(coBoundaries[c[1]+1][c[2]],x->AddReturn(ck,[c[1]+1,x]));
                Append(Chains[k+1],w);
            od;    
        
    od;
    NChains:=StructuralCopy(Chains[N+1]);

    Chains:=[];
    Chains[N+1]:=StructuralCopy(NChains);
    k:=N+1;
    while k>1 do
 
        for i in [1..Length(Chains[k])] do
            x:=StructuralCopy(Chains[k][i]);

            b:=[];
            for j in [1..Length(x)] do
                w:=StructuralCopy(x);
                Remove(w,j);

                if not IsBound(Chains[k-1]) then
                    Chains[k-1]:=[];
                fi;
                if not w in Chains[k-1] then
                    Add(Chains[k-1],w);
                fi;
            od;

        od;
        k:=k-1;
    od;

    return List([1..Length(Chains)],i->Length(Chains[i]));
          
    ##################################################################
#return Objectify(HapNonFreeResolution,
#    rec(
#    Cells:=Cells,
#    Chains:=Chains,
#    coBoundaries:=coBoundaries,
#    homotopy:=fail,
#    properties:=
#    [["length",Maximum(1000,N)],
#    ["characteristic",0],
#    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################

