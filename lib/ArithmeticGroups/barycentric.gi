
#######################################################################
#0
#F  ControlledSubdivision
##

InstallGlobalFunction(HAP_BaryCentricSubdivisionGComplex,
function(C)
local Cells,coBoundaries,N,i,Dims,pos,j,x,y,w,id,t,k,ck,c,s,a,v,g,b,
    Elts,Rep,mult,ListUnion, Chains, IsSameOrbit, AddReturn,
    Orbit, Dimension, StabRec, Action, Stabilizer, Boundary,
    NChains, BoundaryRec, FinalBoundary;

    
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
    IsSameOrbit:=function(a,b)
    local s,w,v,x,y,k,g,i;
        for i in [1..Length(a)] do
            if not a[i][1]=b[i][1] then
                return false; 
            fi;
        od;
        x:=List([1..Length(a)],i->[a[i][1],Cells[a[i][1]+1][a[i][2]]]);
        y:=List([1..Length(b)],i->[b[i][1],Cells[b[i][1]+1][b[i][2]]]);
        for i in [1..Length(x)] do
            if not x[i][2][1]=y[i][2][1] then
                return false; 
            fi;
        od;
        w:=List(C!.stabilizer(x[1][1],x[1][2][1]),i->Elts[y[1][2][2]]*i*Elts[x[1][2][2]]^-1);
        for s in [2..Length(x)] do
            v:=List(C!.stabilizer(x[s][1],x[s][2][1]),i->Elts[y[s][2][2]]*i*Elts[x[s][2][2]]^-1);
            w:=Intersection(w,v);
            if IsEmpty(w) then return false;fi;
        od;
        if not IsEmpty(w) then 
       
    return w[1];fi;
    end;
    ##################################################################
    Dims:=[];
    N:=Length(C);       #Added by Graham 02/2024
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

# Recognizing orbits and compute the boundary of cells.

    Orbit:=[];
    Orbit[N+1]:=[];
    Add(Orbit[N+1],NChains[1]);
    for i in [2..Length(NChains)] do
        for j in [1..Length(Orbit[N+1])] do
            c:=0;
            if not IsSameOrbit(NChains[i],Orbit[N+1][j])=false then
                c:=1;
                break;
            fi;
        od;
        if c=0 then Add(Orbit[N+1],NChains[i]);fi;
    od;

    BoundaryRec:=[];
    k:=N+1;
    while k>1 do
        BoundaryRec[k-1]:=[];
        for i in [1..Length(Orbit[k])] do
            x:=StructuralCopy(Orbit[k][i]);

            b:=[];
            for j in [1..Length(x)] do
                w:=StructuralCopy(x);
                Remove(w,j);

                if not IsBound(Orbit[k-1]) then
                    Orbit[k-1]:=[];
                fi;
                c:=0;
                for s in [1..Length(Orbit[k-1])] do
                    g:=IsSameOrbit(Orbit[k-1][s],w);
                 
                    if not g=false then
                        c:=1; 
                        Add(b,[(-1)^j*s,pos(g)]);
                    fi;
                od; 
                if c=0 then
                    Add(Orbit[k-1],w);
                    Add(b,[(-1)^j*Length(Orbit[k-1]),id]);
                fi;        
            od;
            BoundaryRec[k-1][i]:=b;
        od;
        k:=k-1;
    od;
            
# Find the k-rank
    Dimension:=function(k)
        if k<0 or k>N then return 0;fi;
    return Length(Orbit[k+1]);
    end;   

# Stabilizer subgroup of the representative of the ith-orbit of (k-1)-cells

    StabRec:=[];
    for k in [1..(N+1)] do
        StabRec[k]:=[];
        for i in [1..Dimension(k-1)] do
            a:=Orbit[k][i];
            x:=List([1..Length(a)],w->[a[w][1],Cells[a[w][1]+1][a[w][2]]]);
            w:=ConjugateGroup(C!.stabilizer(a[1][1],x[1][2][1]),Elts[x[1][2][2]]^-1);
            for s in [2..Length(x)] do
                v:=ConjugateGroup(C!.stabilizer(a[s][1],x[s][2][1]),Elts[x[s][2][2]]^-1);
                w:=Intersection(w,v);
            od;
            StabRec[k][i]:=w;
         od;
    od;

    Stabilizer:=function(k,i)
    return StabRec[k+1][i];
    end;

# The cell structure is rigid under the action of G then Action(k,i,j) always be 1.

    Action:=function(k,i,j)
    return 1;
    end;

# Calculate the boundary of the representative of the ith-orbit of k-cells
    
    Boundary:=function(n,k)
    if k>0 then 
        return BoundaryRec[n][k];
    else 
        return NegateWord(BoundaryRec[n][AbsInt(k)]);
    fi;
    end;
          
    ##################################################################
return Objectify(HapNonFreeResolution,
    rec(
    dimension:=Dimension,
    Orbit:=Orbit,
    Cells:=Cells,
    Chains:=Chains,
    boundary:=Boundary,
    coBoundaries:=coBoundaries,
    IsSameOrbit:=IsSameOrbit,
    homotopy:=fail,
    elts:=Elts,
    group:=C!.group,
    stabilizer:=Stabilizer,
    action:=Action,
    properties:=
    [["length",Maximum(1000,N)],
    ["characteristic",0],
    ["type","resolution"]]  ));
end);


################### end of ControlledSubdivision ############################

#########################################################
InstallMethod(BarycentricSubdivision,
"for non-free resolutions",
[IsHapNonFreeResolution],
function(R);
return HAP_BaryCentricSubdivisionGComplex(R);
end);
#########################################################



#########################################################
