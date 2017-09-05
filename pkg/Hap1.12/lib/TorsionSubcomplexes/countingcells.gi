
#######################################################################
#0
#F  Counting cells in a cell complex formatted in the Gcomplex datatype
##  Input: A cell complex formatted in Gcomplex datatype 
##         
##  Output: The structure in which we can extract the list of k-cells  
##          or the number of k-cells
##

InstallGlobalFunction(CountingCellsOfACellComplex,
function(C)
local Cells,N,i,Dims,pos,j,x,y,w,id,t,k,ck,c,s,a,v,g,b,
    Elts,Rep,mult,ListUnion, AddReturn,NCells,
    Orbit, nrCells,R, PWsbgrp;

    
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
#    AddReturn:=function(a,g)
#    local b;
#        b:=StructuralCopy(a);    
#        Add(b,g);
#    return b;
#    end;
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
    id:=pos(One(C!.group));
    for i in [1..N+1] do 
        Cells[i]:=[];

    od;
    
    if IsBound(C!.Partition) then
        Cells[N+1]:=StructuralCopy(C!.Partition);
    else
        for j in [1..Dims[N+1]] do
            Add(Cells[N+1],[j,id]);
        od;
    fi;




# Construct the list of cells and the corresponding coboundary of those cells
    i:=N;
    while i>0 do
        for k in [1..Length(Cells[i+1])] do
            x:=Cells[i+1][k];
            w:=StructuralCopy(C!.boundary(i,AbsInt(x[1])));
            w:=mult(w,x[2]);
            w:=List(w,a->[AbsInt(a[1]),Rep(i-1,AbsInt(a[1]),a[2])]);
            ListUnion(Cells[i],w);


        od;
        i:=i-1;
    od;        

    ##################################################################
    nrCells:=List([1..N+1],i->Length(Cells[i]));          
    ##################################################################
 

return nrCells;
end);


################### end of ControlledSubdivision ############################

InstallGlobalFunction(CountingControlledSubdividedCells,
function(C)
local Cells,N,i,Dims,pos,j,x,y,w,id,t,k,ck,c,s,a,v,g,b,L,intst,
    Elts,Rep,mult,ListUnion, AddReturn,NCells,d,PWsbgrp,bdry,
    Orbit, nrCells,R;

    
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
#    AddReturn:=function(a,g)
#    local b;
#        b:=StructuralCopy(a);    
#        Add(b,g);
#    return b;
#    end;
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
    id:=pos(One(C!.group));
    for i in [1..N+1] do 
        Cells[i]:=[];

    od;
    
    if IsBound(C!.Partition) then
        Cells[N+1]:=StructuralCopy(C!.Partition);
    else
        for j in [1..Dims[N+1]] do
            Add(Cells[N+1],[j,id]);
        od;
    fi;




# Construct the list of cells and the corresponding coboundary of those cells
    i:=N;
    while i>0 do
        for k in [1..Length(Cells[i+1])] do
            x:=Cells[i+1][k];
            w:=StructuralCopy(C!.boundary(i,AbsInt(x[1])));
            w:=mult(w,x[2]);
            w:=List(w,a->[AbsInt(a[1]),Rep(i-1,AbsInt(a[1]),a[2])]);
            ListUnion(Cells[i],w);


        od;
        i:=i-1;
    od; 
    ##################################################################
    PWsbgrp:=[];
    for i in [1..N] do
        PWsbgrp[i]:=[];
        for j in [1..Dims[i+1]] do
            bdry:=C!.boundary(i,j);
            L:=List(bdry,w->Elements(ConjugateGroup(C!.stabilizer(i-1,AbsInt(w[1])),Elts[w[2]]^-1)));
            intst:=Intersection(L);
            Add(PWsbgrp[i],Size(C!.stabilizer(i,j))/Size(intst));
        od;
    od;
 
    d:=[];
    for i in [1..N] do
        d[i]:=0;
        for x in Cells[i+1] do
            d[i]:=d[i]+PWsbgrp[i][AbsInt(x[1])];
        od;
        
    od;
    return d;
end);
#######################################################################################

InstallGlobalFunction(CountingBaryCentricSubdividedCells,
function(C)
local Cells,N,i,Dims,pos,j,x,y,w,id,t,k,ck,c,s,a,v,g,b,L,intst,
    Elts,Rep,mult,ListUnion, AddReturn,NCells,d,nrCells,bdry,
    Orbit, R;

    
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
#    AddReturn:=function(a,g)
#    local b;
#        b:=StructuralCopy(a);    
#        Add(b,g);
#    return b;
#    end;
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
    id:=pos(One(C!.group));
    for i in [1..N+1] do 
        Cells[i]:=[];

    od;
    
    if IsBound(C!.Partition) then
        Cells[N+1]:=StructuralCopy(C!.Partition);
    else
        for j in [1..Dims[N+1]] do
            Add(Cells[N+1],[j,id]);
        od;
    fi;




# Construct the list of cells and the corresponding coboundary of those cells
    i:=N;
    while i>0 do
        for k in [1..Length(Cells[i+1])] do
            x:=Cells[i+1][k];
            w:=StructuralCopy(C!.boundary(i,AbsInt(x[1])));
            w:=mult(w,x[2]);
            w:=List(w,a->[AbsInt(a[1]),Rep(i-1,AbsInt(a[1]),a[2])]);
            ListUnion(Cells[i],w);


        od;
        i:=i-1;
    od; 
    ##################################################################
    nrCells:=[];
    for i in [1..N] do
        nrCells[i]:=[];
        for j in [1..Dims[i+1]] do
            if i=1 then nrCells[i][j]:=2;
            else
                nrCells[i][j]:=0;
                bdry:=C!.boundary(i,j);
                for x in bdry do
                   nrCells[i][j]:=nrCells[i][j]+nrCells[i-1][AbsInt(x[1])];
                od;
            fi;
        od;
    od;
 
    d:=[];
    for i in [1..N] do
        d[i]:=0;
        for x in Cells[i+1] do
            d[i]:=d[i]+nrCells[i][AbsInt(x[1])];
        od;
        
    od;
    return d;
end);

