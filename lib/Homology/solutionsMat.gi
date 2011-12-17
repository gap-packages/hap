#(C) Graham Ellis

#############################################################################
##
#M  SolutionsMatDestructive( <mat>, <matvec> ) . . . .  one solution for each
##							of several equations
##
##  One solution <x> of <x> * <mat> = <vec_i> or `fail' for each row vec_i
##  of <matvec>
##
InstallMethod( SolutionsMatDestructive,
        "generic method",
    [ IsOrdinaryMatrix and IsMutable,
      IsOrdinaryMatrix and IsMutable],
        function( mat, matvec )
    local i,ncols,sem, vno, z,x, row, sol,vec,solmat,z1,sol1,lngm;

    solmat:=[];
    z1 := Zero(mat[1][1]);
    sol1:=ListWithIdenticalEntries(Length(mat),z1);
    lngm:=Length(mat[1]);
    sem := SemiEchelonMatTransformationDestructive(StructuralCopy(mat));
    ncols := Length(matvec[1]);
    

for vec in matvec do
    z := StructuralCopy(z1);
    sol := StructuralCopy(sol1);
    ConvertToVectorRepNC(sol);
    if ncols <> lngm then
        Error("SolutionMat: matrix and vector incompatible");
    fi;
    for i in [1..ncols] do
        vno := sem.heads[i];
        if vno <> 0 then
            x := vec[i];
            if x <> z then
                AddRowVector(vec, sem.vectors[vno], -x);
                AddRowVector(sol, sem.coeffs[vno], x);
            fi;
        fi;
    od;
    if IsZero(vec) then
        Append(solmat,[sol]);
    else
        Append(solmat,[fail]);
    fi;
od;

return solmat;
end);

