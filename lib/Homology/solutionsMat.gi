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
    local i,ncols,sem, vno, z,x, row, sol,vec,solmat;

    sem := SemiEchelonMatTransformationDestructive(mat);
    solmat:=[];

for vec in matvec do
    ncols := Length(vec);
    z := Zero(mat[1][1]);
    sol := ListWithIdenticalEntries(Length(mat),z);
    ConvertToVectorRepNC(sol);
    if ncols <> Length(mat[1]) then
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

