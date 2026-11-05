# A star has a unique center, so its setwise stabilizer is the center's
# stabilizer. Check balance using subgroup generators and cell cosets;
# L!.elts is only a table of elements used in boundaries, not the whole group.
InstallGlobalFunction(PreIsGBalanced, function(L,e,n,k)
local U, StabU, d, orbit, cells, first, stab, h, cell, g;

U:=kStar(L,e,n,k);
StabU:=L!.stabilizer(n,e);
for d in [2..k+1] do
    for orbit in Set(List(U[d],cell->cell[1])) do
        cells:=Filtered(U[d],cell->cell[1]=orbit);
        first:=cells[1][2];
        stab:=cells[1][3];

        # An element carrying h*f to l*f lies in l*Stab(f)*h^-1.
        # Every such element must fix the center, otherwise the translated
        # star overlaps this star without being equal to it.
        for g in GeneratorsOfGroup(stab) do
            if not first*g*Inverse(first) in StabU then
                return [false,StabU];
            fi;
        od;
        for cell in cells do
            if not cell[2]*Inverse(first) in StabU then
                return [false,StabU];
            fi;
        od;

        # Each center-stabilizer generator must permute these finitely
        # many translated cells. With left translates, equality means
        # l^-1*g*h lies in the stabilizer of the orbit representative.
        for g in GeneratorsOfGroup(StabU) do
            for cell in cells do
                h:=cell[2];
                if not ForAny(cells,other->Inverse(other[2])*g*h in stab) then
                    return [false,StabU];
                fi;
            od;
        od;
    od;
od;
return [true,StabU];
end);

InstallGlobalFunction(IsGBalanced, function(L,e,n,k)
return PreIsGBalanced(L,e,n,k)[1];
end);

InstallGlobalFunction(Stabstar, function(L,e,n,k)
return PreIsGBalanced(L,e,n,k)[2];
end);
