InstallGlobalFunction(ForgetG, function(L,U,n)
local FL, Y, k, d, cell, boundarycell, dcell, orbits, grouped, LO, i, g, g1, h, gh, t, r, z, Le, Stab, glist, j;

k:=Length(U)-1;
FL:=[];
FL[1]:=[[1,0]];

for d in [2..k+1] do
    Le:=Length(U[d]);
    FL[d]:=List([1..Le], i -> []);
    for r in [1..Le] do
        cell:=U[d][r];
        t:=0; #t counts how many cells are in the boundary
        FL[d][r][1]:=1; #Platzhalter
        boundarycell:=L!.boundary(d+n-1,cell[1]);
        orbits:=Set(List(U[d-1],x->x[1]));
        grouped:=List(orbits, y -> Filtered(U[d-1], x -> x[1] = y));
        for dcell in boundarycell do
            Stab:=L!.stabilizer(d+n-2,AbsoluteValue(dcell[1]));
            h:=L!.elts[dcell[2]];
            g:=cell[2];
            gh:=g*h;
            LO:=Length(orbits);
            z:=0;
            for i in [1..LO] do
                if i>1 then
                    z:=z+Length(grouped[i-1]); #counts number of cell
                fi;
                if AbsoluteValue(dcell[1])=orbits[i] then
                    glist:=(List(grouped[i],y->y[2]));
                    for j in [1..Length(glist)] do
                        g1:=glist[j];
                        if Inverse(g1)*gh in Stab then #getauscht von g1*Inverse(gh) zu Inverse(g1)*gh
                            Add(FL[d][r],z+j);
                            t:=t+1;
                        fi;
                    od;
                fi;
            od;
        od;
        FL[d][r][1]:=t;
    od;
od;

Add(FL,[]);
#Y:=RegularCWComplex(FL);
#cobnd:=Y!.coboundaries;

#return [FL,cobnd];
return FL;
end);
###################################

