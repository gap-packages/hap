#######################################################
#######################################################
InstallGlobalFunction(RegularCWPolytope,
function(points)
local polytope, H, ind, NrCells, Boundaries, Coboundaries, 
      Properties,dim, n, F, low, hi, HI, LOW, i, j;

polytope:=CreatePolymakeObject();
AppendPointlistToPolymakeObject(polytope,points);
H:=Polymake(polytope,"HASSE_DIAGRAM");
if Length(H.hasse[1][1])=0 then
ind:=H.faceindices;
else
ind:=Reversed(H.faceindices);
fi;
dim:=Length(H.faceindices)-2;
H:=H.hasse{[2..Length(H.hasse)-1]};

################
NrCells:=function(n);
if n<0 or n>Length(ind)-2 then return 0; fi;
return Length(ind[n+2]);
end;
################

Coboundaries:=[];
for n in [1..dim] do
low:=ind[n+2][1];hi:=ind[n+2][Length(ind[n+2])];
F:=Filtered(H,x->x[2][1]>=low and x[2][1]<=hi);
Apply(F,x->x[2]-low+1);
Coboundaries[n]:=F;
od;
Coboundaries[dim+1]:=List([1..NrCells(dim)],i->0);


Boundaries:=[];
Boundaries[1]:=List([1..NrCells(0)],i->[1,0]);
for n in [1..dim] do
Boundaries[n+1]:=List([1..NrCells(n)],i->[]);
for i in [1..NrCells(n-1)] do
for j in Coboundaries[n][i] do
Add(Boundaries[n+1][j],i);
od;
od;
od;

for n in [1..dim] do
for i in [1..NrCells(n)] do
Boundaries[n+1][i]:= Reversed(Boundaries[n+1][i]);
Add(Boundaries[n+1][i],Length(Boundaries[n+1][i]));
Boundaries[n+1][i]:= Reversed(Boundaries[n+1][i]);
od;
od;

for n in [1..dim] do
for i in [1..NrCells(n-1)] do
Coboundaries[n][i]:= Reversed(Coboundaries[n][i]);
Add(Coboundaries[n][i],Length(Coboundaries[n][i]));
Coboundaries[n][i]:= Reversed(Coboundaries[n][i]);
od;
od;
Properties:=[["dimension",dim]];

return Objectify(HapRegularCWComplex,
       rec(
           nrCells:=NrCells,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           orientation:=fail,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           properties:=Properties));
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(RegularCWOrbitPolytope,
function(G,v)
local points;

if IsPermGroup(G) then
points:=Orbit(G,v,Permuted);;
fi;

if IsMatrixGroup(G) then
points:=Orbit(G,v,OnRight);;
fi;

return RegularCWPolytope(points);
end);
#######################################################
#######################################################

