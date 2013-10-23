#######################################################
#######################################################
InstallGlobalFunction(RegularCWPolytope,
function(points)
local polytope, H, ind, Boundaries, n, i, x, Y;

polytope:=CreatePolymakeObject();
AppendPointlistToPolymakeObject(polytope,points);
H:=PolymakeFaceLattice(polytope,true);
ind:=List(H,h->Minimum(Flat(h)));

for n in [1..Length(ind)] do
H[n]:=H[n]-ind[n]+1;
od;
ind:=List(H{[2..Length(H)]},h->Length(h));

##########################################################
Boundaries:=[];
Boundaries[1]:=List([1..ind[1]],x->[1,0]);
for n in [2..Length(ind)] do
Boundaries[n]:=List([1..ind[n]],x->[]);
for x in [1..Length(H[n])] do
for i in H[n][x] do
Add(Boundaries[n][i],x);
od;
od;
Boundaries[n]:=List(Boundaries[n],b->Concatenation([Length(b)],SSortedList(b)));
od;
Boundaries[n+1]:=[Concatenation([ind[n]],[1..ind[n]])];
Boundaries[n+2]:=[];
##########################################################

Y:=RegularCWComplex(Boundaries);
OrientRegularCWComplex(Y);

return Y;
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

