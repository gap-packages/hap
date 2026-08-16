
#####################################################################
#####################################################################
InstallGlobalFunction(GComplexToRegularCWComplex,
function(R,N)
local Y, dim, dim1, Cells,boundaries, StandardForm, 
Tiles, tile, pos,NeighbourGens,TileReps,
w,NewLeaves,Leaves,i,e,f,n,k,g,b,x,T;
#R is a contractible G-complex
#S is a set of elements in G=R!.group
#N is a finite index subgroup of G=R!.group

if not R!.group=N then
T:=RightTransversal(R!.group,N);
else
x:=Group(One(R!.group));
T:=RightTransversal(x,x);
T!.group:=R!.group;
T!.subgroup:=R!.group;
T!.poscan:=function(a) return 1; end;
fi;

#######################################################
dim:=0;
for n in [1..Length(R)] do
if R!.dimension(n)>0 then dim:=dim+1;
else break; fi;
od;
dim1:=dim-1;
#######################################################

########################################
StandardForm:=function(n,x)
local e,g,gc,stab;

e:=AbsInt(x[1]);
g:=x[2];
stab:=R!.stabilizer(n,e);
gc:=g*Elements(stab);
gc:=List(gc, a->T[T!.poscan(a)]);
gc:=Minimum(gc);
return [e,gc];
end;
########################################


Tiles:=[];
for k in [1..R!.dimension(dim)] do
   tile:=1*R!.boundary(dim,k);
   Apply(tile,x->[x[1],R!.elts[x[2]]]);
   Apply(tile,e->StandardForm(dim1,e));
   Add(Tiles,tile);
od;


TileReps:=List([1..Length(Tiles)],i->T);
##################################################################
##################################################################

########################################
########################################

Cells:=List([1..dim+1],i->[]);
## Cells[k+1] is a list of the k-cells present in the complex

for k in [1..R!.dimension(dim)] do
for g in TileReps[k] do
Add(Cells[dim+1],StandardForm(dim,[k,g]));
od;
od;
Cells[dim+1]:=SSortedList(Cells[dim+1]);

for n in Reversed([1..dim]) do
for x in Cells[n+1] do
b:=1*R!.boundary(n,x[1]);
Apply(b,z->[z[1],R!.elts[z[2]]]);
Apply(b,y->[y[1],x[2]*y[2]]);
Apply(b,y->StandardForm(n-1,y));
Append(Cells[n],b);
od;
Cells[n]:=SSortedList(Cells[n]);
od;


boundaries:=List([1..dim+2], n->[]);
boundaries[1]:=List(Cells[1],x->[1,0]);

for n in [1..dim] do
for k in [1..Length(Cells[n+1])] do
x:=Cells[n+1][k];
b:=1*R!.boundary(n,x[1]);
Apply(b,z->[z[1],R!.elts[z[2]]]);
Apply(b,y->[y[1],x[2]*y[2]]);
Apply(b,y->StandardForm(n-1,y));
Apply(b,c->PositionSorted(Cells[n],c));
b:=Concatenation([Length(b)],b);
boundaries[n+1][k]:=b;
od;
od;

Y:=RegularCWComplex(boundaries);

return Y;
end);
#####################################################################
#####################################################################

