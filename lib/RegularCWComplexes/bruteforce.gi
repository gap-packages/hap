#RT:=0;
#####################################################################
#####################################################################
InstallGlobalFunction(GComplexToFiniteRegularCWRegion,
function(R,N)
local Y, dim, dim1, Cells,boundaries, StandardForm, 
Tiles, tile, pos,NeighbourGens,CanonicalLeftCosetElement,TileReps,
w,NewLeaves,Leaves,i,e,f,n,k,g,b,x;
#R is a contractible G-complex

#######################################################
dim:=0;
for n in [1..Length(R)] do
if R!.dimension(n)>0 then dim:=dim+1;
else break; fi;
od;
dim1:=dim-1;
#######################################################
#
if R!.dimension(dim)>1 then
Print("This function currently works only on G-complexes with a single orbit of top-dimensional cells.\n");
return fail;
fi;

########################################
StandardForm:=function(n,x)
local e,g,gc,stab;
e:=AbsInt(x[1]);
g:=x[2];
stab:=R!.stabilizer(n,e);
#RT:=RT-Runtime();   #THIS TAKES MOST OF THE TIME
gc:=CanonicalLeftCosetElement(g,stab);
#RT:=RT+Runtime();
return [e,gc];
end;
########################################

if IsMatrix(R!.elts[1]) then
########################################
CanonicalLeftCosetElement:=function(g,U)
local i,j,n,M,UT,m,L;
if not IsBound(U!.transposes) then
U!.transposes:= List(Elements(U),A->TransposedMat(A));
fi;
n:=Length(U!.transposes[1]);
UT:=U!.transposes;
for i in [1..n] do
for j in [1..n] do
M:=List(UT,u->g[i]*u[j]);
m:=Minimum(M);
L:=Filtered([1..Length(UT)],k->M[k]=m);
UT:=UT{L};
if Length(UT)=1 then return g*TransposedMat(UT[1]); fi;
od;od;
end;
########################################
else
########################################
CanonicalLeftCosetElement:=function(g,U);
#return CanonicalRightCosetElement(U,g^-1)^-1;
#return SortedList(g*Elements(U))[1];
return Minimum(g*Elements(U));
end;
########################################
fi;

Tiles:=[];
for k in [1..R!.dimension(dim)] do
   tile:=1*R!.boundary(dim,k);
   Apply(tile,x->[x[1],R!.elts[x[2]]]);
   Apply(tile,e->StandardForm(dim1,e));
   Add(Tiles,tile);
od;

###################################################################
NeighbourGens:=List([1..Length(Tiles)],x->[]);
for i in [1..Length(Tiles)] do
tile:=Tiles[i];

   for e in tile do
      for f in tile do
         if e[1]=f[1] then
            for g in R!.stabilizer(dim1,f[1]) do
               x:=e[2]*g*f[2]^-1;

               x:=CanonicalLeftCosetElement(x,R!.stabilizer(dim,i));
               if not x in R!.stabilizer(dim,i) then
                  Add(NeighbourGens[i],x);
               fi;
            od;
         fi;
      od;
   od;
od;

Apply(NeighbourGens,x->SSortedList(x));

TileReps:=List([1..Length(Tiles)],i->[One(R!.stabilizer(dim,i))]);
Leaves:=List([1..Length(Tiles)],i->[One(R!.stabilizer(dim,i))]);
for i in [1..Length(Tiles)] do
NewLeaves:=[];
for k in [1..N] do
   for x in Leaves[i] do
      for g in NeighbourGens[i] do
         w:=CanonicalLeftCosetElement(g*x,R!.stabilizer(dim,i));
         if not w in TileReps[i] then
            AddSet(TileReps[i],w);
            Add(NewLeaves,w);
         fi;
      od;
   od;
   Leaves[i]:=NewLeaves;
   NewLeaves:=[];
od;
od;

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
Y!.G_Cells:=Cells;
return Y;
end);
#####################################################################
#####################################################################

