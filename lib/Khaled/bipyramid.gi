#written by Khaled Alzobydi
###################################################
###################################################
InstallGlobalFunction(PoincareBipyramidCWComplex,
function(arg)
local A, G, L, M, N, P, Y, B, F, data, NF, NFF, Vertices, Edges, Faces,
      BipyramidEdges, BipyramidFaces, bool, act, p, q, i, j,
      X, K, KK, DD,SD,TD,x, cnt,pos, C,v,n, nV, nE, nF;
N:=Length(arg);
####################################################
####################################################
#######################################
#Y is the 3-dimensional Bipyramid as a regular CW-complex.

B:=[];
B[1]:=List([1..N+2],i->[1,0]);
B[2]:=[];
  for i in [1..N] do
  Add(B[2], [2,i,N+1]);
  od;
  for i in [1..N-1] do
  Add(B[2], [2,i,i+1]);
  od;
  Add(B[2], [2,1,N]); 
  for i in [1..N] do
  Add(B[2], [2,i,N+2]);
  od;
B[3]:=[];

  for i in [1..N-1] do
  Add(B[3], [3,i,i+1,N+i]);
  od;
  Add(B[3], [3,1,N,2*N]);

  for i in [1..N-1] do
  Add(B[3], [3,2*N+i,2*N+i+1,N+i]);
  od;
  Add(B[3], [3,2*N+1,3*N,2*N]);

B[4]:=[Concatenation([2*N],[1..2*N])];
B[5]:=[];
Y:=RegularCWComplex(B);
K:=BarycentricSubdivision(Y);
KK:=RegularCWComplex(K);
#######################################
#data is a list of N/2 +1 lists (with n even) which describe the identification.
data:=[ [ ], [ ], [ ], [ ] ];

#                       1 ----- 2             
#      N+1------------- .       | ------------N+2          our Bipyramid
#                       .       |              
#                       N ..... 3              

BipyramidEdges:=[];
    for i in [1..N] do
    Add(BipyramidEdges,[i,N+1]);
    od;
    for i in [1..N-1] do
    Add(BipyramidEdges, [i,i+1]);
    od;
    Add(BipyramidEdges, [1,N]); 
    for i in [1..N] do
    Add(BipyramidEdges, [i,N+2]);
    od;

BipyramidFaces:=[];
  for i in [1..N-1] do
    Add(BipyramidFaces, [i,i+1,N+1]);
  od;
  Add(BipyramidFaces, [1,N,N+1]);
  for i in [1..N-1] do
    Add(BipyramidFaces, [i,i+1,N+2]);
  od;
  Add(BipyramidFaces, [1,N,N+2]);

F:=BipyramidFaces;
#Note: Bipyramid Edges/Vertices have the same order as the boundary B above.

###########################################
###########################################
#Normal form for faces
NF:=function(F)
local m;
m:=Minimum(F);
m:=Position(F,m);
return Concatenation(F{[m..Length(F)]},F{[1..m-1]}   );
end;
###########################################

###########################################
#Normal form for faces incorporating orientation
NFF:=function(F)
local T;
T:=NF(F);
if T in BipyramidFaces then return T;
else return NF(Reversed(F));
fi;
end;
###########################################
A:=[];
for x in arg do
Add(A, [Position(F,NFF(x[1])),Position(F,NFF(x[2]))]);
od;

###########################################
#Let's check that the input arguments are genuine faces of our Bipyramid.
#This check is not needed but might catch a programming slip.
bool:=false;
for L in arg do
if not NFF(L[1]) in BipyramidFaces then bool:=true; fi;
if not NFF(L[2]) in BipyramidFaces then bool:=true; fi;
od;
if bool then
Print("Arguments must be pairs of Bipyramid faces.\n"); return fail;
fi;
###########################################
#############################################
###########################################

nV:=N+2;  #number of vertices
nE:=3*N;  #number of edges
nF:=2*N;  #number of faces
cnt:=0;
for X in arg do
Vertices:=[];
Edges:=[];
Faces:=[];
cnt:=cnt+1;

Add(Faces, [Position(BipyramidFaces,NFF(X[1])), Position(BipyramidFaces,NFF(X[2]))]);

for i in [1..Length(X[1])] do
Add(Vertices, [X[1][i], X[2][i]]);
od;

for i in [1..Length(X[1])-1] do
p:=Position(BipyramidEdges,SortedList(X[1]{[i,i+1]}));
q:=Position(BipyramidEdges,SortedList(X[2]{[i,i+1]}));
Add(Edges,([p,q]));
od;
p:=Position(BipyramidEdges,SortedList(X[1]{[1,Length(X[1])]}));
q:=Position(BipyramidEdges,SortedList(X[2]{[1,Length(X[1])]}));
Add(Edges,([p,q]));


SD:=1*K!.simplicesLst[3];
SD:=Filtered(SD,x->A[cnt][1]+nV+nE in x and not 1+nV+nE+nF in x);
TD:=List(SD,x->[x[1],x[2],Faces[1][2]+nV+nE]);


for x in TD do
pos:=PositionProperty(Vertices,v->v[1]=x[1]);
x[1]:=Vertices[pos][2]*1;
od;

for x in TD do
pos:=PositionProperty(Edges,v->v[1]=x[2]-nV);
x[2]:=Edges[pos][2]+nV;
od;

for i in [1..Length(SD)] do
Add(data[3],[Position(K!.simplicesLst[3],SD[i]), Position(K!.simplicesLst[3],TD[i])]);

for j in [1..Length(SD[i])] do
Add(data[1],[Position(K!.simplicesLst[1],SD[i]{[j]}), Position(K!.simplicesLst[1],TD[i]{[j]})]);
od;

for j in [1..Length(SD[i])-1] do
Add(data[2],[Position(K!.simplicesLst[2],SD[i]{[j,j+1]}), Position(K!.simplicesLst[2],TD[i]{[j,j+1]})]);
od;
Add(data[2],[Position(K!.simplicesLst[2],SD[i]{[1,Length(SD[i])]}), Position(K!.simplicesLst[2],TD[i]{[1,Length(SD[i])]})]);


od;

od;

C:=Target(QuotientChainMap(KK,data));

B:=[];
B[1]:=List([1..C!.dimension(0)],i->[1,0]);

for n in [1..3] do
B[n+1]:=[];
for i in [1..C!.dimension(n)] do
v:=C!.boundary(n,i);
v:=Filtered([1..Length(v)],j -> not IsZero(v[j]));
v:=Concatenation([Length(v)],v);
Add(B[n+1],v);
od;

od;
Add(B,[]);

Y:=SimplifiedComplex(RegularCWComplex(B));
#Y:=RegularCWComplex(B);
Y!.BipyramidFacePairings:=arg;
return Y;
end);
###################################################
###################################################
