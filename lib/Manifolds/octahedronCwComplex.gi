

####################################################
####################################################
InstallGlobalFunction(PoincareOctahedronCWComplex,
function(arg)
local A, G, L, M, N, P, Y, B, F, data, NF, NFF, Vertices, Edges, Faces, 
      OctahedronEdges, OctahedronFaces, bool, act, p, q, i,
      X, K, KK, DD,SD,TD,x, cnt,pos, C,v,n;

#######################################
#Y is the 3-dimensional octhedron as a regular CW-complex.
B:=[];
B[1]:=List([1..6],i->[1,0]);
B[2]:=[[2,1,4],[2,3,4],[2,2,3],[2,1,2],
       [2,1,5],[2,4,5],[2,3,5],[2,2,5],
       [2,1,6],[2,4,6],[2,3,6],[2,2,6]];
B[3]:=[[3,1,5,6],[3,2,6,7],[3,3,7,8],[3,4,5,8],
       [3,1,9,10],[3,2,10,11],[3,3,11,12],[3,4,9,12]];
B[4]:=[[8,1,2,3,4,5,6,7,8]];
B[5]:=[];
Y:=RegularCWComplex(B);

K:=BarycentricSubdivision(Y);
KK:=RegularCWComplex(K);
#######################################

#data is a list of three lists which describe the identification.
data:=[ [ ], [ ], [ ], [ ] ];

#           5
#          /|\ .
#         / | \ .
#        /  |  \  .
#       /  2 -- \ -3
#      1 ------- 4 .   Our octahedron
#       \  |    / .
#        \ \   / .
#         \ | / .
#          \|/ .
#           6.

              
OctahedronEdges:=[[1,4],[3,4],[2,3],[1,2],
                  [1,5],[4,5],[3,5],[2,5],
                  [1,6],[4,6],[3,6],[2,6]];
OctahedronFaces:=[[1,4,5],[3,4,5],[2,3,5],[1,2,5],
                  [1,4,6],[3,4,6],[2,3,6],[1,2,6]];
F:=OctahedronFaces;
#Note: Octahedron Edges/Vertices have the same order as the boundary B above.

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
if T in OctahedronFaces then return T;
else return NF(Reversed(F));
fi;
end;
###########################################

if Length(arg)=2 then

#######################################
#Action of an element g of the dihedral group of order 8
#on a list x=[x1,x2,x3,x4] of length 4
act:=function(g,x);
return List([1..3],i->x[i^g]);
end;
#######################################

A:=arg[1];
G:=arg[2];

L:= [ F[A[1][1]], act( G[1], F[A[1][2]] ) ];
M:= [ F[A[2][1]], act( G[2], F[A[2][2]] ) ];
N:= [ F[A[3][1]], act( G[3], F[A[3][2]] ) ];
P:= [ F[A[4][1]], act( G[4], F[A[4][2]] ) ];
fi;

if Length(arg)=4 then

L:= arg[1];
M:= arg[2];
N:= arg[3];
P:= arg[4];

A:=[];
for x in [L,M,N,P] do
Add(A, [Position(F,NFF(x[1])),Position(F,NFF(x[2]))]);
od;

fi;

###########################################
#Let's check that the input arguments are genuine faces of our octahedron.
#This check is not needed but might catch a programming slip.
bool:=false;
if not NFF(L[1]) in OctahedronFaces then bool:=true; fi;
if not NFF(L[2]) in OctahedronFaces then bool:=true; fi;
if not NFF(M[1]) in OctahedronFaces then bool:=true; fi;
if not NFF(M[2]) in OctahedronFaces then bool:=true; fi;
if not NFF(N[1]) in OctahedronFaces then bool:=true; fi;
if not NFF(N[2]) in OctahedronFaces then bool:=true; fi;
if not NFF(P[1]) in OctahedronFaces then bool:=true; fi;
if not NFF(P[2]) in OctahedronFaces then bool:=true; fi;
if bool then
Print("Arguments must be pairs of octahedron faces.\n"); return fail;
fi;
###########################################


cnt:=0;
for X in [L,M,N,P] do
Vertices:=[];
Edges:=[];
Faces:=[];
cnt:=cnt+1;
SD:=1*K!.simplicesLst[3];
SD:=Filtered(SD,x->A[cnt][1]+18 in x and not 27 in x);

Add(Faces, [Position(OctahedronFaces,NFF(X[1])), Position(OctahedronFaces,NFF(X[2]))]);
TD:=List(SD,x->[x[1],x[2],Faces[1][2]+18]);


Add(Vertices, [X[1][1], X[2][1]]);
Add(Vertices, [X[1][2], X[2][2]]);
Add(Vertices, [X[1][3], X[2][3]]);

for x in TD do
pos:=PositionProperty(Vertices,v->v[1]=x[1]);
x[1]:=Vertices[pos][2]*1;
od;

p:=Position(OctahedronEdges,SortedList(X[1]{[1,2]}));
q:=Position(OctahedronEdges,SortedList(X[2]{[1,2]}));
Add(Edges,([p,q]));
p:=Position(OctahedronEdges,SortedList(X[1]{[2,3]}));
q:=Position(OctahedronEdges,SortedList(X[2]{[2,3]}));
Add(Edges,([p,q]));
p:=Position(OctahedronEdges,SortedList(X[1]{[1,3]}));
q:=Position(OctahedronEdges,SortedList(X[2]{[1,3]}));
Add(Edges,([p,q]));


for x in TD do
pos:=PositionProperty(Edges,v->v[1]=x[2]-6);
x[2]:=Edges[pos][2]+6;
od;

for i in [1..Length(SD)] do
Add(data[3],[Position(K!.simplicesLst[3],SD[i]), Position(K!.simplicesLst[3],TD[i])]);

Add(data[1],[Position(K!.simplicesLst[1],SD[i]{[1]}), Position(K!.simplicesLst[1],TD[i]{[1]})]);

Add(data[1],[Position(K!.simplicesLst[1],SD[i]{[2]}), Position(K!.simplicesLst[1],TD[i]{[2]})]);

Add(data[1],[Position(K!.simplicesLst[1],SD[i]{[3]}), Position(K!.simplicesLst[1],TD[i]{[3]})]);

Add(data[2],[Position(K!.simplicesLst[2],SD[i]{[1,2]}), Position(K!.simplicesLst[2],TD[i]{[1,2]})]);

Add(data[2],[Position(K!.simplicesLst[2],SD[i]{[1,3]}), Position(K!.simplicesLst[2],TD[i]{[1,3]})]);

Add(data[2],[Position(K!.simplicesLst[2],SD[i]{[2,3]}), Position(K!.simplicesLst[2],TD[i]{[2,3]})]);


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
Y!.octahedronFacePairings:=[L,M,N,P];
return Y;
end);
###################################################
###################################################




####################################################
####################################################
InstallGlobalFunction(PoincarePrismCWComplex,
function(arg)
local A, G, L, M, N, P, Y, B, F, data, NF, NFF, Vertices, Edges, Faces,
      PrismEdges, PrismFaces, bool, act, p, q, i, j,
      X, K, KK, DD,SD,TD,x, cnt,pos, C,v,n, nV, nE, nF;

N:=2*(Length(arg)-1);
#######################################
#Y is the 3-dimensional octhedron as a regular CW-complex.
B:=[];
B[1]:=List([1..2*N],i->[1,0]);
B[2]:= List([1..N-1],i->[2,i,i+1]);
  Add(B[2],[2,1,N]);
  for i in [1..N-1] do
    Add(B[2], [2,i+N,i+1+N]);
  od;
  Add(B[2], [2, N+1,N+N]);
  for i in [1..N] do
    Add(B[2],[2,i,N+i]);
  od;

B[3]:=
[ Concatenation([N],[1..N]), Concatenation([N],[1..N]+N) ];
  for i in [1..N-1] do
    Add(B[3], [4, i, 2*N+i+1, N+i, 2*N+i]);
  od;
  Add(B[3], [4,N,2*N+N,N+N,2*N+1]);
B[4]:=[Concatenation([N+2],[1..N+2])];
B[5]:=[];
Y:=RegularCWComplex(B);


K:=BarycentricSubdivision(Y);
KK:=RegularCWComplex(K);
#######################################

#data is a list of N/2 +1 lists (with n even) which describe the identification.
data:=List([0..N],i -> []);

#     1 ----- 2             N+1 ----- N+2
#     .       |  .......     .         |      our prism
#     .       |              .         |
#     N ..... 3             N+N ----- N+3 

PrismEdges:=List([1..N-1], i->[i,i+1]);
  Add(PrismEdges,[1,N]);
  for i in [1..N-1] do
    Add(PrismEdges,[N+i,N+i+1]);
  od;
  Add(PrismEdges,[N+1,N+N]);
  for i in [1..N] do
    Add(PrismEdges,[i,N+i]);
  od;

PrismFaces:=[ [1..N], [1..N]+N ];
  for i in [1..N-1] do
    Add(PrismFaces, [i, i+1, N+i+1, N+i]);
  od;
  Add(PrismFaces, [1,N,N+N,N+1]);

F:=PrismFaces;
#Note: Prism Edges/Vertices have the same order as the boundary B above.

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
if T in PrismFaces then return T;
else return NF(Reversed(F));
fi;
end;
###########################################



A:=[];
for x in arg do
Add(A, [Position(F,NFF(x[1])),Position(F,NFF(x[2]))]);
od;

###########################################
#Let's check that the input arguments are genuine faces of our prism.
#This check is not needed but might catch a programming slip.
bool:=false;
for L in arg do
if not NFF(L[1]) in PrismFaces then bool:=true; fi;
if not NFF(L[2]) in PrismFaces then bool:=true; fi;
od;
if bool then
Print("Arguments must be pairs of prism faces.\n"); return fail;
fi;
###########################################

nV:=2*N; #number of vertices
nE:=3*N;  #number of edges
nF:=N+2;  #number of faces
cnt:=0;
for X in arg do
Vertices:=[];
Edges:=[];
Faces:=[];
cnt:=cnt+1;

Add(Faces, [Position(PrismFaces,NFF(X[1])), Position(PrismFaces,NFF(X[2]))]);

for i in [1..Length(X[1])] do
Add(Vertices, [X[1][i], X[2][i]]);
od;

for i in [1..Length(X[1])-1] do
p:=Position(PrismEdges,SortedList(X[1]{[i,i+1]}));
q:=Position(PrismEdges,SortedList(X[2]{[i,i+1]}));
Add(Edges,([p,q]));
od;
p:=Position(PrismEdges,SortedList(X[1]{[1,Length(X[1])]}));
q:=Position(PrismEdges,SortedList(X[2]{[1,Length(X[1])]}));
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
Y!.prismFacePairings:=arg;
return Y;
end);
###################################################
###################################################



