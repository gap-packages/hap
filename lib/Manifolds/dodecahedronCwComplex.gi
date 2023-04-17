

####################################################
####################################################
InstallGlobalFunction(PoincareDodecahedronCWComplex,
function(arg)
local A, G, L, M, N, P, Q,R,S,Y, B, F, data, NF, NFF, Vertices, Edges, Faces, 
      DodecahedronEdges, DodecahedronFaces, bool, act, p, q, i,j,
      nV,nE,nF,X, K, KK, DD,SD,TD,x, cnt,pos, C,v,n;

#######################################
#Y is the 3-dimensional octhedron as a regular CW-complex.
B:=[];
B[1]:=List([1..20],i->[1,0]);
B[2]:=[[2,1,2],[2,2,3],[2,3,4],[2,4,5],[2,1,5],
       [2,1,11],[2,2,12],[2,3,13],[2,4,14],[2,5,15],
       [2,6,7],[2,7,8],[2,8,9],[2,9,10],[2,6,10],
       [2,6,16],[2,7,17],[2,8,18],[2,9,19],[2,10,20],
       [2,11,16],[2,12,16],[2,12,17],[2,13,17],[2,13,18],
       [2,14,18],[2,14,19],[2,15,19],[2,15,20],[2,11,20]];
      
       
B[3]:=[[5,1,2,3,4,5],[5,11,12,13,14,15],
       [5,1,7,22,21,6],[5,2,8,24,23,7],[5,3,9,26,25,8],
       [5,4,10,28,27,9],[5,10,29,30,6,5], 
       [5,21,16,15,20,30],[5,22,23,17,11,16],[5,24,25,18,12,17],
       [5,26,27,19,13,18],[5,28,29,20,14,19]];
B[4]:=[[12,1,2,3,4,5,6,7,8,9,10,11,12]];
B[5]:=[];
Y:=RegularCWComplex(B);

K:=BarycentricSubdivision(Y);
KK:=RegularCWComplex(K);
#######################################

#data is a list of three lists which describe the identification.
data:=[ [ ], [ ], [ ], [ ] ];

              
DodecahedronEdges:=[ [1,2],[2,3],[3,4],[4,5],[1,5],
                   [1,11],[2,12],[3,13],[4,14],[5,15],
                   [6,7],[7,8],[8,9],[9,10],[6,10],
                   [6,16],[7,17],[8,18],[9,19],[10,20],
                   [11,16],[12,16],[12,17],[13,17],[13,18],
                   [14,18],[14,19],[15,19],[15,20],[11,20]];
                   
                  
DodecahedronFaces:=[[1,2,3,4,5],[6,7,8,9,10],
                  [1,2,12,16,11],[2,3,13,17,12],[3,4,14,18,13],
                  [4,5,15,19,14],[1,5,15,20,11],
                  [6,10,20,11,16],[6,16,12,17,7],[7,17,13,18,8],
                  [8,18,14,19,9],[9,19,15,20,10]];
F:=DodecahedronFaces;
#Note: Dodecahedron Edges/Vertices have the same order as the boundary B above.

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
if T in DodecahedronFaces then return T;
else return NF(Reversed(F));
fi;
end;
###########################################

L:= arg[1];
M:= arg[2];
N:= arg[3];
P:= arg[4];
Q:= arg[5];
R:= arg[6];

###########################################
#Let's check that the input arguments are genuine faces of our dodecahedron.
#This check is not needed but might catch a programming slip.
bool:=false;
if not NFF(L[1]) in DodecahedronFaces then bool:=true; fi;
if not NFF(L[2]) in DodecahedronFaces then bool:=true; fi;
if not NFF(M[1]) in DodecahedronFaces then bool:=true; fi;
if not NFF(M[2]) in DodecahedronFaces then bool:=true; fi;
if not NFF(N[1]) in DodecahedronFaces then bool:=true; fi;
if not NFF(N[2]) in DodecahedronFaces then bool:=true; fi;
if not NFF(P[1]) in DodecahedronFaces then bool:=true; fi;
if not NFF(P[2]) in DodecahedronFaces then bool:=true; fi;
if bool then
Print("Arguments must be pairs of dodecahedron faces.\n"); return fail;
fi;
###########################################

nV:=20;  #number of vertices
nE:=30;  #number of edges
nF:=12 ; #number of facets

A:=[];
for x in [L,M,N,P,Q,R] do
Add(A, [Position(F,NFF(x[1])),Position(F,NFF(x[2]))]);
od;

cnt:=0;
for X in arg do
Vertices:=[];
Edges:=[];
Faces:=[];
cnt:=cnt+1;

Add(Faces, [Position(DodecahedronFaces,NFF(X[1])), Position(DodecahedronFaces,NFF(X[2]))]);

for i in [1..Length(X[1])] do
Add(Vertices, [X[1][i], X[2][i]]);
od;

for i in [1..Length(X[1])-1] do
p:=Position(DodecahedronEdges,SortedList(X[1]{[i,i+1]}));
q:=Position(DodecahedronEdges,SortedList(X[2]{[i,i+1]}));
Add(Edges,([p,q]));
od;
p:=Position(DodecahedronEdges,SortedList(X[1]{[1,Length(X[1])]}));
q:=Position(DodecahedronEdges,SortedList(X[2]{[1,Length(X[1])]}));
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
Y!.dodecahedronFacePairings:=[L,M,N,P,Q,R];
return Y;
end);
###################################################
###################################################

