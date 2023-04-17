
####################################################
####################################################
HAP_PoincareCubeManifoldEdgeDegrees:=function(Y)
local s,pos, b;

#Returns a list [[d1,n1], [d2,n2], ..., [dk,nk]]
# where nk is the number of *cube* edges in YY such that
# lie in precisely dk faces.
# NOTE: Y must be constructed using PoincareCubeCWComplexNS. 

s:=List([1..Y!.nrCells(1)],i->Y!.coboundaries[2][i][1]);
s:=Collected(s);;

pos:=PositionProperty(s,x->x[1]=8);;
s[pos][2]:=s[pos][2]-6;;

pos:=PositionProperty(s,x->x[1]=6);;
s[pos][2]:=s[pos][2]-8;;

pos:=PositionProperty(s,x->x[1]=4);;
s[pos][2]:=s[pos][2]-12-24;;

s:=Filtered(s,x->x[2]<>0);
s:=List(s,x->[x[1]/2,x[2]]);
return SortedList(s);;
end;
#####################################################
####################################################


####################################################
####################################################
InstallGlobalFunction(BarycentricallySimplifiedComplex,
function(Y)
local V,W,s,t;
W:=SimplifiedComplex(Y);
s:=Size(W);
t:=0;
while s> t do
V:=W;
s:=Size(V);
W:=SimplifiedComplex(RegularCWComplex(BarycentricSubdivision(V)));
t:=Size(W);
od;
if IsBound(Y!.cubeFacePairings) then
V!.cubeFacePairings:=Y!.cubeFacePairings;
fi;
return V;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(PoincareCubeCWComplex,
function(arg)
local Y, YY;

if Length(arg)=2 then
Y:=PoincareCubeCWComplexNS(arg[1],arg[2]);
else
Y:=PoincareCubeCWComplexNS(arg[1],arg[2],arg[3]);
fi;
YY:=SimplifiedComplex(Y);
YY!.cubeFacePairings:=Y!.cubeFacePairings;
YY!.edgeDegrees:=Y!.edgeDegrees;
return YY;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(PoincareCubeCWComplexNS,
function(arg)
local A, G, L, M, N, Y, B, F, data, NF, NFF, Vertices, Edges, Faces, 
      CubeEdges, CubeFaces, bool, act, p, q, i,
      X, K, KK, DD,SD,TD,x, cnt,pos, C,v,n;

#First input variable is a pairing of the integers 1..6
#such as A:=[ [1,6], [2,3], [4,5]  ];
#Second variable is a list B:=[a,b,c] of three elements of 
#the dihedral group D_8 < Sym([1..4]);

#######################################
#Action of an element g of the dihedral group of order 8
#on a list x=[x1,x2,x3,x4] of length 4
act:=function(g,x);
return List([1..4],i->x[i^g]);
end;
#######################################

#######################################
#Y is the 3-dimensional cube as a regular CW-complex.
B:=[];
B[1]:=List([1..8],i->[1,0]);
B[2]:=[[2,1,2],[2,3,4],[2,7,8],[2,5,6],
       [2,1,4],[2,2,3],[2,6,7],[2,5,8],
       [2,1,5],[2,4,8],[2,3,7],[2,2,6]];
B[3]:=[[4,1,4,9,12],[4,2,3,10,11],[4,1,2,5,6],[4,3,4,7,8],
       [4,5,8,9,10],[4,6,7,11,12]];
B[4]:=[[6,1,2,3,4,5,6]];
B[5]:=[];
Y:=RegularCWComplex(B);

K:=BarycentricSubdivision(Y);
KK:=RegularCWComplex(K);
#######################################

#data is a list of three lists which describe the identification.
data:=[ [ ], [ ], [ ], [ ] ];

#      5 6
#    5 1 2 6
#    8 4 3 7
#      8 7
#      5 6

CubeEdges:= [ [1,2], [3,4], [7,8], [5,6],[1,4],[2,3],[6,7],[5,8],
              [1,5], [4,8], [3,7], [2,6]];
CubeFaces:=[ [1,5,6,2], [3,7,8,4], [1,2,3,4], [5,8,7,6], [1,4,8,5], [2,6,7,3]];
#Note: CubeEdges/Vertices have the same order as the boundary B above.
F:=CubeFaces;

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
if T in CubeFaces then return T;
else return NF(Reversed(F));
fi;
end;
###########################################


if Length(arg)=0 then
Print("This function returns a CW-complex arising as a quotient of a cube with vertices \n [1,5,6,2] in the plane z=1 (face 1) \n [3,7,8,4] in the plane z=0 (face 2)\n [1,2,3,4] in the plane y=0 (face 3)\n [5,8,7,6] in the plane y=1 (face 4)\n [1,4,8,5] in the plane x=0 (face 5)\n [2,6,7,3] in the plane x=1 (face 6).\n\nThe function can be applied to two inputs PoincareCubeCWComplex(A,G) where A is a pairing of the six faces such as A=[[1,2],[3,4],[5,6]] and G is a list of three elements of the dihedral group D8 such as G=[(2,4),(2,4),(2,4)*(1,3)].\n\nAlternatively the function can be applied to three inputs PoincareCubeCWComplex(L,M,N) where each input is a pairing of faces specified by a bijection of face vertices such as L=[[1,2,3,4], [6,7,5,3]].\n\n");
return fail;
fi;

if Length(arg)=2 then
A:=arg[1];
G:=arg[2];

L:= [ F[A[1][1]], act( G[1], F[A[1][2]] ) ];
M:= [ F[A[2][1]], act( G[2], F[A[2][2]] ) ];
N:= [ F[A[3][1]], act(G[3], F[A[3][2]] ) ];
fi;

if Length(arg)=3 then
L:= arg[1];
M:= arg[2];
N:= arg[3];

A:=[];
for x in [L,M,N] do
Add(A, [Position(F,NFF(x[1])),Position(F,NFF(x[2]))]);
od;

fi;

###########################################
#Let's check that the input arguments are genuine faces of our cube.
#This check is not needed but might catch a programming slip.
bool:=false;
if not NFF(L[1]) in CubeFaces then bool:=true; fi;
if not NFF(L[2]) in CubeFaces then bool:=true; fi;
if not NFF(M[1]) in CubeFaces then bool:=true; fi;
if not NFF(M[2]) in CubeFaces then bool:=true; fi;
if not NFF(N[1]) in CubeFaces then bool:=true; fi;
if not NFF(N[2]) in CubeFaces then bool:=true; fi;
if bool then
Print("Arguments must be pairs of cube faces.\n"); return fail;
fi;
###########################################


cnt:=0;
for X in [L,M,N] do
Vertices:=[];
Edges:=[];
Faces:=[];
cnt:=cnt+1;
SD:=1*K!.simplicesLst[3];
SD:=Filtered(SD,x->A[cnt][1]+20 in x and not 27 in x);

Add(Faces, [Position(CubeFaces,NFF(X[1])), Position(CubeFaces,NFF(X[2]))]);
TD:=List(SD,x->[x[1],x[2],Faces[1][2]+20]);


Add(Vertices, [X[1][1], X[2][1]]);
Add(Vertices, [X[1][2], X[2][2]]);
Add(Vertices, [X[1][3], X[2][3]]);
Add(Vertices, [X[1][4], X[2][4]]);

for x in TD do
pos:=PositionProperty(Vertices,v->v[1]=x[1]);
x[1]:=Vertices[pos][2]*1;
od;

p:=Position(CubeEdges,SortedList(X[1]{[1,2]}));
q:=Position(CubeEdges,SortedList(X[2]{[1,2]}));
Add(Edges,([p,q]));
p:=Position(CubeEdges,SortedList(X[1]{[2,3]}));
q:=Position(CubeEdges,SortedList(X[2]{[2,3]}));
Add(Edges,([p,q]));
p:=Position(CubeEdges,SortedList(X[1]{[3,4]}));
q:=Position(CubeEdges,SortedList(X[2]{[3,4]}));
Add(Edges,([p,q]));
p:=Position(CubeEdges,SortedList(X[1]{[1,4]}));
q:=Position(CubeEdges,SortedList(X[2]{[1,4]}));
Add(Edges,([p,q]));

for x in TD do
pos:=PositionProperty(Edges,v->v[1]=x[2]-8);
x[2]:=Edges[pos][2]+8;
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

#Y:=SimplifiedComplex(RegularCWComplex(B));
Y:=RegularCWComplex(B);
Y!.cubeFacePairings:=[L,M,N];
Y!.edgeDegrees:=HAP_PoincareCubeManifoldEdgeDegrees(Y);
return Y;
end);
###################################################
###################################################

