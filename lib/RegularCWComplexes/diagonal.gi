
#######################################
#######################################
InstallGlobalFunction(DiagonalCWMap,
function(X)
local XX, CX, CXX, D, map, chainHomotopy,
      CellDiagonal, 2CellDiagonal, CellDiagRec, q2pX,nn,kk;

CX:=ChainComplexOfRegularCWComplex(X);
XX:=DirectProductOfRegularCWComplexes(X,X,Dimension(X));  #THIS NEEDS A LAZY IMPLEMENTATION
CXX:=ChainComplexOfRegularCWComplex(XX); #THIS ALSO NEEDS A LAZY IMPLEMENTATION
q2pX:=XX!.quad2pair;
CellDiagRec:=List([0..Dimension(X)],i->[]);

###########################
###########################
CellDiagonal:=function(n,k)
local f,F,E,EE,CEE,p2qE,ff,FF,mapff,d, bnd,u,i,
CE, Homotopy, FFmat, x,GG,GGrec;

if IsBound(CellDiagRec[n][k]) then return CellDiagRec[n][k]; fi;

#This function returns the image vector in C(XxX) of the k-th n-cell of X
#where n>1.

#f: Closure(e^n_k) >---> X
f:=CWSubcomplexToRegularCWMap(ClosureCWCell(X,n,k));
#F:C(e^n_k) >--> C(X)
F:=ChainMapOfRegularCWMap(f);
CE:=Source(F);
#E is the closure of e^n_k
E:=Source(f);
#Create a dvf on E using the following command
CriticalCells(E);;
#EE is the direct product E x E with dvf inherited from E
EE:=DirectProductOfRegularCWComplexes(E,E,1+Dimension(E));
EE!.properties[1][2]:=Dimension(E);
CEE:=ChainComplexOfRegularCWComplex(EE);
p2qE:=EE!.pair2quad;

#ff: E x E  >--> X x X
#FF:=C(E x E) >--> C(X x X)
############################
mapff:=function(n,m)
local k,l,i,j, qd;
qd:=p2qE[n+1][m];
k:=qd[1];
l:=qd[2];
i:=f!.mapping(k-1,qd[3]);
j:=f!.mapping(l-1,qd[4]);
return q2pX[k][l][i][j][2];
end;
############################
ff:=Objectify(HapRegularCWMap,
       rec(
           source:=EE,
           target:=XX,
           mapping:=mapff  ));
FF:=ChainMapOfRegularCWMap(ff);

D:=ChainComplexWithChainHomotopy(EE);
Homotopy:=D!.chainHomotopy;  #A contracting homotopy on CEE

#          map
#   C(X)  ------>  C(X x X) 
#    ^               | ^
#    |F            GG| | FF     GG is defined only on the image of FF
#    |               | |
#    |               V |
#   C(E)  ------>  C(E x E) ---
#                      ^      | D
#                      |      |
#                      --------                          

GGrec:=[];
for i in [1..EE!.nrCells(n-1)] do
x:=mapff(n-1,i);
GGrec[x]:=i;
od;

###########################
#THIS FUNCTION NEEDS IMPROVING!!!!!
GG:=function(v)
local u, i;
u:=[1..CEE!.dimension(n-1)]*0;
for i in Filtered([1..CXX!.dimension(n-1)], x->v[x]<>0) do
u[GGrec[i]]:= v[i];
od;
#u:=SolutionMat(FFmat,v);
return u;
end;
###########################

bnd:=1*CE!.boundary(n,1);
bnd:=F!.mapping(bnd,n-1);
u:=map(bnd,n-1);           #POSSIBLE RECURSION
u:=GG(u);
u:=Homotopy(u,n-1);
CellDiagRec[n][k]:= FF!.mapping(u,n);
Unbind(GGrec); Unbind(GG);Unbind(f);Unbind(FF);Unbind(CEE);Unbind(D);
Unbind(EE);Unbind(Homotopy);
return CellDiagRec[n][k];
end;
###########################
###########################

2CellDiagonal:=function(i); return CellDiagonal(2,i); end;
###########################
##THE FOLLOWING NEEDS TO BE COMPLETELY RE-DESIGNED AND BE BASED 
##ON A MORE GENERAL FORMULA THAN THAT OF KRAVATZ -- ONE BASED ON 
##A GENERAL ORDERING OF THE VERTICES OF A POLYGON
######################
if false then    #####
######################
2CellDiagonal:=function(k)
local u, edges, edgebnds, right, EDGES, N,2N,2N1, pos,
a, i, j,l,ORIENT,word,L,s,t;
#k refers to the k-th cell of dimension 2

if IsBound(CellDiagRec[2][k]) then return CellDiagRec[2][k]; fi;

u:=[1..CXX!.dimension(2)]*0;
edges:=X!.boundaries[3][k];
N:=edges[1];
2N:=2*N;
2N1:=2N-1;
edges:=edges{[2..N+1]};
edgebnds:=List(edges,j->X!.boundaries[2][j]{[2,3]});
EDGES:=[edges[1]];
right:=edgebnds[1][2];
ORIENT:=[];
if right=Maximum(edgebnds[1]) then Add(ORIENT,1);
else Add(ORIENT,-1);
fi;
Remove(edges,1);
Remove(edgebnds,1);

while Length(edges)>0 do
pos:=PositionProperty(edgebnds,x->right in x);
Add(EDGES,edges[pos]);
if edgebnds[pos][1]=right then
    right:=edgebnds[pos][2];
else
    right:=edgebnds[pos][1];
fi;
if right=Maximum(edgebnds[pos]) then Add(ORIENT,1);
else Add(ORIENT,-1);
fi;
Remove(edges,pos);
Remove(edgebnds,pos);
od;

Unbind(edges);
Unbind(edgebnds);

word:=[];
for i in [1..N] do
if ORIENT[i]=1 then Append(word,[EDGES[i],-EDGES[i]]);
else Append(word,[-EDGES[i],EDGES[i]]);
fi;
od;

if word[2N]>0 then
word:=Concatenation([word[2N]],word{[1..2N1]});
fi;

pos:=PositionProperty(word,x->x<0);
L:=[];
while true do
pos:=PositionProperty(word,x->x<0);
pos:=PositionProperty(word,x->x>0,pos);
if pos=fail then break; fi;
i:=word[pos];
j:=word[pos-1];
word[pos]:=j;
word[pos-1]:=i;
Add(L,[i,j]);
od;

Print(word,"\n");
Print(L,"\n");
Apply(L,x->q2pX[2][2][AbsInt(x[2])][AbsInt(x[1])][2]);
s:=Minimum(X!.boundaries[2][EDGES[1]]{[2,3]});
Print(L,"\n");
Add(L,q2pX[3][1][k][s][2]);
Add(L,q2pX[1][3][s][k][2]);

for i in L do
u[i]:=u[i]+1;
od;
Print([N,Length(Filtered(u*TransposedMat(BoundaryMatrix(CXX,2)),z->(z mod 2)<>0))],"   ");
CellDiagRec[2][k]:=u;
return u;
end;
###########################
fi; #######################
###########################

###########################
###########################
##map(w,n):C_nX ---> C_n(XxX) is the diagonal approximation mapping 
#in degree n
map:=function(w,n)
local u,k,i,a,b;
u:=[1..CXX!.dimension(n)]*0;

###########################
if n=0 then
for i in [1..CX!.dimension(0)] do
if not w[i]=0 then
k:=q2pX[1][1][i][i][2];
u[k]:=u[k]+w[i];;
fi;;
od;
return u;
fi;
###########################

###########################
if n=1 then
for i in [1..CX!.dimension(1)] do
if not w[i]=0 then
b:=Minimum(X!.boundaries[2][i]{[2,3]});
a:=Maximum(X!.boundaries[2][i]{[2,3]});
k:=q2pX[1][2][a][i][2];
u[k]:=u[k]+w[i];;
k:=q2pX[2][1][i][b][2];
u[k]:=u[k]+w[i];;
fi;
od;
return u;
fi;
###########################

###########################
if n=2 then
for i in [1..CX!.dimension(2)] do
if not w[i]=0 then
u:=u+w[i]*2CellDiagonal(i);
fi;
od;
return u;
fi;
###########################

###########################
for i in [1..CX!.dimension(n)] do
if not w[i]=0 then
u:=u+w[i]*CellDiagonal(n,i);
fi;
od;
return u;
###########################

end;
###########################
###########################

return Objectify(HapChainMap,
       rec(
           source:=CX,
           target:=CXX,
           mapping:=map,
           properties:=[["characteristic", 0],["type","chainMap"]]));

end);
#######################################
#######################################

#######################################
InstallGlobalFunction(RegularCWSimplex,
function(n)
local Y;
Y:=MaximalSimplicesToSimplicialComplex([[1..n+1]]);
return RegularCWComplex(Y);
end);
#######################################

#######################################
InstallGlobalFunction(RegularCWCube,
function(n)
local A,Y,k;

if n=0 then A:=[ [[1,0]], [] ];
return RegularCWComplex(A);
fi;

if n=1 then A:=[ [[1,0],[1,0]], [[2,1,2]], [] ];
return RegularCWComplex(A);
fi;

A:=[1]; 
for k in [1..n-1] do
A:=[A];
od;
Y:=PureCubicalComplex(A);
return RegularCWComplex(Y);
end);
#######################################

#######################################
InstallGlobalFunction(RegularCWPermutahedron,
function(n)
local A,Y,k;

if n=0 then A:=[ [[1,0]], [] ];
return RegularCWComplex(A);
fi;

if n=1 then A:=[ [[1,0],[1,0]], [[2,1,2]], [] ];
return RegularCWComplex(A);
fi;

A:=[1];
for k in [1..n-1] do
A:=[A];
od;
Y:=PurePermutahedralComplex(A);
return RegularCWComplex(Y);
end);
#######################################

#######################################
InstallGlobalFunction(RegularCWPolygon,
function(n)
local A,B,m;
m:=Int(n/2);
A:=[];
A[1]:=List([1..n],i->[1,0]);

if IsEvenInt(n) then
A[2]:=List([1..m],i->[2,2*(i-1),2*i]);
A[2][1][2]:=1;
B:=List([0..m-1],i->[2,2*i+1,2*(i+1)+1]);
B[m][3]:=n;
Append(A[2],B);
else
A[2]:=List([1..m+1],i->[2,2*(i-1),2*i]);
A[2][1][2]:=1;
A[2][m+1][3]:=n;
B:=List([0..m-1],i->[2,2*i+1,2*(i+1)+1]);
Append(A[2],B);
fi;
A[3]:=[Concatenation([n],[1..n])];
A[4]:=[];
return RegularCWComplex(A);

end);
######################################
