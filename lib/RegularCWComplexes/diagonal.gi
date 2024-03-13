
########################################
########################################
HAP_Test:=function(F)
local C, D, v, w, n, k, i, Fv, dFv, dv, Fdv;

C:=Source(F);
D:=Target(F);

for n in [4..Length(C)] do
v:=[1..C!.dimension(n)]*0;
for k in [1..C!.dimension(n)] do
v[k]:=1;
#################################
Fv:=F!.mapping(v,n);
dFv:=[1..D!.dimension(n-1)]*0;
for i in [1..Length(Fv)] do
if not 0=Fv[i] then
dFv:=dFv+D!.boundary(n,i)*Fv[i];
fi;
od;
################################

dv:=C!.boundary(n,k);
Fdv:=F!.mapping(dv,n-1);
if not dFv=Fdv then Print("Fail!\n"); return false; fi;

v[k]:=0;
od;
Print("OK for n= ",n,"\n");
od;


return true;
end;
########################################
########################################

#######################################
#######################################
InstallGlobalFunction(DiagonalChainMap,
function(X)
local XX, CX, CXX, D, diag, firstproj, secondproj, map, chainHomotopy,
      CellDiagonal, 2CellDiagonal, CellDiagRec, fproj,sproj,q2pX,nn,kk;

CX:=ChainComplexOfRegularCWComplex(X);
#XX:=DirectProductOfRegularCWComplexes(X,X,Dimension(X));  #THIS NEEDS A LAZY IMPLEMENTATION
XX:=DirectProductOfRegularCWComplexesLazy(X,X);
#XX:=DirectProductOfRegularCWComplexes(X,X);
CXX:=ChainComplexOfRegularCWComplex(XX); #THIS ALSO NEEDS A LAZY IMPLEMENTATION
q2pX:=XX!.quad2pair;
CellDiagRec:=List([0..Dimension(X)],i->[]);

###########################
###########################
CellDiagonal:=function(n,k)
local f,F,E,EE,CEE,p2qE,ff,FF,mapff,d, bnd,u,v,i,ii,b,j,jj,
CE, Homotopy, FFmat, x,GG,GGrec;

if IsBound(CellDiagRec[n][k]) then 

return CellDiagRec[n][k]; fi;
#This function returns the image vector in C(XxX) of the k-th n-cell of X
#where n>1.
#f: Closure(e^n_k) >---> X
f:=CWSubcomplexToRegularCWMap(ClosureCWCell(X,n,k));
#F:C(e^n_k) >--> C(X)
F:=ChainMapOfRegularCWMap(f);
CE:=Source(F);
#E is the closure of e^n_k
E:=Source(f);

if IsBound(X!.directed) then
E!.directed:=[];
for i in [1..E!.nrCells(1)] do
ii:=f!.mapping(1,i);
ii:=X!.directed[ii];
b:=E!.boundaries[2][i];
j:=b[2]; j:=f!.mapping(0,j);
jj:=b[3];jj:=f!.mapping(0,jj);
if ii=[j,jj] then Add(E!.directed,[b[2],b[3]]);  
else Add(E!.directed,[b[3],b[2]]); fi;
od;
fi;
#Create a dvf on E using the following command
CriticalCells(E);;
if Length(CriticalCells(E))<>1 then 
Print("Error in contracting vector field\n"); return fail; 
fi;
#CocriticalCellsOfRegularCWComplex(E,n);
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

v:=map(bnd,n-1);           #POSSIBLE RECURSION

u:=GG(v);

v:=Homotopy(u,n-1);
CellDiagRec[n][k]:= FF!.mapping(v,n);
Unbind(GGrec); Unbind(GG);Unbind(f);Unbind(FF);Unbind(CEE);Unbind(D);
Unbind(EE);Unbind(Homotopy);
return CellDiagRec[n][k];
end;
###########################
###########################

2CellDiagonal:=function(i); return CellDiagonal(2,i); end;

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

diag:= Objectify(HapChainMap,
       rec(
           source:=CX,
           target:=CXX,
           mapping:=map,
           properties:=[["characteristic", 0],["type","chainMap"]]));


###########################
###########################
fproj := function(w,n)
local v,i,j;

v:=[1..CX!.dimension(n)]*0;
for i in [1..Length(v)] do
for j in [1..CX!.dimension(0)] do
v[i]:=v[i]+w[q2pX[n+1][1][i][j][2]];
od;od;

return v;
end;
###########################
###########################

###########################
###########################
sproj := function(w,n)
local v,i,j;

v:=[1..CX!.dimension(n)]*0;
for i in [1..Length(v)] do
for j in [1..CX!.dimension(0)] do
v[i]:=v[i]+w[q2pX[1][n+1][j][i][2]];
od;od;

return v;
end;
###########################
###########################


firstproj:= Objectify(HapChainMap,
       rec(
           source:=CXX,
           target:=CX,
           mapping:=fproj,
           properties:=[["characteristic", 0],["type","chainMap"]]));
secondproj:= Objectify(HapChainMap,
       rec(
           source:=CXX,
           target:=CX,
           mapping:=sproj,
           properties:=[["characteristic", 0],["type","chainMap"]]));

diag!.firstProjection:=firstproj;
diag!.secondProjection:=secondproj;
return diag;
end);
#######################################
#######################################

#######################################
InstallGlobalFunction(RegularCWSimplex,
function(arg)
local n,Y,k;
n:=arg[1];
Y:=MaximalSimplicesToSimplicialComplex([[1..n+1]]);
Y:=RegularCWComplex(Y);
if Length(arg)=2 then return Y; fi;
Y!.directed:=[];
for k in [1..Y!.nrCells(1)] do
Add(Y!.directed,SSortedList(Y!.boundaries[2][k]{[2,3]}));
od;
return Y;
end);
#######################################

#######################################
InstallGlobalFunction(RegularCWCube,
function(n)
local I,Y;

if n=0 then
I:=[ [[1,0]], [] ];;
return RegularCWComplex(I);;
fi;


I:=[ [[1,0],[1,0]], [[2,1,2]], [] ];;
I:=RegularCWComplex(I);;
Y:=DirectProduct(List([1..n],i->I));
return Y;
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
