
########################################
InstallGlobalFunction(PureCubicalKnot,
function(arg)
local L,dim,K,col,row,i,j,x,y,yy,F,bool,name;

if Length(arg)=1 then
if not IsList(arg[1]) then
Print("Argument should be a pair of positive integers or a list describing a cubical arc presentation.\n");
return fail;
fi; 
L:=arg[1];
name:="unknown";
fi;
if Length(arg)=2 then
if not (IsInt(arg[1]) and IsInt(arg[2])) then
Print("Argument should be a pair of positive integers or a list describing a cubical arc presentation.\n");
return fail;
fi;
if arg[1]<3 then 
Print("There are no knots on ",arg[1]," crossings.\n");
return fail;
fi;
if arg[1]>9 then 
Print("So far only prime knots with fewer than 10 crossings have been stored in HAP.\n");
return fail;
fi;

L:=HAP_Knots[arg[1]-2][arg[2]];
name:=Concatenation("prime knot ",String(arg[2])," with ",String(arg[1])," crossings");
fi;

L:=Reversed(L);
L:=List(L,x->SortedList(x));


dim:=3*Length(L);
K:=[NullMat(dim,dim),NullMat(dim,dim),NullMat(dim,dim)];

for i in [1..Length(L)] do
row:=3*i-2;
for j in [3*L[i][1]-2..3*L[i][2]-2] do
K[3][j][row]:=1;
od;
K[2][3*L[i][1]-2][row]:=1;
K[2][3*L[i][2]-2][row]:=1;
od; 

for i in [1..Length(L)] do
col:=3*i-2;
F:=Filtered(L,a->i in a);
x:=Position(L,F[1]);
yy:=Positions(L,F[2]);
y:=yy[Length(yy)];
for j in [3*x-2..3*y-2] do
K[1][col][j]:=1;
od;
od;

K:=FrameArray(K);

K:=Objectify(HapPureCubicalLink,
           rec(
           binaryArray:=K,
	   name:=name,
           properties:=[
           ["dimension",ArrayDimension(K)],
           ["arraySize",ArrayDimensions(K)]]
           ));
bool:=PathComponentOfPureCubicalComplex(K,0)=1;
Add(K!.properties,["knot",bool]);

return K;

end);
########################################

########################################
InstallGlobalFunction(KnotGroup,
function(K)
local C,Y,G;

C:=ComplementOfPureCubicalComplex(K);
Y:=CubicalComplexToRegularCWComplex(C);
CriticalCellsOfRegularCWComplex(Y);
G:=FundamentalGroup(Y);
return G;

end);
########################################

########################################
InstallGlobalFunction(ViewPureCubicalKnot,
function(L)
local A, B, B1, B2, C, D, row, col, i,ii,jj, j;

B:=L!.binaryArray[2];
row:=10*Length(B);
col:=10*Length(B[1]);
A:=NullMat(row,col);
C:=NullMat(row,col);
D:=NullMat(row,col);

for i in [1..Length(B)] do
for j in [1..Length(B[1])] do
if B[i][j]>0 then 
for ii in [0..9] do
for jj in [0..9] do
A[10*i+ii][10*j+jj]:=1;
od;od;
fi;
od;
od;

B:=L!.binaryArray[4];
for i in [1..Length(B)] do
for j in [1..Length(B[1])] do
if B[i][j]>0 then
for ii in [0..9] do
for jj in [0..9] do
C[10*i+ii][10*j+jj]:=1;
od;od;
fi;
od;
od;

B1:=L!.binaryArray[2];
B1:=PureCubicalComplex(B1);
B2:=L!.binaryArray[4];
B2:=PureCubicalComplex(B2);
B2:=PureCubicalComplexIntersection(B1,B2);
B:=PureCubicalComplex(L!.binaryArray[3]);
B:=PureCubicalComplexDifference(B2,B);
B:=B!.binaryArray;

for i in [1..Length(B)] do
for j in [1..Length(B[1])] do
if B[i][j]>0  then
for ii in [0..10] do
for jj in [0..10] do
D[10*i+ii][10*j+jj]:=1;
od;
od;
fi;
od;
od;

A:=PureCubicalComplex(A);
C:=PureCubicalComplex(C);
D:=PureCubicalComplex(D);
C:=PureCubicalComplexDifference(C,
ThickenedPureCubicalComplex(ThickenedPureCubicalComplex(D)));
A:=PureCubicalComplexUnion(A,C);
ViewPureCubicalComplex(A);
end);
########################################

########################################
InstallGlobalFunction(KnotSum,
function(K,L)
local A,B,M,i,j,x,y,x1,x2,y1,y2;

if not (EvaluateProperty(K,"knot") and EvaluateProperty(K,"knot")) then
Print("The sum is only defined for knots.\n");
return fail;
fi;

A:=K!.binaryArray;
B:=L!.binaryArray;
M:=NullMat(Length(A[1])+Length(B[1])+1,Length(A[1][1])+Length(B[1][1])+1);
M:=List([1..5],i->StructuralCopy(M));

for i in [1..Length(A[2])] do
for j in [1..Length(A[2][1])] do
if A[2][i][j]=1 then M[2][i][j]:=1; fi;
if A[3][i][j]=1 then M[3][i][j]:=1; fi;
if A[4][i][j]=1 then M[4][i][j]:=1; fi;
od;
od;

x:=Length(A[1])+1;
y:=Length(A[1][1])+1;
for i in [1..Length(B[2])] do
for j in [1..Length(B[2][1])] do
if B[2][i][j]=1 then M[2][x+i][y+j]:=1; fi;
if B[3][i][j]=1 then M[3][x+i][y+j]:=1; fi;
if B[4][i][j]=1 then M[4][x+i][y+j]:=1; fi;
od;
od;

x:=Minimum(List(A[2], R->Position(Reversed(R),1)));
x1:=Length(A[2][1])-x+1;
y:=Minimum(List(B[2], R->Position(R,1)));
y1:=y+Length(A[1][1])+1;;

for i in [2..Length(M[2])-4] do
if M[2][x1][i]=1 and M[3][x1][i]=0 then M[2][x1][i]:=0;
else M[2][x1][i]:=1; fi;
od;
for i in [2..Length(M[2])-4] do
if M[2][y1][i]=1 and M[3][y1][i]=0 then M[2][y1][i]:=0;
else M[2][y1][i]:=1; fi;
od;

for j in [x1..y1] do
M[2][j][2]:=1;
M[2][j][Length(M[2])-4]:=1;
od;


M:=Objectify(HapPureCubicalLink,
           rec(
           binaryArray:=M,
           name:=Concatenation(K!.name," + ",L!.name),
           properties:=[
           ["dimension",ArrayDimension(M)],
           ["arraySize",ArrayDimensions(M)]]
           ));
Add(M!.properties,["knot",true]); 
return M;
end);
########################################

###################################################
InstallGlobalFunction(AlexanderMatrix,
function(G)
local R,M,i,j,b,w, F,xx,hom,GG, poly, ExponentSumWord;


R:=ResolutionAsphericalPresentation(G,2);
GG:=R!.group;
M:=[];
#F:=FreeGroup(1); xx:=F.1;
#hom:=GroupHomomorphismByImagesNC(GG,F,GeneratorsOfGroup(GG),
#List(GeneratorsOfGroup(GG),i->xx));

hom:=HAP_NqEpimorphismNilpotentQuotient(GG,1);
F:=Image(hom);
if not IsCyclic(F) then return fail; fi;
xx:=MinimalGeneratingSet(F)[1];

#################
ExponentSumWord:=function(g,xx)
local i;
for i in [0..10000] do
if xx^i=g then return i; fi;
if xx^-i=g then return -i; fi;
od;
end;
#################

#################
poly:=function(m);
if not IsZero(m) then return
SignInt(m)*LaurentPolynomialByCoefficients(FamilyObj(1),[SignInt(m)],m);
else return 1; fi;
end;
#################

for i in [1..R!.dimension(2)] do
M[i]:=[];
b:=R!.boundary(2,i);
for j in [1..R!.dimension(1)] do
w:=Filtered(b,x->AbsInt(x[1])=j);
w:=List(w,x->[SignInt(x[1]),Image(hom,R!.elts[x[2]])]);
w:=List(w,x->[x[1],ExponentSumWord(x[2],xx)]);
w:=List(w,x-> x[1]*poly(x[2]));
w:=Sum(w);
if w=0 then
w:=[[0],1];
else
w:=CoefficientsOfLaurentPolynomial(w);
fi;
w:=LaurentPolynomialByCoefficients(FamilyObj(1),w[1],0);

Add(M[i],StructuralCopy(w));
od;
od;

return M;
end);
#############################################

#############################################
InstallGlobalFunction(AlexanderPolynomial,
function(GG)
local M, row, col, Combs, C,G, dets, A, i, j;

if IsHapPureCubicalLink(GG) then
if PathComponentOfPureCubicalComplex(GG,0)=1 then
G:=KnotGroup(GG);
else
Print("The alexander polynomial is implemented only for knots.\,");
return fail;
fi;
fi;

if IsGroup(GG) then G:=GG; fi;

M:=AlexanderMatrix(G);
if M=fail then return fail; fi;

if Length(M) > Length(M[1]) then
M:=TransposedMat(M);
fi;

row:=Length(M);
col:=Length(M[1]);
Combs:=Combinations([1..col], row);
dets:=[];

for C in Combs do
A:=TransposedMat(M);
A:=A{C};
Add(dets,Determinant(A));
od;

return Gcd(dets);
end);
#############################################


############################################################
InstallGlobalFunction(ReadPDBfileAsPureCubicalComplex,
function(arg)
local  S,T,i,j,char,atoms, acids, cnt, 
 file, scl,  instr, f, avg,  v, A, x1,x2,y1,y2,z1,z2, L, M , B,b;

file:=arg[1];
scl:=2; char:="A"; 

if Length(arg)>1 then 
if IsRat(arg[2]) then scl:=arg[2];  else char:=arg[2]; fi; 
fi;

if Length(arg)>2 then 
if IsRat(arg[3]) then scl:=arg[3];  else char:=arg[3]; fi;
fi;



A:=[];
instr:=InputTextFile(file);
f:=ReadLine(instr);
while not f=fail do
f:=ReadLine(instr);
if f=fail then break; fi;
if Length(f)<4 then break; fi;
if f{[1,2,3,4]}="ATOM" and 
(f{[14,15]}="CA" or f{[15,16]}="CA") then  
if  f{[22]}=char then
v:=[f{[32,33,34,36,37,38]}, f{[40,41,42,44,45,46]}, f{[48,49,50,52,53,54]}]; 
Apply(v,EvalString); 
v[1]:=(scl*v[1]/1000);
v[2]:=(scl*v[2]/1000);
v[3]:=(scl*v[3]/1000);

Add(A,v); 
fi;
fi;
od;

atoms:=Length(A);

if atoms=0 then return fail; fi;


Print("Reading chain containing ",atoms," atoms.\n");

for j in [1..13] do
B:=[];
for i in [1..Length(A)-1] do
B[2*i-1]:=A[i];
v:=(A[i]+A[i+1])/2;
B[2*i]:=v;
od;
B[2*i+1]:=A[i+1];
A:=StructuralCopy(B);
B:=[];
od;

A:=List(A,a->[Int(a[1]),Int(a[2]),Int(a[3])]);

L:=List(A,a->a[1]);
x1:=Minimum(L);
x2:=Maximum(L);
L:=List(A,a->a[2]);
y1:=Minimum(L);
y2:=Maximum(L);
L:=List(A,a->a[3]);
z1:=Minimum(L);
z2:=Maximum(L);
x1:=x1-1;
y1:=y1-1;
z1:=z1-1;

M:=NullMat(y2-y1+1,z2-z1+1);
M:=List([1..x2-x1],i->StructuralCopy(M));

for f in A do
M[f[1]-x1][f[2]-y1][f[3]-z1]:=1;
od;

f:=A[1];
S:=[f[1]-x1,f[2]-y1,f[3]-z1];
f:=A[Length(A)];
T:=[f[1]-x1,f[2]-y1,f[3]-z1];
Add(M,M[1]*0);
Add(M,M[1]*0);
Add(M,M[1]*0);
for i in [S[1]+1..Length(M)] do
if M[i][S[2]][S[3]]=1 then 
Print("Failed to construct simple closed curve.\n"); return fail; fi;
M[i][S[2]][S[3]]:=1;
od;
for i in [T[1]+1..Length(M)] do
if M[i][T[2]][T[3]]=1 then 
Print("Failed to construct simple closed curve.\n"); return fail; fi;
M[i][T[2]][T[3]]:=1;
od;
for i in [1..Length(M[1])] do
for j in [1..Length(M[1][1])] do
M[Length(M)][i][j]:=1;
od;od;

M:=FrameArray(M);
M:=FrameArray(M);

M:=PureCubicalComplex(M);
ContractPureCubicalComplex(M);
M:=ThickenedPureCubicalComplex(M);
return M;


end);
#######################################################


#######################################################
InstallGlobalFunction(ProjectionOfPureCubicalComplex,
function(M)
local A;

A:=M!.binaryArray;
A:=-1*Sum(A);
A:=ArrayToPureCubicalComplex(A,0);

return A;
end);
#######################################################
