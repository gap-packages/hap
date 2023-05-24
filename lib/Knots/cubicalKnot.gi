########################################
InstallGlobalFunction(PureCubicalLink,
function(X)
local L;

if IsString(X) then

if X{[1..4]}="Hopf" or X{[1..4]}="hopf" then
L:=PureCubicalKnot([[2,4],[1,3],[2,4],[1,3]]);
L!.name:="Hopf link";
return L;
fi;

fi;
return PureCubicalKnot(X);
end);
########################################


########################################
InstallGlobalFunction(PurePermutahedralKnot,
function(arg)
local K;

if Length(arg)=2 then
K:=PureCubicalKnot(arg[1],arg[2]);
fi;
if Length(arg)=1 then
K:=PureCubicalKnot(arg[1]);
fi;
K:= ContractedComplex(PurePermutahedralComplex(K!.binaryArray));
return PurePermutahedralComplex(FrameArray(K!.binaryArray));
end);
########################################

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
if arg[1]>11 then 
Print("So far only knots with fewer than 12 crossings have been stored in HAP.\n");
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
local C,Y,G, PureComp, PureCompToReg;

if IsHapPureCubicalComplex(K) then
C:=PureCubicalComplex(FrameArray(K!.binaryArray));
PureCompToReg:=CubicalComplexToRegularCWComplex;
PureComp:=PureCubicalComplex; fi;

if IsHapPurePermutahedralComplex(K) then
C:=PurePermutahedralComplex(FrameArray(K!.binaryArray));
PureCompToReg:=PermutahedralComplexToRegularCWComplex;
PureComp:=PurePermutahedralComplex; fi;

C:=ComplementOfPureComplex(C);
C:=ZigZagContractedPureComplex(C,"one iteration only");
#C:=ContractedComplex(C);
Y:=PureCompToReg(C);
Unbind(C);
Y:=ContractedComplex(Y);
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
local LF, RT, TL, BL, TR, BR, P, Q, x, 
A,B,M,y,x1,x2,y1,y2, bools,  s,t, i, j;

if not (EvaluateProperty(K,"knot") and EvaluateProperty(K,"knot")) then
Print("The sum is only defined for knots.\n");
return fail;
fi;

A:=K!.binaryArray;
B:=L!.binaryArray;

LF:=PositionProperty(Reversed(TransposedMat(A[2])), row-> not IsZero(row));
LF:=Length(A[2][1])+1-LF;
x:=List(A[2],row-> row[LF]);
TL:=Position(x, 1);
BL:=Position(Reversed(x),1);
BL:=Length(x)+1-BL;

RT:=PositionProperty(TransposedMat(B[2]), row-> not IsZero(row));
x:=List(B[2],row-> row[RT]);
TR:=Position(x, 1);
BR:=Position(Reversed(x),1);
BR:=Length(x)+1-BR;

RT:=RT+Length(A[1][1])+3;
TR:=TR+Length(A[1])-1;
BR:=BR+Length(A[1])-1;


P:=Length(A[2])+Length(B[2])-1;
Q:=Length(A[2][1])+Length(B[2][1])+3;
M:=NullMat( P , Q);
M:=List([1..5],i->StructuralCopy(M));

##########
##########INSERT BOTH KNOTS
for i in [1..Length(A[2])] do
for j in [1..Length(A[2][1])] do
if A[2][i][j]=1 then M[2][i][j]:=1; fi;
if A[3][i][j]=1 then M[3][i][j]:=1; fi;
if A[4][i][j]=1 then M[4][i][j]:=1; fi;
od;
od;

x:=Length(A[2])-1; 
y:=Length(A[2][1])+3; 
for i in [1..Length(B[2])] do
for j in [1..Length(B[2][1])] do
if B[2][i][j]=1 then M[2][x+i][y+j]:=1; fi;
if B[3][i][j]=1 then M[3][x+i][y+j]:=1; fi;
if B[4][i][j]=1 then M[4][x+i][y+j]:=1; fi;
od;
od;
##########BOTH KNOTS INSERTED
##########

##########
##########REMOVE FINAL/INITIAL COLUMNS
for i in [TL..BL] do
M[3][i][LF]:=0;
M[4][i][LF]:=0;
od;

for i in [TR..BR] do
M[3][i][RT]:=0;
M[4][i][RT]:=0;
od;
##########BOTH COLUMNS REMOVED
##########

##########
##########EXTEND TOP LEFT ROW
for i in [LF..RT] do
M[2][TL][i]:=1;
od;
##########
##########

##########
##########EXTEND BOTTOM ROW
for i in [LF..RT] do
M[2][BR][i]:=1;
od;
##########
##########

##########
##########ADD CENTRAL CONNECTIONS
M[3][TL][RT]:=1;
M[3][BL][LF]:=1;
M[3][TR][RT]:=1;
M[3][BR][LF]:=1;
##########
##########

##########
##########ADD LEFT COLUMN
for i in [BL..BR] do
M[4][i][LF]:=1;
od;
##########
##########

##########
##########ADD RIGHT COLUMN
for i in [TL..TR] do
M[4][i][RT]:=1;
od;
##########
##########


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

hom:=HAP_NqEpimorphismNilpotentQuotient(GG,1);
F:=Image(hom);
if not IsCyclic(F) then return fail; fi;
xx:=MinimalGeneratingSet(F)[1];

#################
ExponentSumWord:=function(g,xx)
local i;
for i in [0..1000000] do
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

Add(M[i],StructuralCopy(w));
od;
od;

return M;
end);
#############################################

#############################################
InstallGlobalFunction(AlexanderPolynomial,
function(GG)
local M, row, col, Cols, Rows, Combs, C,R,G, dets, gen, A, i, j;

gen:=LaurentPolynomialByCoefficients(FamilyObj(1),[1],1);

if IsHapPureCubicalLink(GG) then
if PathComponentOfPureCubicalComplex(GG,0)=1 then
G:=KnotGroup(GG);
else
Print("The alexander polynomial is implemented only for knots.\,");
return fail;
fi;
fi;

if IsGroup(GG) then G:=GG; fi;

if Length(GeneratorsOfGroup(G))=1 and Length(RelatorsOfFpGroup(G))=0
then return 1;
fi;

M:=(gen^0)*AlexanderMatrix(G);
if M=fail then return fail; fi;


row:=Length(M);
col:=Length(M[1]);
Cols:=Combinations([1..col], col-1);
Rows:=Combinations([1..row], col-1);
dets:=[];

for C in Cols do
for R in Rows do
A:=StructuralCopy(M);
A:=A{R};
A:=TransposedMat(A);
A:=A{C};
if Length(A)>0 then
Add(dets,Determinant(A));
fi;
od;
od;

Apply(dets,w->w*gen^-CoefficientsOfLaurentPolynomial(w)[2]);
A:=Gcd(dets);
return
A*Lcm(List(Flat(CoefficientsOfLaurentPolynomial(A)),x->DenominatorRat(x)));
end);
#############################################

############################################################
InstallGlobalFunction(ReadPDBfileAsPurePermutahedralComplex,
function(arg)
local   file, scl, char, A, B, M, f, g, v, instr, atoms, i, j,
        L, S, T, x1, x2, y1, y2, z1, z2, bool, x, ball;

file:=arg[1];
scl:=2; char:="A";

if Length(arg)>1 then
if IsRat(arg[2]) then scl:=arg[2];  else char:=arg[2]; fi;
fi;

if Length(arg)>2 then
if IsRat(arg[3]) then scl:=arg[3];  else char:=arg[3]; fi;
fi;

ball:=UnitPermutahedralBall(3);
Add(ball,[0,0,0]);
ball:=SSortedList(ball);


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

for j in [1..11+scl] do
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

A:=CubicalToPermutahedralArray(A);

A:=List(A,a->[Int(a[1]),Int(a[2]),Int(a[3])]);
for i in [2..Length(A)] do
if A[i]=A[i-1] then Unbind(A[i-1]); fi;
od;
A:=Filtered(A,a->IsBound(a));
#A:=DuplicateFreeList(A);
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

f:=A[1];
M[f[1]-x1][f[2]-y1][f[3]-z1]:=1;

for i in [2..Length(A)] do
f:=A[i];
M[f[1]-x1][f[2]-y1][f[3]-z1]:=1;

############################
#REPLACE THIS CODE BY A THICKENING AT THE END
#if not f-A[i-1] in ball then
#   g:=f-A[i-1];
#   bool:=true;
#   for x in ball do
#   #if f+x-g in ball then bool:=false; break; fi;
#   if g+x in ball then bool:=false; break; fi;
#   od;
#   if bool then Print("Warning: discontinuity in knot at ", f,"\n"); fi; 
#   g:=f+x;
#   M[g[1]-x1][g[2]-y1][g[3]-z1]:=1;
#fi;
############################
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

Unbind(A);

M:=FrameArray(M);
M:=FrameArray(M);

M:=PurePermutahedralComplex(M);
M:=ThickenedPureComplex(M);
ContractPureComplex(M);
return M;





end);
############################################################



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

for j in [1..11+scl] do
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

for i in [2..Length(A)] do
if A[i]=A[i-1] then Unbind(A[i-1]); fi;
od;
A:=Filtered(A,a->IsBound(a));


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

############################################################
############################################################
InstallGlobalFunction(DisplayPDBfile,
function(File)
local  S,T,i,j,char,atoms, acids, cnt, tmpdir, file,
 scl, instr, f, avg,  v, A, x1,x2,y1,y2,z1,z2, L, M , B,b,
  AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

scl:=2; char:="A";

A:=[];
instr:=InputTextFile(File);
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

A:=List(A,a->[Int(a[1]),Int(a[2]),Int(a[3])]);

M:=Maximum(List(Flat(A), x->AbsInt(x)));
M:=Concatenation( "(", String(M), "," , String(M), "," ,String(M), ")");
A:=List(A,a->Concatenation("(",String(a[1]),",",String(a[2]),",",String(a[3]),")"));

tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "tmp.asy" );

PrintTo(file, "import three;\n\n");
AppendTo(file, "size(500);\n\n");
AppendTo(file, "defaultpen(1.0);\n\n");

AppendTo(file,"path3 g=", A[1]);
for i in [2..Length(A)] do
AppendTo(file, "..",A[i]);
od;
AppendTo(file,";\n\n");

AppendTo(file,"draw(g);\n");

AppendTo(file,"dot(g,red);");

AppendTo( file , "draw(" , A[1] , "..",M,".." , A[Length(A)] , ", blue);" );


Exec( Concatenation( ASY_PATH, "-V ", file) );

RemoveFile(file);
file:=Filename(tmpdir,"");
RemoveFile(file);

end);
#############################################################
#############################################################


###########################################
###########################################
InstallGlobalFunction(ReflectedCubicalKnot,
function(K)
local A, R;

A:=1*K!.binaryArray;
A:=Reversed(A);
Apply(A,a->Reversed(TransposedMat(a)));
A:=ArcPresentation(PureCubicalComplex(A),"anything");
R:=PureCubicalKnot(A);
R!.name:=Concatenation("Reflected( ",K!.name," )");
return R;
end);
################################################
################################################

###########################################
###########################################
InstallGlobalFunction(ArcPresentation,
function(arg)
local K,A, P, L,i, row, m, n;

K:=arg[1];

if Length(arg)=1 and not IsHapPureCubicalLink(K) then
Print("This function applied only to pure cubical knots and links.\n");
return fail;
fi;


A:=K!.binaryArray[4];
A:=TransposedMat(A);

P:=[];


for i in [1..Length(A)] do
row:=A[i];
m:=Position(row, 1);
if not m=fail then
    n:=m;
    while not row[n+1]=0 do
    n:=n+1;
    od;
Add(P,[m,n]);
fi;
od;

P:= DuplicateFreeList(P);
L:=SSortedList(Flat(P));

Apply(P,a->[Position(L,a[1]),Position(L,a[2])]);

return Reversed(P);

end);
################################################
################################################

################################################
################################################
InstallGlobalFunction(NumberOfPrimeKnots,
function(n)
local L;

if n<1 then Print(0,"\n"); fi;
if n>19 then
Print("Only the numberof prime knots up to 19 crossings is currently avaiable.\n");
return fail;
fi;

L:=[0, 0, 1, 1, 2, 3, 7, 21, 49, 165, 552, 2176, 9988, 46972, 253293, 1388705, 8053393, 48266466, 294130458];
return L[n];
end);
################################################
################################################

#############################################
#############################################
InstallGlobalFunction(WirtingerGroup,
function(KK)
local crossing, K, A, i, j, start, toggle, colour, directions,
D,x, G, RELS, aa,bb,cc,dd,a,b,c,d,F,AA,tac,tbd,MX,mx,gens,m,M,
red,subs,subs2, CRS, sft;

##########################
if IsList(KK) then
return WirtingerGroup_gc(KK);
fi;
##########################

crossing:=-1;

K:=KK!.binaryArray;
A:=K[1]*0;;

#############################################
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
A[i][j]:=Maximum(K[2][i][j],K[4][i][j]);
od;od;

for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
if
A[i][j]=1 and A[i-1][j]=1 and A[i+1][j]=1 and A[i][j-1]=1 and A[i][j+1]=1
then A[i][j]:=crossing; crossing:=crossing-1;
fi;
od;od;
#############################################

#############################################
toggle:=true;
for i in [2..Length(A)-1] do
if toggle=false then break; fi;
for j in [2..Length(A[1])-1] do
if
A[i][j]<0 then start:=[i,j]; toggle:=false; break;
fi;
od;od;
#############################################

colour:=2;
directions:=[ [1,0], [-1,0], [0,1], [0,-1] ];
D:=[1,0];;

#############################################
while not D=false do

x:=start+D;

if A[x[1]][x[2]]>1  then  D:=false; fi;

if A[x[1]][x[2]]=1 then
A[x[1]][x[2]]:=colour; start:=x; fi;

if A[x[1]][x[2]]<0  then
start:=x; colour:=colour+1; fi;

if A[x[1]][x[2]]=0  then
D:=false;
for d in directions do
x:=start+d;
if A[x[1]][x[2]]=1 then D:=d; break; fi;
od;
fi;

od;
#####################################################

CRS:=[];
subs:=[];
#####################################################
for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
if
A[i][j]<0 then
m:=Minimum(A[i][j-1],A[i][j+1]);
M:=Maximum(A[i][j-1],A[i][j+1]);
CRS[-A[i][j]]:=[i,j];
A[i][j]:=m;
Add(subs,[m,M]);
fi;
od;od;
#############################################

subs2:=List(subs,x->x[2]);

#############################
#############################
red:=function(k)
local p;
if not k in subs2 then return k; fi;
p:=PositionProperty(subs,x->x[2]=k);
return red(subs[p][1]);
end;
#############################
#############################

#############################################
AA:=A*1;
for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
A[i][j]:=1*red(A[i][j]);
od;od;
#############################################

sft:=SSortedList(Flat(A));
sft:=Filtered(sft,j->not j=0);
sft:=SSortedList(sft);

for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
if A[i][j]>0 then A[i][j]:=Position(sft,A[i][j]); fi;
od;od;

mx:=Maximum(Flat(A));
MX:=Maximum(Flat(AA));
F:=FreeGroup(mx);
gens:=GeneratorsOfGroup(F);

RELS:=[];
for x in CRS do
a:=x+[0,1]; aa:=gens[A[a[1]][a[2]]];
b:=x+[-1,0]; bb:=gens[A[b[1]][b[2]]];
c:=x+[0,-1]; cc:=gens[A[c[1]][c[2]]];
d:=x+[1,0]; dd:=gens[A[d[1]][d[2]]];

tac:=( AA[a[1]][a[2]]>AA[c[1]][c[2]] and (not AA[a[1]][a[2]]=MX) ) or
(AA[a[1]][a[2]]=2 and AA[c[1]][c[2]]=MX) ;
tbd:=( AA[b[1]][b[2]]>AA[d[1]][d[2]]  and (not AA[b[1]][b[2]] =MX)) or
(AA[b[1]][b[2]]=2 and AA[d[1]][d[2]]=MX);


if (tac and tbd) then
Add(RELS,aa^-1*bb^-1*cc*dd);
fi;
if ((not tac) and tbd) then
Add(RELS,aa*bb^-1*cc^-1*dd);
fi;
if (tac and (not tbd)) then
Add(RELS,aa^-1*bb*cc*dd^-1);
fi;
if ((not tac) and (not tbd)) then
Add(RELS,aa*bb*cc^-1*dd^-1);
fi;

od;

G:=F/RELS;

return G;
end);
####################################################
####################################################

#############################################
#############################################
InstallGlobalFunction(GaussCodeOfPureCubicalKnot,
function(KK)
local crossing, start, colour, directions, code, startfn,
      cnt, orientations, K, A, D, d, i, j, x;

#Returns the Gauss code of a pure cubical complex representing an
#arc presentation of a knot or link.

crossing:=-1;
K:=KK!.binaryArray;
A:=K[1]*0;;

#############################################
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
   A[i][j]:=Maximum(K[2][i][j],K[4][i][j]);
od;od;

for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
if A[i][j]=1 and A[i-1][j]=1 and A[i+1][j]=1 and A[i][j-1]=1 and A[i][j+1]=1
   then
   A[i][j]:=crossing; crossing:=crossing-1;
fi;
od;od;

crossing:=Minimum(Flat(A));  #So crossing <=0
orientations:=List([1..AbsInt(crossing)],i->[]);
#############################################

#############################################
startfn:=function();
for i in [2..Length(A)-1] do
for j in [2..Length(A[1])-1] do
   if A[i][j]=1 and A[i][j-1]<0 then D:=[0,1]; return [i,j-1]; fi;
   if A[i][j]=1 and A[i][j+1]<0 then D:=[0,-1]; return [i,j+1]; fi;
   if A[i][j]=1 and A[i-1][j]<0 then D:=[1,0]; return [i-1,j]; fi;
   if A[i][j]=1 and A[i+1][j]<0 then D:=[-1,0]; return [i+1,j]; fi;
od;od;
end;
#############################################

colour:=2;
directions:=[ [1,0], [-1,0], [0,1], [0,-1] ];
code:=[];
cnt:=0;

#############################################
while 1 in Flat(A) do
   cnt:=cnt+1;
   code[cnt]:=[];
   start:=startfn();
   ##########################################
   while not D=false do
      x:=start+D;
      if A[x[1]][x[2]]>1 then D:=false; fi;
      if A[x[1]][x[2]]=1 then A[x[1]][x[2]]:=colour; start:=x; fi;
      if A[x[1]][x[2]]<0  then start:=x; colour:=colour+1;
         Add( orientations[AbsInt(A[x[1]][x[2]])], Filtered(D,b->not b=0));
         if K[2][x[1]-D[1]][x[2]-D[2]]=1 then Add(code[cnt], -A[x[1]][x[2]]);
         else
         Add(code[cnt], A[x[1]][x[2]]); fi;
      fi;

      if A[x[1]][x[2]]=0  then
         D:=false;
         for d in directions do
            x:=start+d;
            if A[x[1]][x[2]]=1 then D:=d; A[x[1]][x[2]]:=colour; start:=x;
               break;
            fi;
         od;
      fi;

   od;
   ##########################################
od;
#############################################

orientations:=List(orientations,x->Product(x));
return [code,orientations];

end);
#############################################
#############################################

########################################################
########################################################
InstallGlobalFunction(HAP_SimplifiedGaussCode,
function(gcc)
local gc, bool, i, j, L;

gc:=1*gcc;
bool:=true;
while bool do
bool:=false;
for i in [1..Length(gc[1])] do
for j in [1..Length(gc[1][i])-1] do
if AbsInt(gc[1][i][j]) = AbsInt(gc[1][i][j+1])  then
gc[2][AbsInt(gc[1][i][j])]:=0;
gc[1][i][j]:=0;
gc[1][i][j+1]:=0; bool:=true; fi;
od;
if AbsInt(gc[1][i][1]) = AbsInt(gc[1][i][Length(gc[1][i])]) then
gc[2][AbsInt(gc[1][i][1])]:=0;
gc[1][i][1]:=0;
gc[1][i][Length(gc[1][i])]:=0;
bool:=true;
fi;
gc[1][i]:=Filtered(gc[1][i],a->not a=0);
gc[2]:=Filtered(gc[2],a->not a=0);
od;
od;

L:=Flat(gc[1]);
L:=SSortedList(List(L,AbsInt));
for i in [1..Length(gc[1])] do
gc[1][i]:=List(gc[1][i],x->SignInt(x)*Position(L,AbsInt(x)));
od;

return gc;
end);
########################################################
########################################################

#############################################
#############################################
InstallGlobalFunction(WirtingerGroup_gc,
function(arg)
local A,GC, F, gens, rels, arcs, orientations, crossings, R, C, c, a, cr,i, x,y,z;

A:=arg[1];
A:=HAP_SimplifiedGaussCode(A);
#GC:=arg[1][1];
#orientations:=arg[1][2];
GC:=A[1];
orientations:=A[2];
crossings:=Flat(GC);
crossings:=List(crossings,AbsInt);
crossings:=SSortedList(crossings);

arcs:=[];

for C in GC do
cr:=Filtered([1..Length(C)],x->C[x]<0);
for i in [1..Length(cr)-1] do
a:=List( C{[cr[i]..cr[i+1]]} , AbsInt);
Add(arcs,a);
od;
x:=C{[cr[Length(cr)]..Length(C)]};
y:=C{[1..cr[1]]};
a:=List(Concatenation(x,y), AbsInt);
Add(arcs,a);
od;

F:=FreeGroup(Length(arcs));
gens:=GeneratorsOfGroup(F);
rels:=[];
for c in crossings do
R:=Filtered(arcs, x-> c in x);
for a in R do
if a[1]=c then x:=a; fi;
if a[Length(a)]=c then y:=a; fi;
if not (a[1]=c or a[Length(a)]=c) then z:=a; fi;
od;
x:=Position(arcs,x);
y:=Position(arcs,y);
z:=Position(arcs,z);
if orientations[c]=1 then
Add(rels,gens[x]^-1*gens[z]*gens[y]*gens[z]^-1);
else
Add(rels,gens[x]^-1*gens[z]^-1*gens[y]*gens[z]);
fi;
od;

return F/rels;
end);
#############################################
#############################################


