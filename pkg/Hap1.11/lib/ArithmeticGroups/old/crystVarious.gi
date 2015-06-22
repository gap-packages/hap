InstallGlobalFunction(IsIntList,
function(list)
local i;
for i in list do
if not IsInt(i) then return false;fi;
od;
return true;
end);
#############################################

InstallGlobalFunction(VectorToCrystMatrix,      
function(v)
local M,n;
v:=Flat(v);
n:=Length(v);
M:=IdentityMat(n+1);
Add(v,1);
Remove(M);
Add(M,v);
return M;
end);
###############################################

InstallGlobalFunction(CrystTranslationMatrixToVector,
function(g)
local n,v;
n:=Length(g);
v:=g[n];
v:=Flat(v);
Remove(v);
return v;
end);
###############################################

InstallGlobalFunction(TranslationSubGroup,
function(G)
local B,SbGrp,trsltbasis;
B:=TranslationBasis(G);
if not IsBound(G!.TranslationBasis) then return false;fi;
B:=G!.TranslationBasis;
trsltbasis:=List(B,w->VectorToCrystMatrix(Flat(w)));
SbGrp:=Group(trsltbasis);
SbGrp!.TranslationBasis:=B;
SetIsCrystTranslationSubGroup(SbGrp,true);
return SbGrp;
end);
###############################################

InstallOtherMethod(\in, 
		"for TranslationSubGroup of a CrystGroup",
	      [IsMatrix,IsCrystTranslationSubGroup],
function(g,G)
local B,v,n;
n:=DimensionSquareMat(g)-1;
if not LinearPartOfAffineMatOnRight(g)=IdentityMat(n) then return false;fi;
B:=G!.TranslationBasis;
v:=CrystTranslationMatrixToVector(g);

return IsIntList(v*TransposedMat(B)^-1);

end);
###############################################

InstallGlobalFunction(IsCrystSameOrbit,
function(arg)
local G,T,H,u,v,B,x,w;
G:=arg[1];
if Length(arg)=3 then
H:=TranslationSubGroup(G);
T:=RightTransversal(G,H);
B:=H!.TranslationBasis;
u:=arg[2];
v:=arg[3];
else
B:=arg[2];
T:=arg[3];
u:=arg[4];
v:=arg[5];
fi;
u:=Flat(u);
v:=Flat(v);
Add(u,1);
Add(v,1);
for x in T do
w:=u*x-v;
w:=Flat(w);
Remove(w);
if IsIntList(w*TransposedMat(B)^-1) then return x*VectorToCrystMatrix(w)^-1;fi;
od;
return false;
end);
#################################################
InstallGlobalFunction(CombinationDisjointSets,
function(arg)
local b,list,n1,i,g,h;
g:=arg[1];
if g=[] then return [[]];fi;
n1:=g[1];
h:=g{[2..Length(g)]};
list:=[];
b:=CombinationDisjointSets(h);
for i in [0..(n1-1)] do
Append(list,List(b,w->AddFirst(w,i)));
od;
return list;
end);
################################################
InstallGlobalFunction(AddFirst,
function(list,g)                 # add g in the first position in list
local w;
w:=[g];
Append(w,list);
return w;
end);
###############################################
InstallGlobalFunction(IsCrystSufficientLattice,
function(B,S)
local v,x,w,B1,c,i,A,Origin,S1,S2;
Origin:=0*B[1];
Origin:=Flat(Origin);
Add(Origin,1);
A:=StructuralCopy(B);
v:=Sum(A)/2;
v:=Flat(v);
Add(v,1);
for x in S do
w:=Flat(v*x-v);
if not IsIntList(w*A^-1) then return false;fi;
od;
for i in [1..Length(A)] do
A[i]:=Flat(A[i]);
Add(A[i],1);
od;
for x in S do

B1:=List(A,w->w*x);
c:=Flat(Origin*x-Origin);
B1:=List(B1,w->w-c);
for x in B1 do
Remove(x);
od;
#S1:=Set(A);
S2:=Set(B1);
for x in B do
if not ((x in S2) or (-x in S2)) then return false;fi;
od;
od;
return true;
end);
###############################################################
InstallGlobalFunction(CrystFinitePartOfMatrix,
function(g)
local x,w,i;
w:=[];
for i in [1..(Length(g)-1)] do
x:=Flat(g[i]);
Remove(x);
Add(w,x);
od;
return w;
end);

###############################################################
InstallGlobalFunction(ResolutionBoundaryOfWordOnRight,
function(R,n,W)
local x, DW, Boundary, Dimension,Elts,pos, ans,H;

Dimension:=R!.dimension;
Boundary:=R!.boundary;
Elts:=R!.elts;
DW:=[];

for x in W do
#x:=[R!.action(n,x[1],x[2])*x[1],x[2]];
ans:=Boundary(n,x[1]);
ans:=List(ans, a->[a[1],Elts[a[2]]]);
ans:=List(ans, a->[a[1],a[2]*Elts[x[2]]]);
Append(DW,ans);
od;

DW:= AlgebraicReduction(DW);
for x in DW do
if not x[2] in Elts then Add(Elts,x[2]);fi;
od;
DW:=List(DW,x->[x[1],Position(Elts,x[2])]);
DW:=List(DW,x->[R!.Sign(n-1,x[1],x[2])*x[1],x[2]]);
for x in DW do
H:=R!.stabilizer(n-1,AbsInt(x[1]));
pos:=Position(R!.elts,CanonicalRightCountableCosetElement(H,R!.elts[x[2]]));
if pos=fail then Add(R!.elts,CanonicalRightCountableCosetElement(H,R!.elts[x[2]])); 
x[2]:=Length(R!.elts);
else x[2]:=pos;
fi;
od;
DW:=List(DW,x->[R!.Sign(n-1,x[1],x[2])*x[1],x[2]]);
DW:= AlgebraicReduction(DW);
return DW;

end);
###############################################################
InstallGlobalFunction(CrystCubicalTiling,
function(n)
local combin,x,w,Til,i;
combin:=Combinations([1..n],2);
Til:=[];
Add(Til,IdentityMat(n));
for i in [1..Length(combin)] do
w:=combin[i];
x:=IdentityMat(n);
x[w[1]][w[1]]:=-1/2;
x[w[1]][w[2]]:=-1/2;
x[w[2]][w[1]]:=1;
x[w[2]][w[2]]:=-1;
Add(Til,x);
od;
return Til;
end);

################################################################
InstallGlobalFunction(AverageInnerProduct,
function(G,u,v)
local i,sum,n,Elts;
n:=Order(G);
Elts:=Elements(G);
sum:=0;
for i in [1..n] do
sum:=sum+(u*Elts[i])*(v*Elts[i]);
od;
sum:=sum/n;
return sum;
end);
#########
InstallGlobalFunction(OrthogonalizeBasisByAverageInnerProduct,
function(B,G)
local Project,i,j,A,n;
A:=StructuralCopy(B);
n:=Length(B);
if not RankMat(B)=Length(B) then 
Print("Input is not a basis");
return fail;
fi;
###
Project:=function(u,v)    #This operator projects the vector v orthogonally 
                          #onto the line spanned by vector u
return (AverageInnerProduct(G,u,v)/AverageInnerProduct(G,u,u))*u;
end;
###
for i in [2..n] do
for j in [1..(i-1)] do
B[i]:=B[i]-Project(B[j],A[i]);
od;od;
for i in [1..n] do
B[i]:=(1/Sqrt(AverageInnerProduct(G,B[i],B[i])))*B[i];
od;
return B;



end);
###################################
InstallGlobalFunction(CrystMatrix,
function(M)
local i,n,x,N;
N:=StructuralCopy(M);
if IsMatrix(N) then
n:=Length(N);
x:=0*N[1];
Add(x,1);
for i in [1..n] do
Add(N[i],0);
od;
Add(N,x);
else
N:=VectorToCrystMatrix(N);
fi;
return N;
end);