#(C) Graham Ellis


###################################################
###################################################
InstallGlobalFunction(ThickeningFiltration,
#function(M,n)
function(arg)
local M, n, F, D, T, C, k, A,step,s;

M:=arg[1];
n:=arg[2];
if Length(arg)>2 then step :=arg[3]; else step:=1; fi;

F:=PureCubicalComplex(M!.binaryArray);
F!.filtration:=StructuralCopy(F!.binaryArray);
T:=M;

for k in [2..n+1] do
   for s in [1..step] do
   T:=ThickenedPureCubicalComplex(T);
   od;
D:=PureCubicalComplexDifference(T,F);
D:=k*D!.binaryArray;
A:=F!.filtration + D;
F:=PureCubicalComplex(T!.binaryArray);
F!.filtration:=A;
od;

return Objectify(HapFilteredPureCubicalComplex, 
                 rec(binaryArray:=F!.binaryArray,
                     filtration:=F!.filtration,
                     properties:=F!.properties));
end);
###################################################
###################################################

###################################################
###################################################
InstallGlobalFunction(Dendrogram,
function(M)
local Betti, F, T, TT, k, phi, i;


if not IsBound(M!.filtrationLength) then
F:=Maximum(Flat(M!.filtration));
M!.filtrationLength:=F;
fi;
F:=M!.filtrationLength;


T:=FiltrationTerm(M,1);;
Betti:=[PathComponentOfPureCubicalComplex(T,0)];
phi:=[];

for k in [1..F-1] do
TT:=FiltrationTerm(M,k+1);
Add(Betti, PathComponentOfPureCubicalComplex(TT,0));

#phi is such that phi[k][i] = path component of TT
#that intersects with i-th path component of T.

phi[k]:=[];

###
for i in [1..Betti[k]] do

Add(phi[k],
ArrayValue(TT!.pathCompBinList, T!.pathReps[i])-1
);

od;
###

T:=TT;
od;

return phi;

end);
###################################################
###################################################

###################################################
###################################################
InstallGlobalFunction(FiltrationTerm,
function(M,t)
local
        B,
        CART, dim,dim1,dims,
        Opp,
        ArrayValueDim,
        ArrayValueDim1,
        x,z;


#############################################
if not IsHapFilteredPureCubicalComplex(M) then
Print("This function must be applied to a filtered pure cubical complex.\n");
return fail;
fi;
#############################################

dim:=Dimension(M);
dim1:=dim-1;
dims:=EvaluateProperty(M,"arraySize");
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
B:=0*M!.binaryArray;
CART:=Cartesian(List([1..dim],a->[1..dims[a]]));

#######################
Opp:=function(y)
local z;
z:=ArrayValueDim1(B,y{[2..Length(y)]});
if ArrayValueDim(M!.filtration,y)<=t and
not ArrayValueDim(M!.filtration,y)=0 then
z[y[1]]:=1;
fi;
end;
########################

for x in CART do
Opp(x);
od;

return Objectify(HapPureCubicalComplex,
           rec(
           binaryArray:=B,
           properties:=[
           ["dimension",Dimension(M)],
           ["arraySize",dims]]
           ));

end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(DendrogramDisplay,
function(M)
local D, V, G,  t, i, shft1, shft2, tmpDir, Loggeddot, Loggedpng;


tmpDir:=DirectoryTemporary();
Loggeddot:=Filename(tmpDir,"Logged.dot");
Loggedpng:=Filename(tmpDir,"Logged.png");


if IsHapFilteredPureCubicalComplex(M) then
D:=Dendrogram(M);
fi;
if IsList(M) then D:=M; fi;
if not ( IsList(M) or IsHapFilteredPureCubicalComplex(M) ) then
Print("The function must be applied to a filtered pure cubical complex or a list representing a denrogram.\n");
return fail;
fi;

V:=List(D,Length);
V:=Sum(V);

shft1:=0;



PrintTo(Loggeddot,"digraph { \n node [shape=circle, style=filled, color=black] \n ");

AppendTo(Loggeddot, "rankdir=LR;\n\n");
AppendTo(Loggeddot, "node [style=filled,shape=point]\n\n");

for t in [1..Length(D)-1] do
shft2:=shft1+Length(D[t]);
for i in [1..Length(D[t])] do
AppendTo(Loggeddot, shft1+i, "->", shft2+D[t][i]," [color=black,arrowhead=none];\n");
od;
shft1:=shft2;
od;

AppendTo(Loggeddot,"}\n");

Exec(Concatenation("dot -Tpng ",Loggeddot," -o ",Loggedpng));
Exec(Concatenation("rm ",Loggeddot));
Exec(Concatenation(DISPLAY_PATH,Loggedpng));
Exec(Concatenation("rm ",Loggedpng));



end);
#################################################################
#################################################################



#################################################################
#################################################################
InstallGlobalFunction(ReadImageAsFilteredCubicalComplex,
function(file,N)
local A, B, C, F, i, j;

A:=ReadImageAsPureCubicalComplex(file,"matrix");
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
A[i][j]:=Int(A[i][j]/N);
od;od;
F:=Maximum(Flat(A));

#C:=List([1..Length(A[1])], i->F);
#C:=List([1..Length(A)], i->StructuralCopy(C));
#A:=C-A;

C:=A*0;
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
if A[i][j]>0 then C[i][j]:=1;fi;
od;od;


return Objectify(HapFilteredPureCubicalComplex,
                 rec(binaryArray:=C,
                     filtration:=A,
                     properties:=[["dimension",2], ["arraySize",ArrayDimensions(A) ]]));



end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(ComplementOfFilteredCubicalComplex,
function(M)
local  F, B, A, T, TT,      D,   k;


F:=Maximum(Flat(M!.filtration));
T:=FiltrationTerm(M,1);
B:=ComplementOfPureCubicalComplex(T);
B:=StructuralCopy(B!.binaryArray);
A:=0*B;

T:=PureCubicalComplex(0*B);


for k in Reversed([1..F]) do
TT:=FiltrationTerm(M, k);
TT:=ComplementOfPureCubicalComplex(TT);


D:=PureCubicalComplexDifference(TT,T);
D:=(F+1-k)*D!.binaryArray;
A:=A + D;
T:=TT;
od;

return Objectify(HapFilteredPureCubicalComplex,
                 rec(binaryArray:=B,
                     filtration:=A,
                     properties:=T!.properties));



end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(DendrogramToPersistenceMat,
function(D)
local B, phi, s, t;

B:=NullMat(Length(D),Length(D));


###########################
phi:=function(s,t)   #only for t>s
local L, x;

L:=StructuralCopy(D[s]);

for x in [s+1..t-1] do
Apply(L, i->D[x][i]);
od;

return Size(SSortedList(L));

end;
###########################

for s in [1..Length(D)] do
B[s][s]:=Length(D[s]);
od;

for s in [1..Length(D)] do
for t in [s+1..Length(D)] do
B[s][t]:=phi(s,t);
od;
od;



return B;
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(FilteredPureCubicalComplexToCubicalComplex,
function(M)
# Converts a filtered pure cubical complex to a filtered cubical complex
local A,F,x,i,dim,dim1,dims,ball,b,dimsM,ArrayValueDim, ArrayValueDim1,
ArrayAssignDim,ArrayIt,dimsSet,Fun, infty;

dim:=Dimension(M);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim-1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
dimsM:=EvaluateProperty(M,"arraySize");
dims:=List(EvaluateProperty(M,"arraySize"),n->2*n+1);
infty:=Maximum(Flat(M!.filtration));
A:=List([1..dims[1]],a->0);
F:=List([1..dims[1]],a->infty);
for i in [2..Dimension(M)] do
A:=List([1..dims[i]],a->StructuralCopy(A));
F:=List([1..dims[i]],a->StructuralCopy(F));
od;

ball:=Cartesian(List([1..dim],i->[-1,0,1]));
dimsSet:=List([1..dim],x->[1..dimsM[x]]);

#########################
Fun:=function(i) local m;
if ArrayValueDim(M!.binaryArray,i)=1 then
for b in ball do
x:=2*i+b;
ArrayAssignDim(A,x,1);
m:=ArrayValueDim(M!.filtration,i);
m:=Minimum(m,ArrayValueDim(F,x));
ArrayAssignDim(A,x,m );
od;
fi;
end;
#########################
ArrayIt(dimsSet,Fun);
return Objectify(HapFilteredCubicalComplex,
           rec(
           binaryArray:=A,
           filtration:=F,
           properties:=[
           ["dimension",dim],
           ["arraySize",dims]]
           ));

end);
#################################################################
#################################################################



