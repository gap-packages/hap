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

phi[F]:=[1..Betti[F]]*0;

return phi;

end);
###################################################
###################################################

###################################################
###################################################
InstallGlobalFunction(FiltrationTerm,
function(M,t)
local
        B, obj,
        CART, dim,dim1,dims,
        Opp,
        ArrayValueDim,
        ArrayValueDim1,
        dimSet, ArrayIt,
        x,z;


#############################################
if not (IsHapFilteredPureCubicalComplex(M) or 
IsHapFilteredCubicalComplex(M)) then
Print("This function must be applied to a filtered cubical or pure cubical complex.\n");
return fail;
fi;
#############################################

if IsHapFilteredPureCubicalComplex(M) then
obj:=HapPureCubicalComplex;
else
obj:=HapCubicalComplex;
fi;

dim:=Dimension(M);
dim1:=dim-1;
dims:=EvaluateProperty(M,"arraySize");
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
B:=0*M!.binaryArray;
#CART:=Cartesian(List([1..dim],a->[1..dims[a]]));

#######################
Opp:=function(y)
local z;
z:=ArrayValueDim1(B,y{[2..Length(y)]});
if ArrayValueDim(M!.filtration,y)<=t and
not ArrayValueDim(M!.binaryArray,y)=0 then
z[y[1]]:=1;
fi;
end;
########################

#for x in CART do
#Opp(x);
#od;

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);
ArrayIt(dimSet,Opp);

return Objectify(obj,
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
Print("The function must be applied to a filtered pure cubical complex or a list representing a dendrogram.\n");
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
local A, B, C, F, i, j, D;
D:=Int(3*255/N);

A:=ReadImageAsPureCubicalComplex(file,"matrix");
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
A[i][j]:=1+Int(A[i][j]/D);
od;od;
#F:=Maximum(Flat(A));

C:=A*0;
C:=C+1;

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


if IsBound(M!.filtrationLength) then
F:=M!.filtrationLength;
else
F:=Maximum(Flat(M!.filtration));
fi;

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
ArrayAssignDim(F,x,m );
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


###############################################################
###############################################################
InstallGlobalFunction(ContractedFilteredPureCubicalComplex,
function(M)
local
        LM, F, i, flen, A,
        dim, dim1, dims, ArrayValueDim, ArrayValueDim1, B, CART,
        Opp, x,t, dimSet, ArrayIt;

LM:=[];
flen:=Maximum(Flat(M!.filtration));

for i in [1..flen] do
Add(LM,FiltrationTerm(M,i));
od;

F:=HomotopyEquivalentMinimalPureCubicalSubcomplex;

###############################
LM[1]:=ContractedComplex(LM[1]);
for i in [2..Length(LM)] do
LM[i]:=F(LM[i],LM[i-1]);
od;
###############################



dim:=Dimension(M);
dim1:=dim-1;
dims:=EvaluateProperty(M,"arraySize");
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
#CART:=Cartesian(List([1..dim],a->[1..dims[a]]));

#######################
Opp:=function(y)
local z;
z:=ArrayValueDim1(B,y{[2..Length(y)]});
if ArrayValueDim(LM[t]!.binaryArray,y)=1
and ArrayValueDim(LM[t-1]!.binaryArray,y)=0
 then
#z[y[1]]:=z[y[1]]+1;
z[y[1]]:=t;
fi;
end;
########################

B:=StructuralCopy(LM[1]!.binaryArray);

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);

for t in [2..Length(LM)] do
#for x in CART do
#Opp(x);
#od;
ArrayIt(dimSet,Opp);
od;

return Objectify(HapFilteredPureCubicalComplex,
                 rec(binaryArray:=LM[flen]!.binaryArray,
                     filtration:=B,
                     properties:=M!.properties));


end);
###############################################################
###############################################################

###############################################################
###############################################################
InstallGlobalFunction(ZigZagContractedFilteredPureCubicalComplex,
function(M)
local
        LM, F, i, flen, A, G,
        dim, dim1, dims, ArrayValueDim, ArrayValueDim1, B, CART,
        Opp, x,t, CM, sz1, sz2, maxB, dimSet, ArrayIt;

LM:=[];
flen:=Maximum(Flat(M!.filtration));

for i in [1..flen] do
Add(LM,FiltrationTerm(M,i));
od;

F:=HomotopyEquivalentMinimalPureCubicalSubcomplex;
G:=HomotopyEquivalentMaximalPureCubicalSubcomplex;
sz1:=infinity;
sz2:=Size(LM[flen]);

while sz2<sz1 do

sz1:=sz2;
LM[flen]:=F(LM[flen],LM[flen-1]);

###############################
CM:=CropPureCubicalComplex(LM[flen]);
for i in Reversed([2..Length(LM)]) do
LM[i-1]:=G(LM[i],LM[i-1]);
od;
###############################

###############################
LM[1]:=ContractedComplex(LM[1]);
for i in [2..Length(LM)] do
LM[i]:=F(LM[i],LM[i-1]);
od;
###############################
sz2:=Size(LM[flen]);

od;
dim:=Dimension(M);
dim1:=dim-1;
dims:=EvaluateProperty(M,"arraySize");
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
#CART:=Cartesian(List([1..dim],a->[1..dims[a]]));



B:=StructuralCopy(LM[1]!.binaryArray);

#######################
Opp:=function(y)
local z;
z:=ArrayValueDim1(B,y{[2..Length(y)]});
if ArrayValueDim(LM[t]!.binaryArray,y)=1
and ArrayValueDim(LM[t-1]!.binaryArray,y)=0
 then
z[y[1]]:=t;
fi;
end;
########################

dimSet:=List([1..dim],x->[1..dims[x]]);
ArrayIt:=ArrayIterate(dim);

for t in [2..Length(LM)] do
#for x in CART do
#Opp(x);
#od;
ArrayIt(dimSet,Opp);
od;


return Objectify(HapFilteredPureCubicalComplex,
                 rec(binaryArray:=LM[flen]!.binaryArray,
                     filtration:=B,
                     filtrationLength:=flen,
                     properties:=
                     [["dimension",Dimension(M)],
                      ["arraySize",ArrayDimensions(M!.binaryArray)]     ]
                      ));


end);
###############################################################
###############################################################

##########################################################
##########################################################
InstallGlobalFunction(ConcentricallyFilteredPureCubicalComplex,
function(M,N)
local CentreOfGravity, Radius, rad, cen, F, B, x,y,z;

##############################
if not IsHapPureCubicalComplex(M) then
Print("Function must be applied to a pure cubical complex.\n");
return fail;
fi;

if not Dimension(M)=3 then
Print("At present this function is only implemented for 3-dimensional pure cubical complexes.\n");
return fail;
fi;
############################


#########################################
CentreOfGravity:=function(M)
local B,x,y,z,V;

V:=[];
B:=M!.binaryArray;

for x in [1..Length(B)] do
for y in [1..Length(B[1])] do
for z in [1..Length(B[1][1])] do
if B[x][y][z]=1 then Add(V,[x,y,z]); fi;
od;od;od;

V:= (1/Length(V))*Sum(V);
V:=[Int(V[1]),Int(V[2]),Int(V[3])];
return V;
end;
#########################################

#########################################
Radius:=function(M)
local B;

B:=M!.binaryArray;
return Maximum(ArrayDimensions(B))/2;
end;
#########################################

rad:=Radius(M);
cen:=CentreOfGravity(M);
F:=M!.binaryArray*0;
B:=M!.binaryArray;

for x in [1..Length(B)] do
for y in [1..Length(B[1])] do
for z in [1..Length(B[1][1])] do
if B[x][y][z]=1 then
F[x][y][z]:= 1+Int(N*EuclideanApproximatedMetric(cen, [x,y,z])/rad) ;  fi;
od;od;od;

return   Objectify(HapFilteredPureCubicalComplex,
                 rec(binaryArray:=M!.binaryArray,
                     filtration:=F,
                     properties:=M!.properties));

end);
########################################################
########################################################



