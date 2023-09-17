
################################################################
################################################################
Basins:=function(L)
local A, M,N, F, i,j,k, min,min2,max,pos,Q,B,P;
#L is a list of filtered pure cubical complexes

M:=Length(L[1]!.binaryArray);
N:=Length(L[1]!.binaryArray[1]);

for i in [1..Length(L)] do
if not (Length(L[i]!.binaryArray)=M
and Length(L[i]!.binaryArray[1])=N)
then Print("Filtered complexes are of differing dimensions.\n");
return fail; fi;
od;
A:=0*L[1]!.binaryArray;

for i in [1..M] do
for j in [1..N] do
F:=List(L,x->x!.filtration[i][j] * x!.binaryArray[i][j]);
max:=Maximum(F);
if max>0 then
min:=Minimum(Filtered(F,z->z>0));
pos:=Position(F,min);
A[i][j]:=pos;
fi;
od;od;

P:=PureCubicalComplex(0*L[1]!.binaryArray);
for k in [1..Length(L)] do
B:=0*L[1]!.binaryArray;
for i in [1..M] do
for j in [1..N] do
if A[i][j]=k then B[i][j]:=1; fi;
od;od;
Q:=PureCubicalComplex(B);
Q:=PureComplexDifference(Q,BoundaryOfPureCubicalComplex(Q));
P:=PureComplexUnion(P,Q);
od;

return P;
return PureCubicalComplex(A);;
end;
################################################################
################################################################


