

####################################
####################################
InstallGlobalFunction(LowDimensionalCupProduct,
function(Y)
local W,G,R, HG2, HG1, H2,cup, B, i, j, k, A,P, newcup;

W:=RegularCWComplex(1*Y!.boundaries);
if not Length(CocriticalCellsOfRegularCWComplex(W,0))=1 then
Print("At present this is implemented for connected spaces only.\n");
return fail; fi;

G:=FundamentalGroup(Y);;
R:=ResolutionAsphericalPresentation(G,3);
HG2:= Cohomology(HomToIntegers(R),2);;
H2:=Cohomology(Y,2);
cup:=CupProduct(G);
if Length(HG2)=Length(H2) then return cup; fi;


HG1:=Cohomology(HomToIntegers(R),1);;

P:=[];
B:=IdentityMat(Length(HG1));
for i in [1..Length(B)] do
for j in [i..Length(B)] do
Add(P,cup(B[i],B[j]));
od;od;
P:=1*TransposedMat(P);

k:=1;
A:=[];
while k<=Length(HG2)-Length(H2)+1 do
if IsZero(P[k]) then Add(A,k);  fi;
k:=k+1;
od;

B:=[1..Length(H2)];
B:=Difference(B,A);
newcup:=function(x,y);
return cup(x,y){B};
end;

return newcup;
end);
####################################
####################################

####################################
####################################
InstallGlobalFunction(CupProductMatrix,
function(YY,p,q)
local Y, A, Hp,Hq, Bp, Bq, cup, i,j;

Y:=YY;
Hp:=Cohomology(Y,p);;
Hq:=Cohomology(Y,q);;
Bp:=IdentityMat(Length(Hp));
Bq:=IdentityMat(Length(Hq));
Y:=ContractedComplex(Y);
#Y:=SimplifiedComplex(Y);
cup:=CupProduct(Y);

A:=[];
for i in [1..Length(Bp)] do
Add(A,[]);
for j in [1..Length(Bq)] do
A[i][j]:=cup(p,q,Bp[i],Bq[j]);
od;od;

return A;
end);
###################################
###################################
