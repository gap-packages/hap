#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IntegralRingGenerators,
function(R,n)
local
		Gens, Cups, cupProd, CB, CohomologyBasis, TR, 
		BasisP, BasisQ, SmithRecord,
		i, p, q, u, v, ln;

TR:=HomToIntegers(R);
if Length(Cohomology(TR,n))=0 then return []; fi;

#####################################################################
CohomologyBasis:=function(Torsion)
local i, v, Basis;
Basis:=[];
for i in [1..Length(Torsion)] do
v:=List([1..Length(Torsion)], j->0);
v[i]:=Torsion[i];
Append(Basis, [v]);
od;
return Basis;
end;
#####################################################################

Cups:=CohomologyBasis(Cohomology(TR,n));

CB:=[];
for p in [1..n] do
CB[p]:=CR_CocyclesAndCoboundaries(R,p,true);
od;

for p in [1..(n-1)] do
q:=n-p;
ln:=Length(Cohomology(TR,p));
BasisP:=CohomologyBasis(List([1..ln],i->1));
ln:=Length(Cohomology(TR,q));
BasisQ:=CohomologyBasis(List([1..ln],i->1));

for u in BasisP do
for v in BasisQ do
Append(Cups,
[ IntegralCupProduct(R,u,v,p,q,CB[p],CB[q],CB[n]) ]);

od;
od;

od;

SmithRecord:=SmithNormalFormIntegerMatTransforms(Cups);

Gens:=[];
for i in [1..Length(Cohomology(TR,n))] do
v:= SmithRecord.normal[i];
if v[i]>1 or v[i]=0 then 
v[i]:=1;
u:=TransposedMat(SmithRecord.coltrans)*v;
Append(Gens,[u]); fi;

od;

return Gens;
end);
#####################################################################

IntegralCohomologyGenerators:=IntegralRingGenerators;
MakeReadOnlyGlobal("IntegralCohomologyGenerators");

