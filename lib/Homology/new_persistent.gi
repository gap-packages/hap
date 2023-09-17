
###################################################
###################################################
InstallGlobalFunction(PersistentBettiNumbersViaContractions,
function(FF,n,p)
local F, W, L, E, chnmap, cwmap, map, ln, k, f, homs;

#F:=ContractedFilteredRegularCWComplex(FF);
F:=FF;
ln:=EvaluateProperty(F,"filtration_length");
W:=[];
for k in [1..ln] do
Add(W,FiltrationTerm(F,k));
od;

L:=[];
for k in [1..ln-1] do
map:=function(n,i); return i; end;
cwmap:=Objectify( HapRegularCWMap,
            rec(
            source:=W[k],
            target:=W[k+1],
            mapping:=map));
Add(L,ChainMapOfRegularCWMap(cwmap));
od;


E:=[];
for k in [1..ln] do
Add(E,ChainComplexEquivalenceOfRegularCWComplex(W[k]));
od;


homs:=[];
for k in [1..ln-1] do
f:=Compose( E[k+1][1], Compose(L[k],E[k][2]) );
Add(homs,HomologyVectorSpace(TensorWithIntegersModP(f,p),n));
od;
return LinearHomomorphismsPersistenceMat(homs);
end);
###################################################
###################################################



