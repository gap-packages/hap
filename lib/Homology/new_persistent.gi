#RT:=0;

###################################################
###################################################
InstallGlobalFunction(PersistentBettiNumbersViaContractions,
function(arg)
local F,N,p,n, W, L, E, chnmap, cwmap, map, ln, k, f, homs,bool;

N:=arg[2];
p:=arg[3];
if Length(arg)=4 then bool:=true; else bool:=false; fi;
if not (IsInt(N) or IsList(N)) then return fail; fi;
if IsInt(N) then n:=N; fi;

F:=arg[1];
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

if IsInt(N) then
homs:=[];
for k in [1..ln-1] do
f:=Compose( E[k+1][1], Compose(L[k],E[k][2]) );
Add(homs,HomologyVectorSpace(TensorWithIntegersModP(f,p),n));
od;

if bool then return homs; fi;
return LinearHomomorphismsPersistenceMat(homs);
fi;

homs:=[];
for n in N do
homs[n+1]:=[];
od;

N:=SortedList(N); 

for k in [1..ln-1] do
f:=Compose( E[k+1][1], Compose(L[k],E[k][2]) );
for n in N do
#RT:=0-Runtime();
Add(homs[n+1],HomologyVectorSpace(TensorWithIntegersModP(f,p),n));
#RT:=RT+Runtime(); Print([k,n,RT],"\n");
od;
od;

if bool then return homs; fi;
homs:=List(N,n-> LinearHomomorphismsPersistenceMat(homs[n+1]) );
return homs;
end);
###################################################
###################################################



