#(C) Graham ellis, 2005-2006

#####################################################################
InstallGlobalFunction(EvaluateProperty,
function(X,name)
local
	x;

if not "properties" in NamesOfComponents(X) then return fail; fi;

x:=First(X!.properties,v->v[1]=name);
if x=fail then return fail;
else return x[2]; fi;

end);
#####################################################################

#####################################################################
InstallGlobalFunction(EvenSubgroup,
function(G)
local
	x,y,gens;

gens:=[];
for x in GeneratorsOfGroup(G) do
for y in GeneratorsOfGroup(G) do
Append(gens,[x*y]);
od;
od;

return Group(gens);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ReduceGenerators,
function(gens,G)    		#I should probably use the
                                #Frattini subgroup!
local x,newgens;

for x in gens do
newgens:=DifferenceLists(gens, [x]);
#if Order(Group(Concatenation(newgens,[Identity(G)])))=Order(G) then
if Group(Concatenation(newgens,[Identity(G)]))=G then
return ReduceGenerators(newgens,G); fi;
od;
return gens;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(AbelianInvariantsToTorsionCoefficients,
function(L)
local
	primes, invariants,max,coeffs,i,x;

primes:=SSortedList(Factors(Product(L)));
invariants:=[];
for i in [1..Length(primes)] do
invariants[i]:=Filtered(L,x->(primes[i] in Factors(x)));
od;

max:=Maximum(List(invariants,x->Length(x)));

for x in invariants do
while Length(x)<max do
Append(x,[1]);
od;
od;

coeffs:=[];
for i in [1..max] do
coeffs[i]:= Product(List([1..Length(primes)],j->invariants[j][i]));
od;

return coeffs;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(TorsionGeneratorsAbelianGroup,
function(G)
local
        L,L1,gens,x,y;

L:=IndependentGeneratorsOfAbelianGroup(G);
L1:=SSortedList(StructuralCopy(L));
gens:=[];

while Length(L1)>0 do
x:=L1[1];
RemoveSet(L1,x);
	for y in L1 do
	if Gcd(Order(x),Order(y))=1 then x:=x*y; 
	RemoveSet(L1,y);
	fi;
	od;
Append(gens,[x]);
od;

return gens;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(BigStepLCS,
function(G,n)
local LCS,BSLCS,i;

LCS:=LowerCentralSeries(G);;
BSLCS:=[LCS[1]];

for i in [2..Length(LCS)] do
if Order(LCS[i])=1 or
Order(BSLCS[Length(BSLCS)])/Order(LCS[i])>n then
Append(BSLCS,[LCS[i]]);
fi;
od;

return BSLCS;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(Coclass,
function(G)
local
	facs;

facs:=SSortedList(Factors(Order(G)));
if Length(facs)>1 then 
Print("G should be a prime-power group. \n");
return fail; 
fi;

return LogInt(Order(G),facs[1])-NilpotencyClassOfGroup(G);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(BoundaryMatrix,
function(C,n)
local
        M,i,j;

M:=[];
for i in [1..C!.dimension(n)] do
M[i]:=C!.boundary(n,i);
od;

return TransposedMat(M);

end);
####################################################################

