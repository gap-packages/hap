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

#####################################################################
InstallGlobalFunction(PrankAlt,
function(G)
local S, p;

p:=SSortedList(Factors(Order(G)));

if not Length(p)=1 then
Print("G must be a finite p-group. \n");
return fail;
fi;

p:=p[1];
S:=LatticeSubgroups(G);
S:=ConjugacyClassesSubgroups(S);
S:=List(S,x->ClassElementLattice(x,1));
S:=Filtered(S,x->IsElementaryAbelian(x));
S:=List(S,x->Order(x));

return Log(Maximum(S),p);

end);
####################################################################

#####################################################################
InstallGlobalFunction(Prank,
function(G)
local
	pPowers,
	xP,
	x,y,
	p,
	AbSubGrps,
	AbSubGrps1,
	AbSubGrps2,
	S, Sx,
	A,
	N,
	N1,
	reps,
	toggle;

###########################################################
p:=SSortedList(Factors(Order(G)));
if not Length(p)=1 then
Print("G must be a finite p-group. \n");
return fail;
fi;
###########################################################

p:=PrimePGroup(G);

###########################################################
if IsAbelian(G) then return
Length(AbelianInvariants(G)); fi;
###########################################################

pPowers:=[];
xP:=[];



AbSubGrps:=[];
AbSubGrps1:=[GeneratorsOfGroup(Center(G))];
AbSubGrps2:=[];
reps:=[Center(G)];
xP:=Concatenation(List(reps,x->Elements(x)));
xP:=SSortedList(xP);

for x in G do
if Order(x)=p and not (x in xP)  then
Add(pPowers,x);
Append(xP,Elements(Group(x)));
fi;
od;


toggle:=true;
while toggle do
toggle:=false;

	for S in AbSubGrps1 do
	for x in pPowers do
	if not x in Group(S) then
	Sx:=Concatenation(S,[x]);
	A:=Group(Sx);
	if IsAbelian(A) and not (A in reps) then
	Add(AbSubGrps2,Sx); AddSet(reps,A); toggle:=true; fi;
	fi;
	od;
	od;
Append(AbSubGrps,AbSubGrps1);
AbSubGrps1:=AbSubGrps2;
AbSubGrps2:=[];

od;


Apply(AbSubGrps, x->Prank(Group(x)));

return Maximum(AbSubGrps);
end);
#####################################################################

####################################################################
#####################################################################
InstallMethod(Compose,
"for group homomorphisms xoy=x(y)",
[IsGroupHomomorphism,IsGroupHomomorphism],
function(x,y)
if not (Source(x)=Range(y)) then return fail; fi;
return GroupHomomorphismByFunction(Source(y),Range(x),
a->Image(x,Image(y,a)));
end);
#####################################################################

#####################################################################
InstallOtherMethod(Compose,
"for FpG-module homomorphisms xoy=x(y)",
[IsHapFPGModuleHomomorphism,IsHapFPGModuleHomomorphism],
function(x,y)
return CompositionOfFpGModuleHomomorphisms(x,y);
end);
#####################################################################

