#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(PrimePartDerivedFunctor,
function(G,R,F,n)
local
	C,P, DCRS, DCRS1, DCRSpruned,L,Y,GroupL,
	X, K, gensK, S, f,fx,
	HP, HK, HPK, HKhomHPK, HPKhomHP, HKhomHP,
	HKx,HPKx, 
	HKxhomHPKx, HPKxhomHP, HKxhomHP, HKhomHKx,  HKhomHP2,
	HPrels, x, y, i,prime, core, conjs, conjelt;

C:=F(R);

P:=R.group;
prime:=Factors(Order(P))[1];
core:=[];
for x in P do
if Order(x)=prime then AddSet(core,x); fi;
od;

DCRS1:=DoubleCosetRepsAndSizes(G,P,P);
DCRS:=[];
for x in DCRS1 do
for y in core do
if x[1]*y*x[1]^-1 in core then Append(DCRS,[x]); break; fi;
od;od;
DCRSpruned:=[];

HP:=GroupHomomorphismByFunction(P,P,x->x);
HP:=EquivariantChainMap(R,R,HP);
HP:=F(HP);
HP:=Homology(HP,n);
HP:=Source(HP);
HPrels:=[Identity(HP)];

if Order(HP)=1 then return []; fi;

conjs:=[];
conjelt:=[];
for x in DCRS do
Y:=(Intersection(P,P^x[1]));
AddSet(conjs,Y);
Append(conjelt,[[x[1],Y]]);
od;

for Y in conjs do
L:=Filtered(conjelt,x->x[2]=Y);
L:=List(L,x->x[1]);
GroupL:=Group(L);
Append(DCRSpruned,[ReduceGenerators(L,GroupL)]); 
od;
DCRSpruned:=Filtered(DCRSpruned,x->Length(x)>0);

for L in DCRSpruned do
K:=Intersection(P,P^L[1]);
gensK:=ReduceGenerators(GeneratorsOfGroup(K),K);

S:=ResolutionFiniteGroup(gensK,n+1);

#S:=ResolutionFiniteSubgroup(R,K);
if not (Homology(F(S),n)=[]) then

f:=GroupHomomorphismByFunction(K,P,x->x);
HKhomHPK:=Homology(F(EquivariantChainMap(S,R,f)),n);
HK:=Source(HKhomHPK);

HPK:=Parent(Image(HKhomHPK));
HPKhomHP:=GroupHomomorphismByImagesNC(HPK,HP,GeneratorsOfGroup(HPK),
                                                  GeneratorsOfGroup(HP));
HKhomHP:=GroupHomomorphismByFunction(HK,HP,x->
Image(HPKhomHP, Image(HKhomHPK,x) ) );

for X in L do
fx:=GroupHomomorphismByFunction(K,P,g->Image(f,g)^(X^-1));
HKxhomHPKx:=Homology(F(EquivariantChainMap(S,R,fx)),n);
HKx:=Source(HKxhomHPKx);
HPKx:=Parent(Image(HKxhomHPKx));
HPKxhomHP:=GroupHomomorphismByImagesNC(HPKx,HP,GeneratorsOfGroup(HPKx),
                                                  GeneratorsOfGroup(HP));
HKxhomHP:=GroupHomomorphismByFunction(HKx,HP,x->
Image(HPKxhomHP, Image(HKxhomHPKx,x) ) );
HKhomHKx:=GroupHomomorphismByImagesNC(HK,HKx,GeneratorsOfGroup(HK),GeneratorsOfGroup(HKx));
HKhomHP2:=GroupHomomorphismByFunction(HK,HP,a->
Image(HKxhomHP, Image(HKhomHKx,a)));

for x in GeneratorsOfGroup(HK) do
Append(HPrels, [Image(HKhomHP,x)*Image(HKhomHP2,x)^-1]);
od;

od;
fi;
od;


return AbelianInvariants(HP/NormalClosure(HP,Group(HPrels)));
end);
#####################################################################

