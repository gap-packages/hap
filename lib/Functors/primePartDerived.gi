#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(PrimePartDerivedFunctor,
function(G,R,F,n)
local
	C,P, DCRS, DCRS1, DCRSpruned,L,Y,GroupL,
	X, K, gensK, S, f,fx, P1, 
	HP, HK, HPK, HKhomHPK, HPKhomHP, HKhomHP,
	HKx,HPKx, 
	HKxhomHPKx, HPKxhomHP, HKxhomHP, HKhomHKx,  HKhomHP2,
	HPrels, x, y, i,prime, core, conjs, conjelt,CentP,
	HPpres,G1,epi,HPP,rho, bool, eqmap;


C:=F(R);
#P:=StructuralCopy(R!.group);
P:=Group(SmallGeneratingSet(R!.group));
HP:=GroupHomomorphismByFunction(P,P,x->x);
HP:=EquivariantChainMap(R,R,HP);
HP!.conjugator:=Identity(P);
HP:=F(HP);
HP:=Homology(HP,n);
HP:=Source(HP);
HPrels:=[Identity(HP)];
if Length(AbelianInvariants(HP))=0 then return []; fi;

P1:=Normalizer(G,P);

prime:=Factors(Order(P))[1];

if not IsPrimeInt(Order(P)) then
DCRS1:=List(DoubleCosetRepsAndSizes(G,P1,P1),x->x[1]);
else
DCRS1:=[];
fi;

if Order(P1)>Order(P) then
Append(DCRS1,Filtered(ReduceGenerators(GeneratorsOfGroup(P1),P1),
x->not x in P));
fi;

core:=[];
for x in P do
if Order(x)=prime then AddSet(core,x); fi;
od;

DCRS:=[];
for x in DCRS1 do  #I've forgotten what all this means!!
for y in core do
if x*y*x^-1 in core then Add(DCRS,x); break; fi;
od;od;
DCRSpruned:=[];

conjs:=[];
conjelt:=[];
for x in DCRS do
Y:=Intersection(P,P^x);
AddSet(conjs,Y);
Append(conjelt,[[x,Y]]);  #An improvement would be to not save all Y (twice!).
od;

for Y in conjs do
L:=Filtered(conjelt,x->x[2]=Y);
L:=List(L,x->x[1]);
GroupL:=Group(L);
Add(DCRSpruned,ReduceGenerators(L,GroupL)); 
od;
DCRSpruned:=Filtered(DCRSpruned,x->Length(x)>0);


for L in DCRSpruned do
K:=Intersection(P,P^L[1]);
gensK:=ReduceGenerators(GeneratorsOfGroup(K),K);
if not Length(gensK)=0 then

G1:=Group(gensK);

if Order(G1)<64 and n<4 then 	##NEED TO FIND AN "OPTIMAL" CHOICE HERE
S:=ResolutionFiniteGroup(gensK,n+1);
S!.group:=Group(SmallGeneratingSet(S!.group));
else
S:=ResolutionNormalSeries(LowerCentralSeries(G1),n+1);
S!.group:=Group(SmallGeneratingSet(S!.group));
fi;


if not (Homology(F(S),n)=[]) then

f:=GroupHomomorphismByFunction(K,P,x->x);

eqmap:=EquivariantChainMap(S,R,f);
eqmap!.conjugator:=Identity(S!.group);
#HKhomHPK:=Homology(F(EquivariantChainMap(S,R,f)),n);
#Print("F  ",F,"\n");
#Print("eqmap  ",eqmap,"\n");
HKhomHPK:=Homology(F(eqmap),n);
#Print("YES\n");
#################################rho##################
if "twist" in NamesOfComponents(F(R)) then
rho:=F(R)!.twist; 
else
rho:=function(x) return 1; end;
fi;
#################################rho done#############

HK:=Source(HKhomHPK);

HPK:=Range(HKhomHPK);
HPKhomHP:=GroupHomomorphismByImagesNC(HPK,HP,GeneratorsOfGroup(HPK),
                                                  GeneratorsOfGroup(HP));
HKhomHP:=GroupHomomorphismByFunction(HK,HP,x->
Image(HPKhomHP, Image(HKhomHPK,x) ) );

for X in L do
fx:=GroupHomomorphismByFunction(K,P,g->Image(f,g)^(X^-1));
eqmap:=EquivariantChainMap(S,R,fx);
eqmap!.conjugator:=X^-1;
#HKxhomHPKx:=Homology(F(EquivariantChainMap(S,R,fx)),n);
HKxhomHPKx:=Homology(F(eqmap),n);
HKx:=Source(HKxhomHPKx);
HPKx:=Parent(Range(HKxhomHPKx));
HPKxhomHP:=GroupHomomorphismByImagesNC(HPKx,HP,GeneratorsOfGroup(HPKx),
                                                  GeneratorsOfGroup(HP));
HKxhomHP:=GroupHomomorphismByFunction(HKx,HP,x->
Image(HPKxhomHP, Image(HKxhomHPKx,x) )^rho(X) );
HKhomHKx:=GroupHomomorphismByImagesNC(HK,HKx,GeneratorsOfGroup(HK),GeneratorsOfGroup(HKx));
HKhomHP2:=GroupHomomorphismByFunction(HK,HP,a->
Image(HKxhomHP, Image(HKhomHKx,a)));

for x in GeneratorsOfGroup(HK) do
Append(HPrels, [Image(HKhomHP,x)*Image(HKhomHP2,x)^-1]);
od;

od;
fi;
fi;
od;

if IsPcpGroup(HP) or IsPcGroup(HP)  then 
HPP:=HP/Group(SSortedList(HPrels));
else
epi:=EpimorphismNilpotentQuotient(HP,1);
HPP:=Range(epi)/Group(SSortedList(List(HPrels,x->Image(epi,x))));
#HPP:=HP/Group(SSortedList(HPrels));
fi;

return AbelianInvariants(HPP);
end);
#####################################################################

