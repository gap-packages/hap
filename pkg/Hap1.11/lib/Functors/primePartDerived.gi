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
	HPpres,G1,epi,HPP,rho;

C:=F(R);
P:=StructuralCopy(R!.group);
P1:=Normalizer(G,P);
prime:=Factors(Order(P))[1];

core:=[];
for x in P do
if Order(x)=prime then AddSet(core,x); fi;
od;

if not IsPrimeInt(Order(P)) then
DCRS1:=List(DoubleCosetRepsAndSizes(G,P1,P1),x->x[1]);
else
DCRS1:=[];
fi;

if Order(P1)>Order(P) then
Append(DCRS1,Filtered(ReduceGenerators(GeneratorsOfGroup(P1),P1),
x->not x in P));
fi;

DCRS:=[];
for x in DCRS1 do  #I've forgotten what all this means!!
for y in core do
if x*y*x^-1 in core then Add(DCRS,x); break; fi;
od;od;
DCRSpruned:=[];

HP:=GroupHomomorphismByFunction(P,P,x->x);
HP:=EquivariantChainMap(R,R,HP);
HP:=F(HP);


HP:=Homology(HP,n);

HP:=Source(HP);

HPrels:=[Identity(HP)];

if Length(AbelianInvariants(HP))=0 then return []; fi;

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
Append(DCRSpruned,[ReduceGenerators(L,GroupL)]); 
od;
DCRSpruned:=Filtered(DCRSpruned,x->Length(x)>0);

for L in DCRSpruned do
K:=Intersection(P,P^L[1]);
gensK:=ReduceGenerators(GeneratorsOfGroup(K),K);

G1:=Group(Concatenation(gensK,[Identity(P)]));
if Order(G1)<64 and n<4 then 	##NEED TO FIND AN "OPTIMAL" CHOICE HERE
S:=ResolutionFiniteGroup(gensK,n+1);
else
S:=ResolutionNormalSeries(LowerCentralSeries(G1),n+1);
fi;
#S:=TietzeReducedResolution(S);


if not (Homology(F(S),n)=[]) then

f:=GroupHomomorphismByFunction(K,P,x->x);

HKhomHPK:=Homology(F(EquivariantChainMap(S,R,f)),n);

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
HKxhomHPKx:=Homology(F(EquivariantChainMap(S,R,fx)),n);
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
od;

if IsPcpGroup(HP)  then 
HPP:=HP/Group(SSortedList(HPrels));
else
epi:=EpimorphismNilpotentQuotient(HP,1);
HPP:=Range(epi)/Group(SSortedList(List(HPrels,x->Image(epi,x))));
#HPP:=HP/Group(SSortedList(HPrels));
fi;

return AbelianInvariants(HPP);
end);
#####################################################################

