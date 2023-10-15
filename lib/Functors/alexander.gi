#####################################################################
InstallGlobalFunction(HAP_PrimePartModified,
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
HKhomHPK:=Homology(F(eqmap),n);

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
HKxhomHPKx:=Homology(F(eqmap),n);
HKx:=Source(HKxhomHPKx);
HPKx:=Parent(Range(HKxhomHPKx));
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
fi;
od;

f:=NormalClosure(HP,Group(SSortedList(HPrels)));
return NaturalHomomorphismByNormalSubgroup(HP,f);
end);
#####################################################################

#####################################################################
InstallGlobalFunction(HAP_SylowSubgroups,
function(f,p)
local H,G, SH,SG,P,Q;
H:=Source(f);
G:=Target(f);
SH:=SylowSubgroup(H,p);
SG:=SylowSubgroup(G,p);
P:=Image(f,SH);
for Q in SG^G do
if IsSubgroup(Q,P) then return [Q,SH]; fi;
od;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(PrimePartDerivedFunctorHomomorphism,
function(arg)
local f,RG,RH,F,n,SH,SG,P,Q,ef,A,B,ans,G,H,x,iso1,iso2,hom,cf,hf;
f:=arg[1];
RG:=arg[2];  #resolution for sylow subgroup of G
RH:=arg[3];  #resolution for sylow subgroup of H
F:=arg[4];
n:=arg[5];


H:=Source(f);
G:=Target(f);
SG:=RG!.group;
SH:=RH!.group;
A:=HAP_PrimePartModified(G,RG,F,n);
B:=HAP_PrimePartModified(H,RH,F,n);

if Length(arg)=6 then cf:=arg[6];
else
ef:=EquivariantChainMap(RH,RG,f);
cf:=F(ef);
fi;
hf:=Homology(cf,n);

iso1:=GroupHomomorphismByImages(Source(B),Source(hf),GeneratorsOfGroup(Source(B)),GeneratorsOfGroup(Source(hf)));
iso2:=GroupHomomorphismByImages(Target(hf),Source(A),GeneratorsOfGroup(Target(hf)),GeneratorsOfGroup(Source(A)));


hom:=GroupHomomorphismByFunction(Source(B),Source(A),x->
   Image(iso2,Image(hf,Image(iso1,x)))  );


return
GroupHomomorphismByFunction(Target(B),Target(A),
x->Image(A,Image(hom,PreImagesRepresentative(B,x))));
end);
#####################################################################

#####################################################################
InstallGlobalFunction(ModPCohomologyRing_alt,
function(G,R)
local P,A,B,S,row,k,Bases,mat,hom,f,F,n,b,x,ef,cf,p;

P:=R!.group;
p:=PrimePGroup(P);
A:=ModPCohomologyRing(R);
B:=Basis(A);
Bases:=List([1..Length(R)],i->[]);
for b in B{[2..Length(B)]} do
Add(Bases[A!.degree(b)],b);
od;
F:=function(R); return TensorWithIntegersModP(R,p); end;
f:=GroupHomomorphismByFunction(P,G,x->x);
ef:=EquivariantChainMap(R,R,f);
cf:=F(ef);

S:=[B[1]];

for n in [1..Length(R)-1] do
hom:=PrimePartDerivedFunctorHomomorphism(f,R,R,F,n,cf);;
mat:=HomomorphismAsMatrix(hom);
mat:=TransposedMat(mat);
for row in mat do
   x:=0*B[1];
   for k in [1..Length(row)] do
      x:=x+row[k]*Bases[n][k];
   od;
   Add(S,x);
od;
od;

S:=Subalgebra(A,S);
S!.degree:=A!.degree;
return S;
end);
#####################################################################

