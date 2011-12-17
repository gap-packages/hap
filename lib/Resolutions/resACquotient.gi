#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAlmostCrystalQuotient,
function(G,K,KK)
local
	GhomP,P,T,Derived,i, RGD,
	pcpGD, GD, GhomGD, GDhomG, GDhomP, TD, gensP, RP,RTD;

if not (IsAlmostCrystallographic(G) and IsPcpGroup(G)) then
Print("This function can only be applied to Almost Crystallographic pcp groups.  \n"); return fail;
fi;


GhomP:=NaturalHomomorphismOnHolonomyGroup(G);
P:=Image(GhomP);
T:=Kernel(GhomP);
Derived:=T;

for i in [2..KK] do
Derived:=CommutatorSubgroup(G,Derived);
od;

pcpGD:=Pcp(G,Derived);
GhomGD:=NaturalHomomorphism(G,Derived);
GD:=Range(GhomGD);
GDhomG:=GroupHomomorphismByImagesNC(GD,G,
GeneratorsOfGroup(GD),GeneratorsOfPcp(pcpGD));

gensP:=List(GeneratorsOfPcp(pcpGD),x->Image(GhomP,x));
GDhomP:=GroupHomomorphismByImagesNC(GD,P,
GeneratorsOfGroup(GD),gensP);

TD:=Kernel(GDhomP);

RP:=ResolutionFiniteGroup(P,K);
RTD:=ResolutionNilpotentGroup(TD,K);

RGD:=ResolutionExtension(GDhomP,RTD,RP,"Don't Test Finiteness");

RGD.quotientHomomorphism:=GhomGD;

return RGD;
end);
#####################################################################
