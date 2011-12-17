#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAlmostCrystalGroup,
function(G,K)
local
	GhomP,P,T,RP,RT,RG;

if not (IsAlmostCrystallographic(G) and IsPcpGroup(G)) then
Print("This function can only be applied to Almost Crystallographic pcp groups. \n"); return fail;
fi;

GhomP:=NaturalHomomorphismOnHolonomyGroup(G);
P:=Image(GhomP);
T:=Kernel(GhomP);

RP:=ResolutionFiniteGroup(P,K);
RT:=ResolutionNilpotentGroup(T,K);

return ResolutionExtension(GhomP,RT,RP,"Don't Test Finiteness");



end);
#####################################################################
