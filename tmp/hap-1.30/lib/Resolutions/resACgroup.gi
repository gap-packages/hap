#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAlmostCrystalGroup,
function(G,K)
local
	GhomP,P,T,RP,RT,RG,PC,PhomPC;

if not (IsAlmostCrystallographic(G) and IsPcpGroup(G)) then
Print("This function can only be applied to Almost Crystallographic pcp groups. \n"); return fail;
fi;

GhomP:=NaturalHomomorphismOnHolonomyGroup(G);
P:=Image(GhomP);
T:=Kernel(GhomP);

PhomPC:=IsomorphismPcGroup(P);
PC:=Image(PhomPC); #For some reason, good resolutions of finite pc groups
		   #seem to be easier to construct than for finite pcp groups.

RP:=ResolutionFiniteGroup(PC,K);

RP!.group:=P;
RP!.elts:=List(RP!.elts,x->PreImageElm(PhomPC,x));
RT:=ResolutionNilpotentGroup(T,K);

return ResolutionExtension(GhomP,RT,RP,"Don't Test Finiteness");



end);
#####################################################################
