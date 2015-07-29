#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionAlmostCrystalQuotient,
function(arg)
local
	G,K,KK,bool,
	GhomP,P,T,Derived,i, RGD,
	pcpGD, GD, GDhomPCGD,PCGD,GhomGD, GDhomG, GDhomP, TD, gensP, RP,RTD,
        iso, isoTD;

G:=arg[1];
K:=arg[2];
KK:=arg[3];
if Length(arg)>3 then bool:=false; else bool:=true; fi; 


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
GD:=Image(GhomGD);

if IsNilpotent(GD) and bool then 

	if IsFinite(GD) then
	GDhomPCGD:=IsomorphismPermGroup(GD);
	PCGD:=Image(GDhomPCGD);
	RGD:=ResolutionNormalSeries(BigStepLCS(PCGD,4),K);
					###MIGHT WANT TO VARY THIS 4
	RGD!.quotientHomomorphism:=
	GroupHomomorphismByFunction(G,PCGD,
		x->Image(GDhomPCGD,Image(GhomGD,x)));
	return RGD;

	else
	RGD:=ResolutionNilpotentGroup(GD,K);
	fi;


else

GDhomG:=GroupHomomorphismByImagesNC(GD,G,    
GeneratorsOfGroup(GD),GeneratorsOfPcp(pcpGD));

gensP:=List(GeneratorsOfPcp(pcpGD),x->Image(GhomP,x));
GDhomP:=GroupHomomorphismByImagesNC(GD,P,          
GeneratorsOfGroup(GD),gensP);

TD:=Kernel(GDhomP);

RP:=ResolutionFiniteGroup(P,K);
if 0 in AbelianInvariants(TD) then
RTD:=ResolutionNilpotentGroup(TD,K);
else
RTD:=ResolutionFiniteGroup(TD,K);
#TD:=Group(TorsionGeneratorsAbelianGroup(TD));
#iso:=IsomorphismPcGroup(TD);
#isoTD:=Image(iso);
#RTD:=ResolutionFiniteGroup(isoTD,K);
#RTD!.elts:=List(RTD!.elts,x->PreImagesRepresentative(iso,x));
fi;


RGD:=ResolutionExtension(GDhomP,RTD,RP,"Don't Test Finiteness");
fi;

RGD!.quotientHomomorphism:=GhomGD;

return RGD;
end);
#####################################################################

#####################################################################
InstallGlobalFunction(RelativeCentralQuotientSpaceGroup,
function(G,c)
local Q;
if not IsSpaceGroup(G) then
Print("This function must be applied to a space group.\n"); return fail; fi;
Q:=Image(IsomorphismPcpGroup(G));;
Q:=ResolutionAlmostCrystalQuotient(Q,0,c)!.group;;
return Q;
end);
#####################################################################

