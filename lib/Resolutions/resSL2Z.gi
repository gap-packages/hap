#(C) Graham Ellis and Bui Anh Tuan

################################################
################################################
InstallGlobalFunction(ResolutionSL2Z,
function(p,n)

local C, Rsl, Rsl2, RH, G, H, sl, sl2, 
      Hhomsl, Hhomsl2, RHRsl, RHRsl2, D, K; 

if not IsPrime(p) and not p=0 then
Print("At present this function is only defined for SL(2,Z[1/p]) with p prime or SL(2,Z) with p=0.\n\n");
return fail;
fi;

##################################
C:=ContractibleGcomplex("SL(2,Z)");
Rsl:=ResolutionGTree(C,n);
Rsl!.group:=SL(2,Integers);
sl:=Rsl!.group;
##################################

if p=0 then return
ResolutionGTree(C,n);
fi;

##################################
H:=CongruenceSubgroupGamma0(p);
SetName(H,"Gamma");
RH:=ResolutionFiniteSubgroup(Rsl,H);
RH!.group:=H;
##################################


###################################
sl2:=SL2Z(p);
G:=SL2Z(1/p);
Rsl2:=ConjugatedResolution(Rsl,[[1,0],[0,p]]);
Rsl2!.group:=sl2;
##################################


###################################
Hhomsl:=GroupHomomorphismByFunction(H,sl,x->x);
RHRsl:=EquivariantChainMap(RH,Rsl,Hhomsl);
Hhomsl2:=GroupHomomorphismByFunction(H,sl2,x->x);
RHRsl2:=EquivariantChainMap(RH,Rsl2,Hhomsl2);
D:=[Rsl,Rsl2,[RHRsl,RHRsl2]];
K:=TreeOfResolutionsToContractibleGcomplex(D,G);
###################################

return FreeGResolution(K,n);

end);

###############################################################
###############################################################



