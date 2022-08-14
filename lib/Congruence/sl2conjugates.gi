
###################################################################
###################################################################
InstallGlobalFunction(HAP_GenericSL2ZConjugatedSubgroup,
function()
local type, G;

    type := NewType( FamilyObj([[[1,0],[0,1]]]),
                     IsGroup and
                     IsAttributeStoringRep and
                     IsFinitelyGeneratedGroup and
                     IsMatrixGroup and
                     IsHapSL2ZConjugatedSubgroup);

G:=rec(
    sl2Zsubgroup:=fail,
    conjugator:=fail,
    conjugatorInverse:=fail,
    membership:= fail,
    generators:= fail,
    cosetRep:= fail,
    cosetPos:= fail,
    name:="Conjugated congruence subgroup");

ObjectifyWithAttributes(G, type,
DimensionOfMatrixGroup, 2,
IsMatrixGroup, true,
IsFinite, false,
IsHapSL2ZConjugatedSubgroup, true,
IsHapSL2ConjugatedSubgroup, true);

return G;
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod(\in,
"membership test for Hap_SL2Subgroups",
[IsMatrix, IsHapSL2ConjugatedSubgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(x,G)
local g,H;
H:=G!.sl2Zsubgroup;
g:=G!.conjugatorInverse;
return  H!.membership(x^g);
end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallMethod(GeneratorsOfGroup,
"Generating set for Hap_SL2ZSubgroups",
[IsHapSL2ConjugatedSubgroup and IsGroup],
1000000,  #There must be a better way to ensure this method is used!
function(G)
local H,g;
H:=G!.sl2Zsubgroup;
g:=G!.conjugator;
return  List(GeneratorsOfGroup(H),x->x^g);
end);
###################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_ConjugatedCongruenceSubgroupGamma0,
function(n,g)
local H, G, gg, membership;

G:=HAP_GenericSL2ZConjugatedSubgroup();
gg:=g^-1;
H:=HAP_CongruenceSubgroupGamma0(n);

###################################################
membership:=function(x);
return x^gg in H;
end;
###################################################

G!.membership:=membership;
G!.sl2Zsubgroup:=H;
G!.conjugator:=g;
G!.conjugatorInverse:=gg;

G!.name:="ConjugatedCongruenceSubgroupGamma0";
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(HAP_ConjugatedCongruenceSubgroup,
function(H,g)
local  G, gg, membership;

G:=HAP_GenericSL2ZConjugatedSubgroup();
gg:=g^-1;
#H:=HAP_CongruenceSubgroupGamma0(n);

###################################################
membership:=function(x);
return x^gg in H;
end;
###################################################

G!.membership:=membership;
G!.sl2Zsubgroup:=H;
G!.conjugator:=g;
G!.conjugatorInverse:=gg;

G!.name:="ConjugatedCongruenceSubgroup";
return G;

end);
###################################################################
###################################################################

###################################################################
###################################################################
InstallGlobalFunction(ResolutionSL2ZConjugated,
function(n,g)
local R, G, gg, elts, gelts, eltfun, posfun;

R:=ResolutionSL2Z_alt(n);
G:=HAP_ConjugatedCongruenceSubgroupGamma0(1,g);
gg:=g^-1;
elts:=R!.elts;
gelts:=ElementsLazy(G);


R!.elts:=gelts;
R!.group:=G;

##############################
eltfun:=function(i)
local x;
x:=elts[i];
if x=fail then return x;
else
return elts[i]^g;
fi;
end;
##############################

##############################
posfun:=function(x);
return elts!.posfun( x^gg);
end;
##############################


gelts!.eltfun:=eltfun;
gelts!.posfun:=posfun;
return R;
end);
###################################################################
###################################################################


############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsHapSL2ZConjugatedSubgroup,IsMatrixGroup],
100000,
function(G,GG)
local R, R2, H, K,gg,g,poscan;

H:=G!.sl2Zsubgroup;
g:=G!.conjugator;
gg:=G!.conjugatorInverse;
R:=RightTransversal(H,GG^gg);

R2:=List(R,x->x^g);

##########################################
poscan:=function(x);
return PositionCanonical(R,x^gg);
end;
##########################################
return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := GG,
               cosets:=R2,
               poscan:=poscan ));
end);
############################################################
############################################################



