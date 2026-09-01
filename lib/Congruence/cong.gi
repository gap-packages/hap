MakeReadWriteGlobal("HAP_transversalmethodtoggle");
HAP_transversalmethodtoggle:=true;;
MakeReadOnlyGlobal("HAP_transversalmethodtoggle");

############################################################
############################################################
InstallMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsMatrixGroup,IsMatrixGroup],
function(H,HH);
if not ( IsHapSL2ZSubgroup(H) or IsHapSL2ZSubgroup(HH)) then 
   TryNextMethod(); 
fi;
if H=SL(2,Integers) then
   if HAP_transversalmethodtoggle then
       return HAP_TransversalCongruenceSubgroupInAmbientGroup(H,HH);
   else
       return HAP_RightTransversalSLnZSubgroups_alt(H,HH); #Works well!
   fi;
else
   return HAP_RightTransversalSLnZSubgroups(H,HH);
fi;
end);
############################################################
############################################################

####################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroupInAmbientGroup,
function(G,H)
local poscan, bool;

bool:=IsBound(H!.transversal);
if bool then bool:= not (H!.transversal = fail); fi;

#if not IsBound(H!.transversal) then
if not bool then 
return HAP_TransversalCongruenceSubgroupInAmbientGroupSlow(G,H);
fi;

poscan:=function(g); return H!.cosetPos(g^-1); end;

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSLnZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=H!.transversal,
               poscan:=poscan ));

end);
####################################################################
###################################################################


###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroupInAmbientGroupSlow,
function(G,H)
local tree,InH,v,p,g,s,n,q,vv,gens,
      nodes, nodesinv, leaves, ambientGenerators, InLowDegreeNodesModH,
      one, poscan, nind;

ambientGenerators:=H!.ambientGenerators;
one:=ambientGenerators[1]^0;
tree:=[1 ];
cnt:=1;
leaves:=NewDictionary(one,true,G);
nodes:=[one];
AddDictionary(leaves,one,1);

InH:=H!.membershipLight;

###########################################
InLowDegreeNodesModH:=function(g)
local x,gg,B1,B2;

gg:=g^-1;

for x in nodes do
if InH(x*gg) then return x; fi;
od;

return false;
end;
###########################################



############Tree Construction########################
while Size(leaves)>0 do
vv:=leaves!.entries[1];
v:=vv[1];
    for s in [1..Length(ambientGenerators)] do
        g:=v*ambientGenerators[s];
        q:=InLowDegreeNodesModH(g);
        if q=false then
         Add(nodes,g);
         AddDictionary(leaves,g,Length(nodes));
            p:=LookupDictionary(leaves,v);
            Add(tree,[p, s]);
        fi;
    od;
RemoveDictionary(leaves,v);
od;
#####################################################

nodes:=Filtered(nodes,g-> g in G);
nodesinv:=List(nodes,g->g^-1);
nind:=[1..Length(nodes)];

####################################################
poscan:=function(x)
local i;

for i in nind do
if InH(  x*nodesinv[i]  ) then return i; fi;
od;
return fail;
end;
####################################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSLnZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=nodes,
               poscan:=poscan ));
end);
###################################################################
###################################################################


############################################################
############################################################
InstallGlobalFunction(HAP_RightTransversalSLnZSubgroups,
function(H,HH)
local F, RHH, poscan, G; 

G:=SL(2,Integers);
#RHH:=HAP_RightTransversalSLnZSubgroups(HAP_CongruenceSubgroupGamma0(1),HH);
RHH:=RightTransversal(HAP_CongruenceSubgroupGamma0(1),HH);
F:=Filtered(RHH,x->x in H);

poscan:=function(g)
local i, gg;
gg:=g^-1;
for i in [1..Length(F)] do
if F[i]*gg in HH then return i; fi;
od;
end;

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSLnZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := H,
               subgroup := HH,
               cosets:=F,
               poscan:=poscan ));
end);
############################################################
############################################################

###########################################################
###########################################################
InstallOtherMethod(PositionCanonical,
"for HapRightTransversals of subrougs in SL(2,Z) or SL(2,O)",
[IsHapRightTransversalSLnZSubgroup,IsObject],
function(R,x)
return R!.poscan(x);
end);
###########################################################
###########################################################



################################################
################################################
InstallOtherMethod(IndexNC,
"index for HapSLOSubgroups",
[IsMatrixGroup,IsHapSL2ZSubgroup],
function(G,H);
return Length(RightTransversal(G,H));
end);
################################################
################################################

################################################
################################################
InstallOtherMethod(IndexNC,
"index for HapSLOSubgroups",
[IsMatrixGroup,IsHapSL2OSubgroup],
function(G,H);
return Length(RightTransversal(G,H));
end);
################################################
################################################


