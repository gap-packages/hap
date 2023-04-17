#############################
InstallGlobalFunction(HAP_MyIsInfiniteFpGroup,
function(G)
local L;
if GeneratorsOfGroup(G)=[] then return false; fi;
#L:=LowIndexSubgroupsFpGroup(G,3);    #switch 6 to 3 for speed
L:=LowIndexSubgroupsFpGroup(G,6);
Apply(L,AbelianInvariants);
if 0 in Flat(L) then return true;
else return false; fi;    #Careful: false should really be fail
#else return fail; fi;
end);
#############################

#############################
InstallGlobalFunction(HAP_MyIsFiniteFpGroup,
function(G);
if HAP_MyIsInfiniteFpGroup(G)=true then return false; fi;
return IsFinite(G);
end);
#############################

#############################
InstallGlobalFunction(HAP_MyIsBieberbachFpGroup,
function(F)
local L;
L:=LowIndexSubgroupsFpGroup(F,6); #6 suffices for 3D groups
L:=Filtered(L,g->IsNormal(F,g));
L:=List(L,g->Image(IsomorphismFpGroup(g)));
L:=Filtered([1..Length(L)],i->AbelianInvariants(L[i])=[0,0,0] and Length(RelatorsOfFpGroup(L[i]))>=3);
if Length(L)>0 then return true; fi;;
return false;
end);
#############################

#############################
InstallGlobalFunction(ManifoldType,
function(M)
local F;
if Dimension(M)<>3 then
Print("Argument must be a 3-manifold.\n"); return fail; fi;
if not IsClosedManifold(M) then
Print("Argument must be a 3-manifold.\n"); return fail; fi;

F:=FundamentalGroup(M);
if HAP_MyIsInfiniteFpGroup(F)=true then
  if HAP_MyIsBieberbachFpGroup(F) then return "euclidean";
  else return "other"; fi;
fi;
if HAP_MyIsFiniteFpGroup(F)=true then return "spherical";
else return fail; fi;
end); 
#############################
