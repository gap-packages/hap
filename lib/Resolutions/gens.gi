#####################################################
#####################################################
InstallMethod(Generators,
"Group generators derived froma free resolution",
[IsHapResolution],
function(R)
local P;
P:=PresentationOfResolution(R);
return R!.elts{P.gens};
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(Generators,
"Group generators, as 2x2 matrices, derived from a free resolution of a subgroup of GL(2,R) with R a ring of quatratic integers ",
[IsHapResolution,IsBool],
function(R,bool)
local P, G, d, gens;
G:=R!.group;
P:=PresentationOfResolution(R);
gens:= R!.elts{P.gens};
if not (bool and IsBound(G!.bianchiInteger) and IsMatrixGroup(G)) then 
return gens; fi;
d:=G!.bianchiInteger;
return List(gens,x->HAP_4x4MatTo2x2Mat(x,d) );
end);
#####################################################
#####################################################

