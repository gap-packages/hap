IsHapCongruenceSubgroup:=NewFilter("IsHapCongruenceSubgroup");;
DeclareGlobalFunction("HAP_GenericCongruenceSubgroup");
DeclareOperation("HAPCongruenceSubgroupGamma0",[IsInt,IsInt]);
DeclareOperation("HAPCongruenceSubgroupTree",[IsHapCongruenceSubgroup]);

InstallMethod( ViewObj,
"for HapCongruenceSubgroup",
[IsHapCongruenceSubgroup and IsGroup],
10000000,  #Ensures that this method is chosen
function(G)
Print(G!.name," of ",G!.fam);
 end);

InstallMethod( PrintObj,
"for HapCongruenceSubgroup",
[IsHapCongruenceSubgroup and IsGroup],
100000000, #Ensures that this method is chosen
function(G)
Print(G!.name," of ",G!.fam);
 end);

