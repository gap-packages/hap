IsHapCongruenceSubgroup:=NewFilter("IsHapCongruenceSubgroup");;
IsHapRightTransversalSL3ZSubgroup:=NewFilter("IsHapRightTransversalSL3ZSubgroup");;
DeclareGlobalFunction("HAP_GenericCongruenceSubgroup");
DeclareOperation("FiniteProjectiveLine",[IsInt]);
DeclareOperation("FiniteProjectivePlane",[IsInt]);
DeclareOperation("FiniteProjectiveLine_alt",[IsInt]);
DeclareGlobalFunction("HAP_FiniteProjectiveLineIntegers");
DeclareGlobalFunction("HAP_FiniteProjectiveLineIntegers_alt");
DeclareGlobalFunction("HAP_FiniteProjectivePlaneIntegers");
DeclareGlobalFunction("HAP_SL3ZSubgroupTree_fast");
DeclareOperation("HAPCongruenceSubgroupGamma0",[IsInt,IsInt]);
DeclareOperation("HAPCongruenceSubgroupGamma0_alt",[IsInt,IsInt]);
DeclareOperation("HAPIntersectionConjugatedCongruenceSubgroup",[IsHapCongruenceSubgroup,IsMatrix]);
DeclareOperation("HAPCongruenceSubgroupTree",[IsHapCongruenceSubgroup]);
DeclareGlobalFunction("Gamma0_SL3ZTopRationalHomology");

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

