IsHapCongruenceSubgroup:=NewFilter("IsHapCongruenceSubgroup");;
IsHapRightTransversalSL3ZSubgroup:=NewFilter("IsHapRightTransversalSL3ZSubgroup");;
DeclareGlobalFunction("HAP_GenericCongruenceSubgroup");
DeclareOperation("FiniteProjectiveLine",[IsInt]);
DeclareOperation("FiniteProjectivePlane",[IsInt]);
DeclareOperation("FiniteProjectiveLine_alt",[IsInt]);
DeclareOperation("FiniteProjectivePlane_alt",[IsInt]);
DeclareOperation("FiniteProjectiveSpace3",[IsInt]);
DeclareOperation("FiniteProjectiveSpace3_alt",[IsInt]);
DeclareOperation("FiniteProjectiveSpace",[IsInt,IsInt]);
DeclareOperation("FiniteProjectiveSpace_alt",[IsInt,IsInt]);
DeclareGlobalFunction("HAP_FiniteProjectiveLineIntegers");
DeclareGlobalFunction("HAP_FiniteProjectiveLineIntegers_alt");
DeclareGlobalFunction("HAP_FiniteProjectivePlaneIntegers");
DeclareGlobalFunction("HAP_FiniteProjectivePlaneIntegers_alt");
DeclareGlobalFunction("HAP_FiniteProjectiveSpace3Integers");
DeclareGlobalFunction("HAP_FiniteProjectiveSpace3Integers_alt");
DeclareGlobalFunction("HAP_FiniteProjectiveSpaceIntegers");
DeclareGlobalFunction("HAP_FiniteProjectiveSpaceIntegers_alt");
DeclareGlobalFunction("HAP_MatrixInSLn_Hermite");
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

