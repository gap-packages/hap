#(C) Graham Ellis, 2012


#####################################################################
InstallOtherMethod(Homology,"Integral homology of cat-1-groups",
[IsHapCatOneGroup, IsInt],
function(C,n)
local N,K;

N:=NerveOfCatOneGroup(QuasiIsomorph(C),n+1);
K:=ChainComplexOfSimplicialGroup(N);
return Homology(K,n);

end );
#####################################################################
