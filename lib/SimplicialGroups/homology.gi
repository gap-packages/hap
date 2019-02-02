#############################################################################
#0
#O	Homology
##	Input:	A crossed module XC and an integer number n
##	Output:	The integral homology H_n(XC,Z)
##
InstallOtherMethod(Homology, "Homology of crossed modules",
[IsHapCrossedModule,IsInt], function(XC,n)
local  C,D,N,K;

	C:=CatOneGroupByCrossedModule(XC);
	D:=QuasiIsomorph(C);
	N:=NerveOfCatOneGroup(D,n+1);
	K:=ChainComplexOfSimplicialGroup(N);
	return Homology(K,n);
end);

#############################################################################
#0
#O	Homology
##	Input:	A cat-1-group C and an integer number n
##	Output:	The integral homology H_n(C,Z)
##
InstallOtherMethod(Homology, "Homology of cat-1-groups",
[IsHapCatOneGroup,IsInt], function(C,n)
local  D,N,K;

	D:=QuasiIsomorph(C);
	N:=NerveOfCatOneGroup(D,n+1);
	K:=ChainComplexOfSimplicialGroup(N);
	return Homology(K,n);
end);

#############################################################################
#0
#O	Homology
##	Input:	A simplicial group G and an integer number n
##	Output:	The integral homology H_n(G,Z)
##
InstallOtherMethod(Homology, "Homology of simplicial groups",
[IsHapSimplicialGroup,IsInt], function(G,n)
local  K;

	K:=ChainComplexOfSimplicialGroup(G);
	return Homology(K,n);
end);

	