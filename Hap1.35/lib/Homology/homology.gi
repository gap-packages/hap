#(C) Graham Ellis, 2005-2006


#####################################################################
#InstallGlobalFunction(Homology,

InstallMethod(Homology,"homology calculations",
[IsHapChain, IsInt],0,
function(X,N)

if EvaluateProperty(X,"characteristic")=0 then
return IntegralHomology(X,N); fi;

if 	IsPrimeInt(EvaluateProperty(X,"characteristic"))
	or (EvaluateProperty(X,"characteristic")=-1/2
	    and 
	    EvaluateProperty(X,"type")="chainComplex")
then
	return ModularHomology(X,N); 
else
	Print("Induced morphisms in rational homology are not yet implemented\n");
	return fail;
fi;

end );
#####################################################################

#####################################################################
InstallMethod(Homology,"homology calculations",
[IsHapSparseChainComplex and IsHapChain, IsInt],0,
function(C,n)
return Homology( SparseChainComplexToChainComplex(C) ,n);
end);
#####################################################################
