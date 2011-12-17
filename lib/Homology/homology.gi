#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(Homology,
function(X,N)

if EvaluateProperty(X,"characteristic")=0 then
return IntegralHomology(X,N); fi;

if 	IsPrimeInt(EvaluateProperty(X,"characteristic"))
#	and
#	EvaluateProperty(X,"type")="chainComplex"
then
	return ModularHomology(X,N); 
else
	Print("Induced morphisms in mod p homology are not yet implemented\n");
	return fail;
fi;

end );
#####################################################################
