#(C) Graham Ellis, 2005-2006

#####################################################################
#InstallGlobalFunction(Cohomology,

InstallMethod(Cohomology,"cohomology calculations",
[IsHapCochain, IsInt],0,

function(X,N)

if EvaluateProperty(X,"characteristic")=0 then
return IntegralCohomology(X,N); 
fi;

if      IsPrimeInt(EvaluateProperty(X,"characteristic"))
        or (EvaluateProperty(X,"characteristic")=-1/2
        and
        EvaluateProperty(X,"type")="cochainComplex")
	then
        return ModularCohomology(X,N);
	else
        Print("Induced morphisms in modular or rational homology are not yet implemented\n");
        return fail;
	fi;
end );
#####################################################################

