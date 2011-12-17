#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(Cohomology,
function(X,N)

if EvaluateProperty(X,"characteristic")=0 then
return IntegralCohomology(X,N); 

else
Print("Cohomology in characteristic p>0 is not yet implemented\n");
	return fail;
fi;

end );
#####################################################################
