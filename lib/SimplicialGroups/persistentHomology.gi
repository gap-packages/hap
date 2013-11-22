#############################################################################
#0
#O	PersistentHomologyOfCrossedModule
##	Input:	A crossed module X and an integer number n
##	Output:	The matrix of persistent Betti numbers of X at degree n
##
InstallGlobalFunction(PersistentHomologyOfCrossedModule, function(X,n)
local 
	p,Maps,
	PrimeOne,PrimeTwo,PrimeOneTwo;
	   
	PrimeOne:=PrimeDivisors(Size(HomotopyGroup(X,1)));
	PrimeTwo:=PrimeDivisors(Size(HomotopyGroup(X,2)));
	PrimeOneTwo:=Set(Concatenation(PrimeOne,PrimeTwo));
	if Length(PrimeOneTwo) <>1 then
		return fail;
	fi;
	
	p:=PrimeOneTwo[1];
	Maps:=HomotopyLowerCentralSeriesOfCrossedModule(X);
	Maps:=CatOneGroupByCrossedModule(Maps);
	Maps:=NerveOfCatOneGroup(Maps,n+1);
	Maps:=ChainComplexOfSimplicialGroup(Maps);
	Maps:=List(Maps,f->TensorWithIntegersModP(f,p));
	Maps:=List(Maps,f->HomologyVectorSpace(f,n));
	return LinearHomomorphismsPersistenceMat(Maps);
end);
##
################### end of PersistentHomologyOfCrossedModule ################

#############################################################################
#0
#O	PersistentHomology
##	Input:	A crossed module X and an integer number n
##	Output:	The matrix of persistent Betti numbers of X at degree n
##
InstallOtherMethod(PersistentHomology, "Persistent homology of crossed modules",
[IsHapCrossedModule,IsInt], function(X,n)

	return PersistentHomologyOfCrossedModule(X,n);
end);
##
################### end of PersistentHomology ###############################

#############################################################################
#0
#O	PersistentHomology
##	Input:	A cat-1-group C and an integer number n
##	Output:	The matrix of persistent Betti numbers of C at degree n
##
InstallOtherMethod(PersistentHomology, "Persistent homology of cat-1-groups",
[IsHapCatOneGroup,IsInt], function(C,n)
local  XC;
	
	XC:= CrossedModuleByCatOneGroup(C);
	return PersistentHomologyOfCrossedModule(XC,n);
end);
##
################### end of PersistentHomology ###############################