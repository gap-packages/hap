#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(GroupHomology,
function(arg)
local
		G, gensG, N, p, D,
		HomologyPrimePowerGroup,
		HomologyGenericGroup,
		HomologyGenericGroupModP,
		HomologySmallGroup,
		HomologyArtinGroup,
		HomologyAbelianGroup;

############################### INPUT DATA ##########################
if IsList(arg[1]) then D:=arg[1]; G:=false; 
else
	if IsGroup(arg[1]) then G:=arg[1]; 
	gensG:=ReduceGenerators(GeneratorsOfGroup(G),G); 
	D:=false;
	else
	Print("ERROR: first variable must be a group or a Coxeter diagram. \n");
	fi;
fi;

if IsInt(arg[2]) then N:=arg[2]; else
Print("ERROR: second variable must be a positive integer.\n");
fi;

if Length(arg)>2 then p:=arg[3]; else p:=0; fi;
if not (p=0 or IsPrimeInt(p)) then 
Print("ERROR: third variable should be a prime integer. \n");
fi;

############################## DATA INPUT ###########################


#####################################################################
#####################################################################
HomologySmallGroup:=function()
local R;

R:=ResolutionFiniteGroup(gensG,N+1,false,p);
if p=0 then
return Homology(TensorWithIntegers(R),N);
else
return Homology(TensorWithIntegersModP(R,p),N);
fi;

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyPrimePowerGroup:=function()
local
		R, C, L, LL, x;

if Order(G)<32 then R:=ResolutionFiniteGroup(gensG,N+1,false,p); 
else
L:=LowerCentralSeries(G);

R:=ResolutionNormalSeries(L,N+1,false,p);
fi;

if p=0 then
C:=TensorWithIntegers(R);
return Homology(C,N);
else
C:=TensorWithIntegersModP(R,p);
return Homology(C,N);
fi;

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyAbelianGroup:=function()
local L,R,C;

L:=AbelianInvariantsToTorsionCoefficients(AbelianInvariants(G));
R:=ResolutionAbelianGroup(L,N+1);

if p=0 then
C:=TensorWithIntegers(R);
else
C:=TensorWithIntegersModP(R,p);
fi;

return Homology(C,N);
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyGenericGroup:=function()
local
	primes, q, S, gens, functor, R, H, TorsionCoeffs;

TorsionCoeffs:=[];
primes:= SSortedList(Factors(Order(G)));


for q in primes do

S:=SylowSubgroup(G,q);
gens:=GeneratorsOfGroup(S);
gens:=ReduceGenerators(gens,S);
R:=ResolutionFiniteGroup(gens,N+1,false,p);

H:=Homology(TensorWithIntegers(R),N);

if 
Order(Centralizer(G,S))=Order(Center(S))*Order(G)/Order(S)
or 
Length(H)=0 
then Append(TorsionCoeffs,H); 
else
Append(TorsionCoeffs,
PrimePartDerivedFunctor(G,R,TensorWithIntegers,N));
fi;
od;

return TorsionCoeffs;
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HomologyArtinGroup:=function()
local R;

R:=ResolutionArtinGroup(D,N+1);;
R:=TensorWithIntegers(R);;
return Homology(R,N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyGenericGroupModP:=function()
local
        S, gens, R, H, functor;

#####################################################################
functor:=function(X);
return TensorWithIntegersModP(X,p);
end;
#####################################################################

	S:=SylowSubgroup(G,p);
	gens:=GeneratorsOfGroup(S);
	gens:=ReduceGenerators(gens,S);
	R:=ResolutionFiniteGroup(gens,N+1,false,p);

	H:=Homology(functor(R),N);

	if
	Order(Centralizer(G,S))=Order(Center(S))*Order(G)/Order(S)
	or
	H=0
	then return H;
	else
	return	
	PrimePartDerivedFunctor(G,R,functor,N);
	fi;
	
end;
#####################################################################
#####################################################################



if IsList(D) then
return HomologyArtinGroup(); fi;

if IsAbelian(G) then
return HomologyAbelianGroup(); fi;

if IsPrimePowerInt(Order(G)) then
return HomologyPrimePowerGroup(); fi;

if Order(G)<16 then
return HomologySmallGroup(); fi;

if p=0 then
return HomologyGenericGroup(); fi;

return HomologyGenericGroupModP();

end);
