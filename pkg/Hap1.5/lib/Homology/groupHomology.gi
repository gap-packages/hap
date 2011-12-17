#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(GroupHomology,
function(arg)
local
		G, gensG, N, p, D, 
		Functor,
		HomologyPrimePowerGroup,
		HomologyGenericGroup,
		HomologySmallGroup,
		HomologyArtinGroup,
		HomologyAbelianGroup,
		HomologyNilpotentPcpGroup;


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
return fail;
fi;

############################## DATA INPUT ###########################

if IsPrime(p) then
Functor:=function(R); return TensorWithIntegersModP(R,p); end;
else
Functor:=TensorWithIntegers;
fi;

#####################################################################
#####################################################################
HomologySmallGroup:=function()
local R;

R:=ResolutionFiniteGroup(gensG,N+1,false,p);
return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyPrimePowerGroup:=function()
local
		R, L, LL, x;

if Order(G)<32 then R:=ResolutionFiniteGroup(gensG,N+1,false,p); 
else
#L:=LowerCentralSeries(G);
L:=BigStepLCS(G,9);

R:=ResolutionNormalSeries(L,N+1,false,p);
fi;

return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyAbelianGroup:=function()
local L,R;

L:=AbelianInvariantsToTorsionCoefficients(AbelianInvariants(G));
R:=ResolutionAbelianGroup(L,N+1);

return Homology(Functor(R),N);
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

H:=Homology(Functor(R),N);
if IsInt(H) then
	if H=0 then H:=[]; else H:=[H]; fi;
fi;

if 
Order(Centralizer(G,S))=Order(Center(S))*Order(G)/Order(S)
or 
Length(H)=0 
then Append(TorsionCoeffs,H); 
else
Append(TorsionCoeffs,
PrimePartDerivedFunctor(G,R,Functor,N));
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
return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HomologyNilpotentPcpGroup:=function()
local R;

R:=ResolutionNilpotentGroup(G,N+1);;
return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################



if IsList(D) then
return HomologyArtinGroup(); fi;

if IsPcpGroup(G) then
	if IsNilpotentGroup(G) then
	return HomologyNilpotentPcpGroup();
	else
	Print("Only nilpotent pcp groups are handled by this function. \n");
	return fail;
	fi;
fi;

if IsAbelian(G) then
return HomologyAbelianGroup(); fi;

if IsPrimePowerInt(Order(G)) then
return HomologyPrimePowerGroup(); fi;

if Order(G)<16 then
return HomologySmallGroup(); fi;

return HomologyGenericGroup(); 


end);
