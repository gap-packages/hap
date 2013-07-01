#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(GroupHomology,
function(arg)
local
        GroupCohomologyOriginal, answer;

#####################################################################
GroupCohomologyOriginal:=function()
	local
		G, gensG, N, p, D, 
		Functor,
		HomologyGraphOfGroups,
		HomologyPrimePowerGroup,
		HomologyGenericGroup,
		HomologySmallGroup,
		HomologyCoxeterGroup,
		HomologyArtinGroup,
		HomologyAbelianGroup,
		HomologyNilpotentPcpGroup,
		HomologySpaceGroup;
		


############################### INPUT DATA ##########################
if IsList(arg[1]) then 
   if IsString(arg[1][1]) then D:=arg[1][2]; G:=false; 
   else D:=arg[1]; G:=false; fi; 
else
	if IsGroup(arg[1]) then G:=arg[1]; 
	if Order(G)<infinity then
	gensG:=ReduceGenerators(GeneratorsOfGroup(G),G); fi; 
	D:=false;
	else
	Print("ERROR: first variable must be a group or a Coxeter diagram or a graph of groups. \n");
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

if N=0 and p=0 then return [0]; fi;
if N=0 and p>0 then return [p]; fi;

############################## DATA INPUT ###########################

if IsPrime(p) then
Functor:=function(R); return TensorWithIntegersModP(R,p); end;
else
Functor:=TensorWithIntegers;
fi;

#####################################################################
#####################################################################
HomologyGraphOfGroups:=function()
local R;

R:=ResolutionGraphOfGroups(D,N+1);
return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################


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
HomologySpaceGroup:=function()
local R, G1;

if IsMatrixGroup(G) then
G1:=Image(IsomorphismPcpGroup(G)); fi;

R:=ResolutionAlmostCrystalGroup(G1,N+1);
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
L:=LowerCentralSeries(G);
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

if Order(S)>=128 or N>2 then
R:=ResolutionNormalSeries(LowerCentralSeries(S),N+1);
else
gens:=GeneratorsOfGroup(S);
gens:=ReduceGenerators(gens,S);
R:=ResolutionFiniteGroup(gens,N+1,false,p);
fi;
#R:=TietzeReducedResolution(R);
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
HomologyCoxeterGroup:=function()
local R;

if N+1>Length(GeneratorsOfGroup(CoxeterDiagramFpArtinGroup(D)[1])) then
Print("At present this function only works in dimensions < ",Length(GeneratorsOfGroup(CoxeterDiagramFpArtinGroup(D)[1]))," for the given Coxeter group.\n");
return fail;
fi;

R:=ResolutionCoxeterGroup(D,N+1);;
return Homology(Functor(R),N);

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
if GraphOfGroupsTest(D) then
return HomologyGraphOfGroups();
fi;
if IsString(arg[1][1])   then
   if ('c' in arg[1][1]) or ('C' in arg[1][1]) then
   return HomologyCoxeterGroup(); fi;
fi;
return HomologyArtinGroup();
fi;

if "CrystCatRecord" in KnownAttributesOfObject(G) or
   "AlmostCrystallographicInfo" in KnownAttributesOfObject(G) 
   then
return HomologySpaceGroup(); fi;

if IsPcpGroup(G) then
	if IsNilpotentGroup(G) then
	return HomologyNilpotentPcpGroup();
	else
	Print("Only nilpotent pcp groups are handled by this function. \n");
	return fail;
	fi;
fi;

if IsFinite(G) then

if IsAbelian(G) then
return HomologyAbelianGroup(); fi;

if IsPrime(p) and IsPrimePowerInt(Order(G)) and
Order(G)<257 then
if p=PrimePGroup(G) then
return List([1..RankPrimeHomology(G,N)(N)],i->p);
else
return [];
fi;
fi;

if IsPrimePowerInt(Order(G)) then
return HomologyPrimePowerGroup(); fi;

if Order(G)<16 then
return HomologySmallGroup(); fi;

return HomologyGenericGroup(); 
fi;

end;
#####################################################################

answer:= GroupCohomologyOriginal();

if IsList(answer) then return answer; fi;

if IsInt(answer) then return
ListWithIdenticalEntries(answer,arg[3]); fi;


end);
