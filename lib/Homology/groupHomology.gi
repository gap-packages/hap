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
		Functor,Functor2,
		HomologyGraphOfGroups,
		HomologyPrimePowerGroup,
		HomologyGenericGroup,
		HomologySmallGroup,
		HomologyCoxeterGroup,
		HomologyArtinGroup,
		HomologyAbelianGroup,
		HomologyNilpotentPcpGroup,
		HomologySpaceGroup,
                HomologyAlmostCrystallographicGroup,
                HomologyOfGroupHomomorphism;
		


############################### INPUT DATA ##########################
if IsList(arg[1]) then 
   if IsString(arg[1][1]) then D:=arg[1][2]; G:=false; 
   else D:=arg[1]; G:=false; fi; 
else
	if IsGroup(arg[1]) then G:=arg[1]; 
	if Order(G)<infinity then
        if IsMatrixGroup(G) then G:=Image(IsomorphismPermGroup(G),G); fi;
	gensG:=ReduceGenerators(GeneratorsOfGroup(G),G); fi; 
	D:=false;
	else
        if IsGroupHomomorphism(arg[1]) then D:=false; G:=false;
        else
	Print("ERROR: first variable must be a group or a Coxeter diagram or a graph of groups or a group homomorphism. \n");
        fi;
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

if N=0 and p=0 and not IsGroupHomomorphism(arg[1]) then return [0]; fi;
if N=0 and p>0 and not IsGroupHomomorphism(arg[1]) then return [p]; fi;

############################## DATA INPUT ###########################

if IsPrime(p) then
Functor:=function(R); return TensorWithIntegersModP(R,p); end;
else
Functor:=function(R); return ContractedComplex(TensorWithIntegers(R)); end;
fi;

if IsPrime(p) then
Functor2:=function(R); return TensorWithIntegersModP(R,p); end;
else
Functor2:=TensorWithIntegers;
fi;

#####################################################################
#####################################################################
HomologyOfGroupHomomorphism:=function()
local f,H,G,RH,RG,PH,PG,ans,hlgy,k,A;
f:=arg[1];
H:=Source(f);
G:=Target(f);

if (not IsFinite(H)) or (not IsFinite(G)) then 
Print("This function is currently implemented only for homomorphisms of finite groups.\n");
return fail;
fi;

hlgy:=GroupHomology(H,N,p);
if Length(hlgy)=0 then 
A:=AbelianGroup(GroupHomology(G,N,p));
return GroupHomomorphismByFunction(SmallGroup(1,1),A,x->One(A)); fi;
hlgy:=Product(hlgy);
hlgy:=Factors(hlgy);
hlgy:=SSortedList(hlgy);
hlgy:=Filtered(hlgy,i->not i=1);
ans:=[];
for k in hlgy do
PH:=SylowSubgroup(H,k);
PG:=SylowSubgroup(G,k);
RG:=ResolutionGenericGroup(PG,N+1);
RH:=ResolutionGenericGroup(PH,N+1);
Add(ans, PrimePartDerivedFunctorHomomorphism(f,RG,RH,Functor2,N) );
od;
if Length(ans)>0 then ans:=DirectProduct(ans); fi;
return ans;
end;
#####################################################################
#####################################################################

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
HomologyAlmostCrystallographicGroup:=function(G,N)
local R ;

R:=ResolutionAlmostCrystalGroup(G,N+1);
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
if N<=1 then
R:=ResolutionNilpotentGroup(G,N+1);
else
L:=LowerCentralSeries(G);
L:=BigStepLCS(G,6);    #changed from 9, 22/12/2022

R:=ResolutionNormalSeries(L,N+1,false,p);
fi;
fi;

return Homology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyAbelianGroup:=function()
local L,R;

#L:=AbelianInvariantsToTorsionCoefficients(AbelianInvariants(G));
#R:=ResolutionAbelianGroup(L,N+1);

R:=ResolutionAbelianGroup(G,N+1);

return Homology(Functor(R),N);
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyGenericGroup:=function()
local
	primes, q, S, gens, functor, R, H, TorsionCoeffs;

if IsPrime(p) then
S:=SylowSubgroup(G,p);
R:=ResolutionPrimePowerGroup(S,N+1);
return PrimePartDerivedFunctor(G,R,Functor2,N);
fi;



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
PrimePartDerivedFunctorViaSubgroupChain(G,R,Functor2,N));
#PrimePartDerivedFunctor(G,R,Functor,N));
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

if IsGroupHomomorphism(arg[1]) then 
return HomologyOfGroupHomomorphism(); fi;

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

#if "CrystCatRecord" in KnownAttributesOfObject(G) or
#   "AlmostCrystallographicInfo" in KnownAttributesOfObject(G) 
#   then
#return HomologySpaceGroup(); fi;

if IsPcpGroup(G) then
        if IsAlmostCrystallographic(G) then 
        return HomologyAlmostCrystallographicGroup(G,N);
        fi;
	if IsNilpotentGroup(G) then
	return HomologyNilpotentPcpGroup();
	else
	Print("Only nilpotent pcp groups are handled by this function. \n");
	return fail;
	fi;
fi;

if "CrystCatRecord" in KnownAttributesOfObject(G) or
   "AlmostCrystallographicInfo" in KnownAttributesOfObject(G)
   then
return HomologySpaceGroup(); fi;

if IsFinite(G) then

if p>0 and not p in Factors(Order(G)) then return []; fi;

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

if not IsInt(answer) then return answer; fi;

if IsInt(answer) then return
ListWithIdenticalEntries(answer,arg[3]); fi;


end);
#####################################################################
#####################################################################

#######################################################
#######################################################
InstallGlobalFunction(RelativeGroupHomology,
function(G,N,kk)
local i, GhomQ, Q, RQ, RN, T, C, D, K,
      newdimension, newboundary, newproperties;

GhomQ:=NaturalHomomorphismByNormalSubgroup(G,N);
Q:=Range(GhomQ);
RQ:=ResolutionFiniteGroup(Q,kk+1);
RN:=ResolutionFiniteGroup(N,kk+1);
T:=ResolutionExtension(GhomQ,RN,RQ);

C:=TensorWithIntegers(T);
D:=TensorWithIntegers(RQ);

######
newdimension:=function(n);
return C!.dimension(n)-D!.dimension(n);
end;
######

######
newboundary:=function(n,k);
return C!.boundary(n,k+D!.dimension(n));
end;
######

######
newproperties:=[["length",kk+1],
                ["connected",true],
                ["type","chainComplex"],
                ["characteristic",0]];
######

K:=rec(dimension:=newdimension,
       boundary:=newboundary,
       properties:=newproperties);
K:=Objectify(HapChainComplex,K);

return Homology(K,kk);
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(DirectProductOfGroupHomomorphisms,
function(f,g)
local S,T;
S:=DirectProduct(Source(f),Source(g));
T:=DirectProduct(Target(f),Target(g));
return GroupHomomorphismByFunction(S,T,x->
Image(Embedding(T,1), Image(f,Image(Projection(S,1),x)))
*
Image(Embedding(T,2),Image(g,Image(Projection(S,2),x))) 
);
end);
#######################################################
#######################################################
