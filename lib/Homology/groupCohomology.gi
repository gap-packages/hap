#(C) Graham Ellis, 2007

#####################################################################
#####################################################################
InstallGlobalFunction(GroupCohomology,
function(arg)
local
	GroupCohomologyOriginal, answer;

#####################################################################
GroupCohomologyOriginal:=function()
local
		G, gensG, N, p, D, 
		Functor,
		CoHomologyGraphOfGroups,
		CoHomologyArtinGroup,
		CoHomologyNilpotentPcpGroup,
		CoHomologySpaceGroup;


############################### INPUT DATA ##########################
if IsList(arg[1]) then D:=arg[1]; G:=false; 
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

############################## DATA INPUT ###########################

if IsPrime(p) then
Functor:=function(R); return HomToIntegersModP(R,p); end;
else
Functor:=HomToIntegers;
fi;

#####################################################################
#####################################################################
CoHomologyGraphOfGroups:=function()
local R;

R:=ResolutionGraphOfGroups(D,N+1);
return Cohomology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
CoHomologySpaceGroup:=function()
local R, G1;

if IsMatrixGroup(G) then
G1:=Image(IsomorphismPcpGroup(G)); fi;

R:=ResolutionAlmostCrystalGroup(G1,N+1);
return Cohomology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
CoHomologyArtinGroup:=function()
local R;

R:=ResolutionArtinGroup(D,N+1);;
return Cohomology(Functor(R),N);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
CoHomologyNilpotentPcpGroup:=function()
local R;

R:=ResolutionNilpotentGroup(G,N+1);;
return Cohomology(Functor(R),N);

end;
#####################################################################
#####################################################################

if IsList(D) then
if GraphOfGroupsTest(D) then
return CoHomologyGraphOfGroups();
else
return CoHomologyArtinGroup(); 
fi;
fi;

if "CrystCatRecord" in KnownAttributesOfObject(G) or
   "AlmostCrystallographicInfo" in KnownAttributesOfObject(G) 
   then
return CoHomologySpaceGroup(); fi;

if IsPcpGroup(G) then
	if IsNilpotentGroup(G) then
	return CoHomologyNilpotentPcpGroup();
	else
	Print("Only nilpotent pcp groups are handled by this function. \n");
	return fail;
	fi;
fi;

if IsFinite(G) then

if p=0 and N=0 then return [0];fi;
if p=0 and N=1 then return []; fi;
if p=0 then return GroupHomology(G,N-1); fi;
if IsPrimeInt(p) then return
#Length(Filtered(GroupHomology(G,N-1),i->GCD_INT(p,i)>1))  +
#Length(Filtered(GroupHomology(G,N),i->GCD_INT(p,i)>1));
GroupHomology(G,N,p);
fi;

fi;

end;
#####################################################################

answer:= GroupCohomologyOriginal();

if IsList(answer) then return answer; fi;

if IsInt(answer) then return
ListWithIdenticalEntries(answer,arg[3]); fi;

end);
