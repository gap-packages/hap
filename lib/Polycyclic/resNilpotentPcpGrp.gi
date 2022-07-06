#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionNilpotentGroup,
function(arg)
local  
	EE,K,Test,
	N,E,G, 
	EEhomGG,
	EhomG,
	GmapE,
	RN,
	RG,
	RE,
	T;

SetInfoLevel(InfoWarning,0);
HAPconstant:=5;              #THIS CAUSES PROBLEMS !!!
EE:=arg[1];
K:=arg[2];
if Length(arg)>2 then Test:=arg[3]; else Test:="NoTest"; fi;
if IsAbelian(EE) then 
###                return ResolutionAbelianGroup(EE,K);
if IsPcpGroup(EE) then 
		return ResolutionAbelianPcpGroup(EE,K);
	else
		return ResolutionAbelianGroup_alt(EE,K);
	fi;
fi;

N:=Center(EE);

if IsPcpGroup(EE) then
	EEhomGG:=NaturalHomomorphism(EE,N);
else
	EEhomGG:=NaturalHomomorphismByNormalSubgroupNC(EE,N); 
fi;

E:=Source(EEhomGG);
G:=Range(EEhomGG);			#####here
if Test="TestFiniteness" then
RN:=ResolutionFiniteGroup(N,K);
else
RN:=ResolutionNilpotentGroup(N,K);;
fi;
RG:=ResolutionNilpotentGroup(G,K);;


T:= ResolutionExtension(EEhomGG,RN,RG,Test);
HAPconstant:=2;
SetInfoLevel(InfoWarning,1);
return T;
end);
#####################################################################

