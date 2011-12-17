#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionNilpotentGroup,
function(arg)
local  
	EE,K,Test,
	N,E,G, 
	EEhomGG,
	NhomE,
	EhomG,
	GmapE,
	RN,
	RG,
	RE;

EE:=arg[1];
K:=arg[2];
if Length(arg)>2 then Test:=arg[3]; else Test:="NoTest"; fi;
if IsAbelian(EE) then 
	if IsPcpGroup(EE) then 
		return ResolutionAbelianPcpGroup(EE,K);
	else
		return ResolutionAbelianGroup(EE,K);
	fi;
fi;

N:=Center(EE);

if IsPcpGroup(EE) then
	EEhomGG:=NaturalHomomorphism(EE,N);
else
	EEhomGG:=NaturalHomomorphismByNormalSubgroupNC(EE,N); 
fi;

E:=Source(EEhomGG);
G:=Image(EEhomGG);
if Test="TestFiniteness" then
RN:=ResolutionFiniteGroup(N,K);
else
RN:=ResolutionNilpotentGroup(N,K);;
fi;
RG:=ResolutionNilpotentGroup(G,K);;

return ResolutionExtension(EEhomGG,RN,RG,Test);
end);
#####################################################################

