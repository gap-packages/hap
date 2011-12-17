#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionSubnormalSeries,
function(arg)
local  
	L,Ltrunc,K,RN,RQ,N,E,Q,EhomQ,gensE,gensQ;

L:=arg[1];
K:=arg[2];

if Length(L)=1 then return
ResolutionFiniteGroup(L[1],K);
fi;

Ltrunc:=List([1..Length(L)-1],i->L[i]);

RN:=ResolutionSubnormalSeries(Ltrunc,K);
N:=RN!.group;
E:=L[1];
EhomQ:=NaturalHomomorphismByNormalSubgroup(E,N);
Q:=Image(EhomQ);
gensE:=GeneratorsOfGroup(E);
gensQ:=List(gensE,x->Image(EhomQ,x));

RQ:=ResolutionFiniteGroup(gensQ,K);


return ResolutionFiniteExtension(gensE,gensQ,RQ,K,false,RN);
end);
#####################################################################

