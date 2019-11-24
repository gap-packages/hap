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

if not Order(L[Length(L)])=1 then     #
Add(L,Group(Identity(L[1]))); fi;     # 
Ltrunc:=List([2..Length(L)],i->L[i]); #changed here 7/03/2007

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

