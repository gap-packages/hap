#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionFiniteExtension,
function(arg)
local 
		GensE,GensG,R,n,tietze,
                G, EltsG,
		E, EltsE, Mult, InvE, 
		N, EltsN, GensN, GensNfirst, GensNsecond,
		EhomG, EhomGfirst, EhomGsecond,
		GmapE, GmapEfirst, GmapEsecond,
		NhomE, NhomEfirst, NhomEsecond,
		NEhomN, NEhomNfirst,
		S,
		k,p;

GensE:=StructuralCopy(arg[1]);

GensG:=StructuralCopy(arg[2]);
R:=arg[3];
n:=arg[4];
if Length(arg)>4 then tietze:=arg[5]; else tietze:=false; fi;
G:=R!.group;
EltsG:=R!.elts;
E:=Group(GensE);
EltsE:=Elements(E);

#####################################################################
Mult:=function(i,j);
return Position(EltsE,EltsE[i]*EltsE[j]);
end;
#####################################################################

#####################################################################
InvE:=function(i);
return Position(EltsE,EltsE[i]^-1);
end;
#####################################################################

EhomGfirst:=GroupHomomorphismByImagesNC(E,G,GensE,GensG);

if Order(E)<10^4 or n> 5 then
EhomGsecond:=List([1..Size(E)],i->Position(EltsG,Image(EhomGfirst,EltsE[i])));

#####################################################################
EhomG:=function(i);
return EhomGsecond[i];
end;
#####################################################################

else
EhomGsecond:=List([1..Size(E)],i->0);

#####################################################################
EhomG:=function(i);
if EhomGsecond[i]=0 then
EhomGsecond[i]:= Position(EltsG,Image(EhomGfirst,EltsE[i]));
fi;
return EhomGsecond[i];
end;
#####################################################################
fi;
#GmapEfirst:=GroupHomomorphismByImagesNC(G,E,GensG,GensE);
#GmapEsecond:=List([1..Size(G)],i->Position(EltsE,Image(GmapEfirst,EltsG[i])));

GmapEsecond:=List([1..Size(G)],i->Position(EltsE,
PreImagesRepresentative(EhomGfirst,EltsG[i])));

#####################################################################
GmapE:=function(i);
return GmapEsecond[i];
end;
#####################################################################

N:=Kernel(EhomGfirst);
if IsAbelian(N) then GensN:=TorsionGeneratorsAbelianGroup(N);
else 
GensN:=ReduceGenerators(GeneratorsOfGroup(N),N);
fi;

if Order(N)=1 then GensN:=[Identity(N)];
else
if Length(GensN) > 1 then
   GensNfirst:=StructuralCopy(GensN);
	for k in GensNfirst do
	GensNsecond:=SSortedList(GensN);
	RemoveSet(GensNsecond,k);
	if Order(Group(GensNsecond))=Order(N) then GensN:=GensNsecond; fi;
	od;
   fi;
fi;

if Length(arg)>5 then
S:=arg[6];
else
#if IsAbelian(N) then			#This should always work but it doesn't! 
#S:=ResolutionFiniteGroup(GensN,n,tietze);  #June 2022
S:=ResolutionGenericGroup(Group(GensN),n);
fi;

EltsN:=S!.elts;
NhomEfirst:=GroupHomomorphismByImagesNC(N,E,GensN,GensN);
#NhomEfirst:=GroupHomomorphismByFunction(N,E,x->x);

NhomEsecond:=List([1..Size(N)],i->Position(EltsE,Image(NhomEfirst,EltsN[i])));

#####################################################################
NhomE:=function(i);
return NhomEsecond[i];
end;
#####################################################################

NEhomNfirst:=List([1..Size(E)], k->Position(NhomEsecond,k));

#This next function can produce a fail when incorrectly used!!
#####################################################################
NEhomN:=function(i);
return NEhomNfirst[i];
end;
#####################################################################

return
TwistedTensorProduct(R,S,EhomG,GmapE,NhomE,NEhomN,EltsE,Mult,InvE);

end);
#####################################################################



