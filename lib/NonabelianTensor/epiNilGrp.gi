#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(UpperEpicentralSeries,
function(G,cls)
local
	GhomFpG, FpGhomG, FpG, U, UU,F, rels, newrels, x, r,
	epi,upperZ,c,i, UUhomG, UhomG,L,epicentre;

if "EpiCentre" in NamesOfComponents(G) then
return G!.EpiCentre; fi;

if not IsNilpotent(G) then
Print("Group must be nilpotent \n"); return fail;
fi;
    
GhomFpG:=IsomorphismFpGroup(G);
FpGhomG:=InverseGeneralMapping(GhomFpG);
FpG:=Image(GhomFpG);
F:=FreeGroupOfFpGroup(FpG);
rels:=RelatorsOfFpGroup(FpG);

for i in [1..cls] do
newrels:=[];
for x in GeneratorsOfGroup(F) do
for r in rels do
c:=Comm(x,r);
if not c=Identity(F) then
Append(newrels,[c]); fi;
od;
od;
rels:=newrels;
od;

U:=F/rels;
if Length(GeneratorsOfGroup(F))>1 then
epi:=HAP_NqEpimorphismNilpotentQuotient(U);
else
epi:=HAP_NqEpimorphismNilpotentQuotient(U,1);
fi;
UU:=Image(epi);

upperZ:=Reversed(UpperCentralSeries(UU));
upperZ:=upperZ[cls+1];

if IsPcpGroup(G) then

UUhomG := GroupHomomorphismByImagesNC(
          UU,G,
          List(GeneratorsOfGroup(U), x->Image(epi,x)),
          List(GeneratorsOfGroup(FpG),x->Image(FpGhomG,x))
 );

else
UhomG:=GroupHomomorphismByImagesNC(U,G,GeneratorsOfGroup(U),
List(GeneratorsOfGroup(FpG),x->PreImagesRepresentative(GhomFpG,x)));

UUhomG:=GroupHomomorphismByFunction(UU,G,x->
Image(UhomG,PreImagesRepresentative(epi,x)));
fi;

epicentre:=Subgroup(G,List(GeneratorsOfGroup(upperZ),x->Image(UUhomG,x)));

G!.EpiCentre:=epicentre;
return epicentre;
end);
#####################################################################
