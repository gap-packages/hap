#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(BaerInvariant,
function(arg)
local
	G,cls,GhomFpG, FpGhomG, FpG, U, UU,F, rels, newrels, x, r,
	epi,lowerUU,c,i, UUhomG, UhomG,L,bool;

G:=arg[1];
cls:=arg[2];
if Length(arg)>2 then bool:=true;
else bool:=false; fi;

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
epi:=HAP_NqEpimorphismNilpotentQuotient(U);
UU:=Image(epi);

lowerUU:=LowerCentralSeries(UU);
lowerUU:=lowerUU[cls+1];


if IsPcpGroup(G) then

UUhomG := GroupHomomorphismByImagesNC(
          UU,G,
          List(GeneratorsOfGroup(U), x->Image(epi,x)),
          List(GeneratorsOfGroup(FpG),x->Image(FpGhomG,x))
 );

if not bool then
 return AbelianInvariants(Intersection(Kernel(UUhomG), lowerUU));
else
 return
 List(GeneratorsOfGroup(TorsionSubgroupPcpGroup(Intersection(Kernel(UUhomG), lowerUU))),
 x->[PreImagesRepresentative(epi,x),Order(x)]);
 fi;

else
UhomG:=GroupHomomorphismByImagesNC(U,G,GeneratorsOfGroup(U),
List(GeneratorsOfGroup(FpG),x->PreImagesRepresentative(GhomFpG,x)));

UUhomG:=GroupHomomorphismByFunction(lowerUU,G,x->
Image(UhomG,PreImagesRepresentative(epi,x)));
fi;


return AbelianInvariants(Kernel(UUhomG));
end);
#####################################################################
