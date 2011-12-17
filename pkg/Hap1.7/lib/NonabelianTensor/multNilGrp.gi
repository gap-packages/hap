#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(BaerInvariant,
function(G,cls)
local
	GhomFpG, FpGhomG, FpG, U, UU,F, rels, newrels, x, r,
	epi,lowerUU,c,i, UUhomG, UhomG,L;

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
epi:=NqEpimorphismNilpotentQuotient(U);
UU:=Image(epi);

lowerUU:=LowerCentralSeries(UU);
lowerUU:=lowerUU[cls+1];


if IsPcpGroup(G) then

UUhomG := GroupHomomorphismByImages(
          UU,G,
          List(GeneratorsOfGroup(U), x->Image(epi,x)),
          List(GeneratorsOfGroup(FpG),x->Image(FpGhomG,x))
 );

 return AbelianInvariants(Intersection(Kernel(UUhomG), lowerUU));


else
UhomG:=GroupHomomorphismByImagesNC(U,G,GeneratorsOfGroup(U),
List(GeneratorsOfGroup(FpG),x->PreImagesRepresentative(GhomFpG,x)));

UUhomG:=GroupHomomorphismByFunction(lowerUU,G,x->
Image(UhomG,PreImagesRepresentative(epi,x)));
fi;


return AbelianInvariants(Kernel(UUhomG));
end);
#####################################################################
