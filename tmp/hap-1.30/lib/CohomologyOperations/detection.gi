####################################################
####################################################
InstallGlobalFunction(InducedSteenrodHomomorphisms,
function(f,N)
local G,H,AG,AH, RG, RH, eqmap, BasH, BasG, VH, VG, hom, MAT, gens, imgs, M, i, r, s, x;

G:=Range(f);
H:=Source(f);
RG:=ResolutionPrimePowerGroup(G,N);
AG:=ModPCohomologyRing(RG);;
AH:=Mod2SteenrodAlgebra(H,N);
RH:=AH!.res;

eqmap:=EquivariantChainMap(RH,RG,f);

BasH:=[];
VH:=[];
BasG:=[];
VG:=[];
hom:=[];

for i in [0..N] do
BasH[i+1]:=Filtered(Basis(AH),x->AH!.degree(x)=i);
VH[i+1]:=Subspace(AH,BasH[i+1]);
BasG[i+1]:=Filtered(Basis(AG),x->AG!.degree(x)=i);
VG[i+1]:=Subspace(AG,BasG[i+1]);

MAT:=NullMat(RH!.dimension(i),RG!.dimension(i));;
for r in [1..RH!.dimension(i)] do
M:=eqmap!.mapping([[r,1]],i);
M:=List(M,x->x[1]);
M:=Collected(M);
for x in M do
MAT[r][x[1]]:=x[2] mod 2;
od;
od;

imgs:=List(BasH[i+1],x->x)*MAT;

gens:=BasG[i+1];
hom[i+1]:=LeftModuleGeneralMappingByImages(VG[i+1],VH[i+1], gens,imgs);
od;

return [AG, AH, hom];

end);
####################################################
####################################################



