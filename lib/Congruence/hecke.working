
###########################################################
###########################################################
InstallGlobalFunction(HeckeOperatorWeight2,
function(H,p,N)
local g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3,
map32, map321,t;

#g:=[[1/p,0],[0,p]];;
g:=[[1,0],[0,p]];;   #This one seems to agree with Stein's book!!

Hg:=H^g;;
HH:=Intersection(H,Hg);;
gg:=g^-1;;
gHH:=HH^gg;;

R:=ResolutionSL2Z_alt(N+1);;
RH:=ResolutionFiniteSubgroup(R,H);;
TH:=TietzeReducedResolution(RH);;
RgHH:=ResolutionFiniteSubgroup(TH,gHH);;
TgHH:=TietzeReducedResolution(RgHH);;
#TgHH:=RgHH;

HHhomgHH:=GroupHomomorphismByFunction(HH,gHH,x->x^gg);
gHHhomH:=GroupHomomorphismByFunction(gHH,H,x->x);

map1:=TransferChainMap(TH,HH);;
THH:=map1!.resH;


eqmap2:=EquivariantChainMap(THH,TgHH,HHhomgHH);
map2:=TensorWithIntegers(eqmap2);
eqmap3:=EquivariantChainMap(TgHH,TH,gHHhomH);
map3:=TensorWithIntegers(eqmap3);

map32:=Compose(map3,map2);
map321:=Compose(map32,map1);
map321:=HomToIntegers(map321);

t:=Cohomology(map321,N);
return t;
t:=List(GeneratorsOfGroup(Source(t)),x->Exponents(Image(t,x)));
return [t,Homology(Source(map321),1)];
end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(HomomorphismAsMatrix,
function(t);
return
List(GeneratorsOfGroup(Source(t)),x->Exponents(Image(t,x)));
end);
###########################################################
###########################################################

