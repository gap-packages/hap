

###########################################################
###########################################################
InstallGlobalFunction(HeckeOperatorHigherWeight,
function(H,p,N,k)
local g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3,
map23, map123,t,A, AA,  TietzeReduced, hom;

TietzeReduced:=function(R) return R; end;
TietzeReduced:=TietzeReducedResolution;

#g:=[[1/p,0],[0,p]];;
g:=[[1,0],[0,p]];;   #This one seems to agree with Stein's book!!

Hg:=H^g;;
HH:=Intersection(H,Hg);;
gg:=g^-1;;
gHH:=HH^gg;;

R:=ResolutionSL2Z_alt(N+1);;
RH:=ResolutionFiniteSubgroup(R,H);;
TH:=TietzeReduced(RH);;
RgHH:=ResolutionFiniteSubgroup(TH,gHH);;
TgHH:=TietzeReduced(RgHH);;

HHhomgHH:=GroupHomomorphismByFunction(HH,gHH,x->x^gg);
gHHhomH:=GroupHomomorphismByFunction(gHH,H,x->x);

A:=HomogeneousPolynomials(H,k);
map1:=TransferCochainMap(TH,HH,A);;  #THH-->TH
THH:=map1!.resH;


eqmap2:=EquivariantChainMap(THH,TgHH,HHhomgHH);
map2:=HomToIntegralModule(eqmap2,A);  #TgHH-->THH  NEEDS FIXING HERE
eqmap3:=EquivariantChainMap(TgHH,TH,gHHhomH);
map3:=HomToIntegralModule(eqmap3,A);   #TH-->TgHH

map23:=Compose(map2,map3);  #TH-->THH
map123:=Compose(map1,map23); #TH-->TH


t:=Cohomology(map123,N);
return t;
t:=List(GeneratorsOfGroup(Source(t)),x->Exponents(Image(t,x)));
return [t,Homology(Source(map123),1)];
end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(HeckeOperatorWeight2,
function(H,p,N)
local g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3,
map32, map321,t, TietzeReduced;

TietzeReduced:=function(R) return R; end;

#g:=[[1/p,0],[0,p]];;
g:=[[1,0],[0,p]];;   #This one seems to agree with Stein's book!!

Hg:=H^g;;
HH:=Intersection(H,Hg);;
gg:=g^-1;;
gHH:=HH^gg;;

R:=ResolutionSL2Z_alt(N+1);;
RH:=ResolutionFiniteSubgroup(R,H);;
TH:=TietzeReduced(RH);;
RgHH:=ResolutionFiniteSubgroup(TH,gHH);;
TgHH:=TietzeReduced(RgHH);;

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

