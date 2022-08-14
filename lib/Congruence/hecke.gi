
###########################################################
###########################################################
InstallGlobalFunction(HeckeOperator,
function(H,pp,k)
local N,p, q, g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH, gH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3, map4,
ag, A, RR, facts,Q;

N:=1;
if IsInt(pp) then
   facts:=Factors(pp);
   p:=facts[1];
   g:=[[p,0],[0,1]];;
   if Length(facts)=1 then q:=1;
   else q:=Product(facts{[2..Length(facts)]});
   fi;
fi;
if IsMatrix(pp) then
   g:=pp; q:=1;
fi;

ag:=[[ g[2][2]  , -g[1][2] ],
     [ -g[2][1] , g[1][1]  ]];;   

Hg:=H^g;
HH:=Intersection(H,Hg);;
gg:=g^-1;;
gH:=H^gg;
gHH:=HH^gg;;

R:=ResolutionSL2Z_alt(N+1);;
RH:=ResolutionFiniteSubgroup(R,H);;
TH:=TietzeReducedResolution(RH);;
RgHH:=ResolutionFiniteSubgroup(TH,gHH);;
TgHH:=TietzeReducedResolution(RgHH);;

HHhomgHH:=GroupHomomorphismByFunction(HH,gHH,x->x^gg);
gHHhomH:=GroupHomomorphismByFunction(gHH,H,x->x);

A:=HomogeneousPolynomials(H,k-2);


if q=1 then
map1:=TransferCochainMap(TH,HH,A);;      
THH:=map1!.resH;   #resolution for HH
map1:=HomToIntegersModP(map1,-1/2);
map1:=HomologyVectorSpace(map1,N);
map1:=List(Basis(Source(map1)),x->Image(map1,x));
map1:=TransposedMat(map1);
else
Print("Second argument is not prime.\n"); return fail;
map1:=HeckeOperator(HAP_CongruenceSubgroupGamma0(p*H!.level),q,N,k);;
THH:=TransferCochainMap(TH,HH,A);;    
THH:=THH!.resH;
fi;


eqmap2:=EquivariantChainMap(THH,TgHH,HHhomgHH);
map2:=HomToIntegralModule(eqmap2,A,ag);  #TgHH-->THH  
map2:=HomToIntegersModP(map2,-1/2);
map2:=HomologyVectorSpace(map2,N);
map2:=List(Basis(Source(map2)),x->Image(map2,x));
map2:=TransposedMat(map2);


eqmap3:=EquivariantChainMap(TgHH,TH,gHHhomH);
map3:=HomToIntegralModule(eqmap3,A);   #TH-->TgHH
map3:=HomToIntegersModP(map3,-1/2);
map3:=HomologyVectorSpace(map3,N);
map3:=List(Basis(Source(map3)),x->Image(map3,x));
map3:=TransposedMat(map3);


return Product([map3,map2,map1]);

end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(HeckeOperatorWeight2,
function(H,p,N)
local g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3,
map32, map321,t;

if not IsPrimeInt(p) then
Print("Second argument is not prime.\n");
return fail;
fi;
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

