
###########################################################
###########################################################
InstallGlobalFunction(HeckeOperator,
function(arg)
local  H, pp, k, N, p, g, gg, HH, R, RH, TH, map1, THH, map3, bool,
ag, A, facts, T,E, Sfacts,Hecke ,fn, ConjugateIntersection,HHhomH; 
#TietzeReducedResolution;

H:=arg[1];
pp:=arg[2];
k:=arg[3];
#TietzeReducedResolution:=function(R); return R; end;
#There is a bug in the contracting homotopy of Tietze reducedResolution for
#degree >1 but we only use degree 1 here.

if IsMatrix(pp) then
return HeckeComponent(H,pp,k);
fi;

bool:=true;
if "name" in NamesOfComponents(H) then
if H!.name="CongruenceSubgroupGamma0" then bool:=false;
fi;
fi;
if bool then 
Print("Warning: This function is designed for \"CongruenceSubgroupGamma0\"\n"); 
return fail; fi;

N:=1;
if IsInt(pp) then
   facts:=Factors(pp);
   if Length(facts)=1 then 
   p:=facts[1];
   g:=[[1,0],[0,p]];
   else g:=fail;
   fi;
   if Gcd(pp,H!.level)>1 then 
   #Print("The second argument is not coprime to the level.\n");
   #return fail;
   fi;
fi;


########################
fn:=function(A,n);
if n=0 then return p^(k-2)*IdentityMat(Length(A)); fi;
if n=1 then return A; fi;
if n=2 then return A^2 - p^(k-1)*IdentityMat(Length(A)); fi;

return A*fn(A,n-1) - p^(k-1)*fn(A,n-2);

end;
########################

T:=[];
E:=[];
Hecke:=[];
Sfacts:=SSortedList(facts);
for p in Sfacts do
T[p]:=HeckeComponent(H,[[1,0],[0,p]],k);
E[p]:=Length(Filtered(facts,x->x=p));
if Gcd(p,H!.level)>1 then
Add(Hecke,  T[p]^E[p]);
else
Add(Hecke,  fn(T[p],E[p]));
fi;
od;

if Hecke[1]=[] then return []; fi;
return Product(Hecke);
end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(HeckeComponentWeight2,
function(H,p,N)
local g, gg, Hg, HH, gHH, R, RH, TH, RgHH, TgHH,
HHhomgHH, gHHhomH, map1, THH, eqmap2, map2, eqmap3, map3,
map32, map321,t;
#TietzeReducedResolution;

#TietzeReducedResolution:=function(R) return R; end;

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
if N=1 then
TH:=TietzeReducedResolution(RH);;
else
TH:=RH;
fi;
RgHH:=ResolutionFiniteSubgroup(TH,gHH);;
if N=1 then
TgHH:=TietzeReducedResolution(RgHH);;
else
TgHH:=RgHH;
fi;

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
function(t)
local A, MyExponents;

#######
MyExponents:=function(w)
local e;
e:=Exponents(w);
if e=[] then e:=[0]; fi;
return e;
end;
#######
A:=List(GeneratorsOfGroup(Source(t)),x->MyExponents(Image(t,x)));
if Length(A)=0 then A:=[[0]]; fi;
return A;
end);
###########################################################
###########################################################


###########################################################
###########################################################
InstallGlobalFunction(HeckeComponent,
function(arg)
local  H, pp, k, N, p, g, gg, HH, R, RH, TH, map1, THH, map3, 
ag, A, facts, T,E, Sfacts,Hecke ,fn, ConjugateIntersection,HHhomH;
#TietzeReducedResolution;

H:=arg[1];
pp:=arg[2];
k:=arg[3];
#TietzeReducedResolution:=function(R); return R; end;
#There is a bug in the contracting homotopy of Tietze reducedResolution for
#degree >1 but we only use degree 1 here.


N:=1;
   g:=pp; p:=0;
ag:=[[ g[2][2]  , -g[1][2] ],
     [ -g[2][1] , g[1][1]  ]];;

#######################################
ConjugateIntersection:=function(H,g)
local Hg, HH;

#THIS FUNCTION NEEDS TO BE PROPERLY IMPLEMENTED USING THE OBVIOUS
#COSET REPS FOR GAMMA0. THIS SHOULD THEN WORK MUST FASTER

Hg:=H^g;
HH:=Intersection(H,Hg);;
return HH;
end;
#######################################

#HH:=ConjugateIntersection(H,g);
HH:=HAP_IntersectionConjugatedCongruenceSubgroup(H,g);
gg:=g^-1;;


if Length(arg)<4 then
R:=ResolutionSL2Z_alt(N+1);;
R!.group:=HAP_CongruenceSubgroupGamma0(1);
R:=TietzeReducedResolution(R);;  #This just invokes a recording of homotopies.
RH:=ResolutionFiniteSubgroup(R,H);;
TH:=TietzeReducedResolution(RH);;
else
TH:=arg[4];
fi;

HHhomH:=GroupHomomorphismByFunction(HH,H,x->x^gg);

A:=HomogeneousPolynomials(H,k-2);

#Print(Runtime(),"\n");
map1:=TransferCochainMap(TH,HH,A);;
THH:=map1!.resH;   #resolution for HH
map1:=HomToIntegersModP(map1,-1/2);
#if false then

map1:=HomologyVectorSpace(map1,N);


map1:=List(Basis(Source(map1)),x->Image(map1,x));
map1:=TransposedMat(map1);
#Print(Runtime(),"\n");;
#fi;

map3:=EquivariantChainMap(THH,TH,HHhomH);
map3:=HomToIntegralModule(map3,A,ag);
map3:=HomToIntegersModP(map3,-1/2);
#return [map3,map1];

map3:=HomologyVectorSpace(map3,N);


map3:=List(Basis(Source(map3)),x->Image(map3,x));
map3:=TransposedMat(map3);
#Print(Runtime(),"\n");;


if map1=[] or map3=[] then return []; fi;
return Product([map3,map1]);
end);
#####################################################
#####################################################


