
#################################################
#################################################
InstallGlobalFunction(CrossedInvariant,
function(G,CC)
local C, eta, pi1, pi2, homs, liftedhoms, delta, gensF, gensG, imgens, 
S,T, M, L, P, F, f, ff, r, m, x,g,Nf, Md, relsG, cnt, cnthom, cntrel, cntden,U;

if not IsFpGroup(G) then 
Print("First argument must be an fp group.\n");
return fail;
fi;

if not (IsHapCatOneGroup(CC) or IsHapCrossedModule(CC)) then
Print("Second argument must be a crossed module or a cat-1-group.\n");
return fail;
fi;


if IsHapCrossedModule(CC) then 
C:=CatOneGroupByCrossedModule(CC);;
else
C:=CC;
fi;

S:=C!.sourceMap; #S:U-->U
T:=C!.targetMap; #T:U-->U
U:=Source(S);
P:=Image(S);
M:=Kernel(S);
delta:=GroupHomomorphismByFunction(M,P,x->Image(T,x));
eta:=NaturalHomomorphismByNormalSubgroup(P,Image(delta));
pi1:=Target(eta);
pi2:=HomotopyGroup(C,2);

homs:=AllHomomorphisms(G,pi1);
liftedhoms:=[];
F:=FreeGroupOfFpGroup(G);
gensF:=GeneratorsOfGroup(F);
gensG:=GeneratorsOfGroup(G);
relsG:=RelatorsOfFpGroup(G);
cnt:=0;

Md:=List([1..Length(gensF)],i->pi2);
Md:=Cartesian(Md);


for f in homs do
imgens:=List(gensG, x->Image(f,x));
imgens:=List(imgens,x->PreImagesRepresentative(eta,x));
ff:=GroupHomomorphismByImages(F,P,gensF,imgens);
Add(liftedhoms,ff);
od;


for f in liftedhoms do
######################
cnthom:=[];
for r in relsG do
cntrel:=0;
for m in M do
if Image(delta,m)=Image(f,r) then cntrel:=cntrel+1;  fi;
od;
Add(cnthom,cntrel);
od;
Nf:=Product(cnthom);# number of morphisms G-->C over f:F-->P

cntden:=0;
for x in Md do
imgens:=List([1..Length(gensF)],i->x[i]*Image(f,gensF[i]));
g:=GroupHomomorphismByImages(F,U,gensF,imgens);
L:=List(relsG,r->Image(g,r));
L:=List(L,y->y*Image(S,y)^-1);
L:=List(L,Order);
if SSortedList(L)=[1] then  cntden:=cntden+1;  fi;
od;
if not IsInt(Nf/cntden) then Print("PROBLEM\n\n"); fi;
cnt:=cnt+(Nf/cntden);
######################
od;


return  cnt;

end);
#################################################
#################################################

