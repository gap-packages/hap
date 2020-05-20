#(C) Graham Ellis, 2005-2006
cnt:=0;
#####################################################################
InstallGlobalFunction(PrimePartDerivedFunctorViaSubgroupChain,
function(GG,R,F,n)
local
	G,C,P,P1, prime, AscChn, HP, HPrels, AddRels, Q,
        DCRS, L, S, f,fx, imfx, bool, dcrs,
	HK, HPK, HKhomHPK, HPKhomHP, HKhomHP, HKx,HPKx, 
	HKxhomHPKx, HPKxhomHP, HKxhomHP, HKhomHKx,  HKhomHP2,
	x, y, i, Cent, hh, HPpres, ord;


####################################
####################################
P:=R!.group;
prime:=PrimePGroup(P);
C:=F(R);

if IsGroup(GG) then G:=GG; 
P1:=Normalizer(G,P);

AscChn:=AscendingChain(G,P1);
fi;
if IsList(GG) then G:=GG[Length(GG)]; 
AscChn:=GG;
P1:=Normalizer(G,P);

fi;


HP:=GroupHomomorphismByFunction(P,P,x->x);
HP:=EquivariantChainMap(R,R,HP);
HP:=F(HP);
HP:=Homology(HP,n);   #This takes too much time!!

HP:=Source(HP);
if Length(AbelianInvariants(HP))=0 then return []; fi;
HPrels:=[Identity(HP)];
####################################
####################################


#########################################
#########################################
AddRels:=function(Q,L)  #Here P < Q < G where P=Syl_p(G)
local i, hh, Lhh, gg, g, gg1, h, sylQQ, QQ, xx;

QQ:=Intersection(Q,Q^L); 
sylQQ:=SylowSubgroup(QQ,prime);

gg:=false;
for g in Q do
if IsSubgroup(P,sylQQ^g) then gg:=g;  break; fi;
od;

if not Order(sylQQ)>1 then return; fi;

#########################################
#########################################
if Order(P)/Order(sylQQ)>1 then  #NEED TO OPTIMIZETHIS CHOICE!!

S:=ResolutionGenericGroup(sylQQ,n+1);

else

S:=ResolutionFiniteSubgroup(R,sylQQ^gg);

S!.group:=sylQQ;
gg1:=gg^-1;
S!.elts:=List(S!.elts,x->x^(gg1));
fi;
#########################################
#########################################

hh:=Homology(F(S),n);
if IsInt(hh) then hh:=List([1..hh],i->0); fi;
if not Length(hh)>0 then return; fi;
f:=GroupHomomorphismByFunction(sylQQ,P,x->x^gg);
xx:=F(EquivariantChainMap(S,R,f));;
HKhomHPK:=Homology(xx,n);
HK:=Source(HKhomHPK);
HPK:=Range(HKhomHPK);
HPKhomHP:=GroupHomomorphismByImagesNC(HPK,HP,GeneratorsOfGroup(HPK),
                                                  GeneratorsOfGroup(HP));
HKhomHP:=GroupHomomorphismByFunction(HK,HP,x->
Image(HPKhomHP, Image(HKhomHPK,x) ) );


fx:=GroupHomomorphismByFunction(sylQQ,Q,g->g^(L^-1));
imfx:=Image(fx);
hh:=false;
for h in Q do
if IsSubgroup(P,imfx^h) then hh:=h; break; fi;
od;
Lhh:=L^-1*hh;
#fx:=GroupHomomorphismByFunction(sylQQ,P,g->(g^(L^-1))^hh);
fx:=GroupHomomorphismByFunction(sylQQ,P,g->g^Lhh);

xx:=F(EquivariantChainMap(S,R,fx));
HKxhomHPKx:=Homology(xx,n);
HKx:=Source(HKxhomHPKx);
HPKx:=Parent(Image(HKxhomHPKx));
HPKxhomHP:=GroupHomomorphismByImagesNC(HPKx,HP,GeneratorsOfGroup(HPKx),
                                                  GeneratorsOfGroup(HP));
HKxhomHP:=GroupHomomorphismByFunction(HKx,HP,x->
Image(HPKxhomHP, Image(HKxhomHPKx,x) ) );
HKhomHKx:=GroupHomomorphismByImagesNC(HK,HKx,GeneratorsOfGroup(HK),GeneratorsOfGroup(HKx));
HKhomHP2:=GroupHomomorphismByFunction(HK,HP,a->
Image(HKxhomHP, Image(HKhomHKx,a)));

for x in GeneratorsOfGroup(HK) do
Append(HPrels, [Image(HKhomHP,x)*Image(HKhomHP2,x)^-1]);
od;

end;
#############################################
#############################################

####################################
####################################

ord:=function(x,y); return Order(x)<Order(y); end; 
if Order(P1)>Order(P) then 
DCRS:=SmallGeneratingSet(P1);
   for L in DCRS do
   AddRels(P,L);
   od;
fi;
for i in [2..Length(AscChn)] do

DCRS:=List(DoubleCosetRepsAndSizes(AscChn[i],AscChn[i-1],AscChn[i-1]),
x->x[1]);
Cent:=Centralizer(AscChn[i],AscChn[i-1]);


Sort(DCRS,ord);
DCRS:=Filtered(DCRS,a->not a in Cent);

############################################
#dcrs:=[];
#for x in DCRS do
#bool:=true;
#for y in dcrs do
#if y^-1*x in Cent or y^-1*x in P or x in Group(y) then bool:=false; break; fi;
#od;
#if bool then Add(dcrs,x); fi;
#od;
#for L in dcrs do  #This causes problems.  Need to think more!!!
############################################
   for L in DCRS do
   cnt:=cnt+1;
   AddRels(AscChn[i-1],L);
   od;
od;
#########################################
#########################################



return AbelianInvariants(HP/NormalClosure(HP,Group(HPrels)));
end);
#####################################################################

