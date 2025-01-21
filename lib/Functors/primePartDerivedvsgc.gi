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
	charac, x, y, i, Cent, hh, HPpres, ord, Pone, RPone;


####################################
####################################
P:=R!.group;
prime:=PrimePGroup(P);
C:=F(R);

if IsGroup(GG) then G:=GG; 
P1:=Normalizer(G,P);

AscChn:=AscendingChain(G,P1 : refineIndex:=10);  #Added refineIndex, December 2024
fi;
if IsList(GG) then G:=GG[Length(GG)]; 
AscChn:=GG;
P1:=Normalizer(G,P);
fi;

charac:=EvaluateProperty(F(R),"characteristic");
if charac=0 then
x:=IntegralHomology("HomologyAsFpGroup",n);       #Modified December 2024
else
x:=ModularHomology("HomologyAsFpGroup",n);       #Modified December 2024
fi;
HP:=x(F(R),n);
HP:=HP.fpgroup;
if Length(AbelianInvariants(HP))=0 then return []; fi;
HPrels:=[Identity(HP)];
####################################
####################################


#########################################
#########################################
AddRels:=function(chn,L)  #chn=[P = Q1 < <Q2 < ... < Qk=Q] is a chain of 
                          #k subgroups in G with P=Syl_p(G)
local Q, NS, i, hh, Lhh, ggg, gg, g, gg1, h, sylQQ, QQ, xx, RC, gens, bool, Q0, Q1, S, SS;
Q:=chn[Length(chn)];
QQ:=Intersection(Q,Q^L); 
sylQQ:=SylowSubgroup(QQ,prime);
if not Order(sylQQ)>1 then return; fi;
gg:=One(G);

for i in [1..Length(chn)-1] do
################
################  Changed January 2025
S:=sylQQ^gg;
gens:=SmallGeneratingSet(S);
Q0:=chn[Length(chn)-i+1];
Q1:=chn[Length(chn)-i];
NS:=Normalizer(Q0,Q1);
RC:=RightTransversal(Q0,NS); 
RC:=List(RC,x->x^-1);

Unbind(ggg);
for g in RC do             
bool:=true;
   for xx in gens do
      if not xx^g in Q1 then bool:=false; break; fi;
   od;
if bool then ggg:=g; break; fi; 
od;
################
################
gg:=gg*ggg;

od;



#########################################
#########################################
if Order(P)/Order(sylQQ)>1 then  #NEED TO OPTIMIZETHIS CHOICE!!

S:=ResolutionGenericGroup(sylQQ,n+1);

else

#S:=ResolutionFiniteSubgroup(R,sylQQ^gg);   #WITH THIS!!
S:=R;

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

hh:=One(G);

for i in [1..Length(chn)-1] do
################
################  Changed January 2025
SS:=imfx^hh;
gens:=SmallGeneratingSet(SS);
Q0:=chn[Length(chn)-i+1];
Q1:=chn[Length(chn)-i];
NS:=Normalizer(Q0,Q1);
RC:=RightTransversal(Q0,NS);
RC:=List(RC,x->x^-1);

Unbind(ggg);
for h in RC do
bool:=true;
   for xx in gens do
      if not xx^h in Q1 then bool:=false; break; fi;
   od;
if bool then ggg:=h; break; fi;
od;
################
################
hh:=hh*ggg;

od;


#######################
#######################
#RC:=RightCosets(Q,Normalizer(Q,imfx));
#RC:=List(RC,Representative);
#for h in RC do                   #December 2024 changed from Q to RC
#if IsSubgroup(P,imfx^h) then hh:=h; break; fi;
#od;
#######################
#######################


Lhh:=L^-1*hh;
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
   AddRels([P],L);
   od;
fi;
for i in [2..Length(AscChn)] do

DCRS:=List(DoubleCosetRepsAndSizes(AscChn[i],AscChn[i-1],AscChn[i-1]),
x->x[1]);
Cent:=Centralizer(AscChn[i],AscChn[i-1]);


Sort(DCRS,ord);
DCRS:=Filtered(DCRS,a->not a in Cent);  #This does not achieve much
#DCRS:=Classify(DCRS,x->Cent*x);        #And this achieves nothing!
#DCRS:=List(DCRS,x->x[1]);              #

   for L in DCRS do
   cnt:=cnt+1;
   AddRels(AscChn{[1..i-1]},L);
   od;
od;
#########################################
#########################################



return AbelianInvariants(HP/NormalClosure(HP,Group(HPrels)));
end);
#####################################################################

