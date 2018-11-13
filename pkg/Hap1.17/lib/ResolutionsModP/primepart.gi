#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(PoincareSeriesPrimePart,
function(G,p,n)
local
	PrimePart,
	PSeries,
	Resolutions,ResIndex,
        EqChnMps,EqChnMpsIndex,
	DCRSpruned;
ResIndex:=[];
Resolutions:=[];
EqChnMps:=[];
EqChnMpsIndex:=[];

#####################################################################
PrimePart:=
function(G,R,F,n)
local
	C,P, DCRS, DCRS1, L,Y,GroupL,
	X, K, gensK, S, f,fx, P1,
	HP, HK, HPK, HKhomHPK, HPKhomHP, HKhomHP,
	HKx,HPKx, 
	HKxhomHPKx, HPKxhomHP, HKxhomHP, HKhomHKx,  HKhomHP2,
	HPrels, x, y, i,prime, core, conjs, conjelt,CentP,
	HPpres,G1,epi,HPP,
        Pos;

HP:=HomologyVectorSpace(F(R),n);
HPrels:=[Zero(HP)];
if Dimension(HP)=0 then return []; fi;

P:=StructuralCopy(R!.group);
prime:=Factors(Order(P))[1];

if not IsBound(DCRSpruned) then
###############################################
P1:=Normalizer(G,P);
core:=[];
for x in P do
if Order(x)=prime then AddSet(core,x); fi;
od;

if not IsPrimeInt(Order(P)) then
DCRS1:=List(DoubleCosetRepsAndSizes(G,P1,P1),x->x[1]);
else
DCRS1:=[];
fi;
if Order(P1)>Order(P) then
Append(DCRS1,Filtered(ReduceGenerators(GeneratorsOfGroup(P1),P1),
x->not x in P));
fi;

DCRS:=[];
for x in DCRS1 do
for y in core do
if x*y*x^-1 in core then Append(DCRS,[x]); break; fi;
od;od;
DCRSpruned:=[];

conjs:=[];
conjelt:=[];
for x in DCRS do
Y:=Intersection(P,P^x);
AddSet(conjs,Y);
Append(conjelt,[[x,Y]]);
od;

for Y in conjs do
L:=Filtered(conjelt,x->x[2]=Y);
L:=List(L,x->x[1]);
GroupL:=Group(L);
Append(DCRSpruned,[ReduceGenerators(L,GroupL)]); 
od;
DCRSpruned:=Filtered(DCRSpruned,x->Length(x)>0);
############################################################
fi;


for L in DCRSpruned do
K:=Intersection(P,P^L[1]);
gensK:=ReduceGenerators(GeneratorsOfGroup(K),K);

G1:=Group(Concatenation(gensK,[Identity(P)]));

Pos:=Position(ResIndex,G1);
if Pos=fail then
S:=ResolutionPrimePowerGroup(G1,n+1);
Add(ResIndex,G1);
Add(Resolutions,S);
Pos:=Length(ResIndex);
fi;
S:=Resolutions[Pos];

if not Dimension(HomologyVectorSpace(F(S),n))=0 then

f:=GroupHomomorphismByFunction(K,P,x->x);

Pos:=Position(EqChnMpsIndex,[IdGroup(S!.group),IdGroup(R!.group),f]);
if Pos=fail then
Add(EqChnMpsIndex,[IdGroup(S!.group),IdGroup(R!.group),f]);
Add(EqChnMps,F(ModularEquivariantChainMap(S,R,f)));
Pos:=Length(EqChnMps);
fi;
HKhomHPK:=HomologyVectorSpace(EqChnMps[Pos],n);

HK:=Source(HKhomHPK);

HPK:=Range(HKhomHPK);

HPKhomHP:=LeftModuleGeneralMappingByImages(HPK,HP,
#HPKhomHP:=LeftModuleHomomorphismByImagesNC(HPK,HP,
Elements(Basis(HPK)),
Elements(Basis(HP)));

HKhomHP:=LeftModuleGeneralMappingByImages(HK,HP,
#HKhomHP:=LeftModuleHomomorphismByImagesNC(HK,HP,
Elements(Basis(HK)),
List(Elements(Basis(HK)),x->
ImagesRepresentative(HPKhomHP, ImagesRepresentative(HKhomHPK,x) )) );

for X in L do
fx:=GroupHomomorphismByFunction(K,P,g->Image(f,g)^(X^-1));


Pos:=Position(EqChnMpsIndex,[IdGroup(S!.group),IdGroup(R!.group),f,X]);
if Pos=fail then
Add(EqChnMpsIndex,[IdGroup(S!.group),IdGroup(R!.group),f,X]);
Add(EqChnMps,F(ModularEquivariantChainMap(S,R,fx)));
Pos:=Length(EqChnMps);
fi;
HKxhomHPKx:=HomologyVectorSpace(EqChnMps[Pos],n);


HKx:=Source(HKxhomHPKx);
HPKx:=Parent(Range(HKxhomHPKx));

HPKxhomHP:=LeftModuleGeneralMappingByImages(HPKx,HP,
#HPKxhomHP:=LeftModuleHomomorphismByImagesNC(HPKx,HP,
Elements(Basis(HPKx)),
Elements(Basis(HP)));

HKxhomHP:=LeftModuleGeneralMappingByImages(HKx,HP,
#HKxhomHP:=LeftModuleHomomorphismByImagesNC(HKx,HP,
Elements(Basis(HKx)),
List(Elements(Basis(HKx)),x->
ImagesRepresentative(HPKhomHP, ImagesRepresentative(HKxhomHPKx,x) )) );

HKhomHKx:=LeftModuleGeneralMappingByImages(HK,HKx,
#HKhomHKx:=LeftModuleHomomorphismByImagesNC(HK,HKx,
Elements(Basis(HK)),
Elements(Basis(HKx))  );

HKhomHP2:=LeftModuleGeneralMappingByImages(HK,HP,
#HKhomHP2:=LeftModuleHomomorphismByImagesNC(HK,HP,
Elements(Basis(HK)),
List(Elements(Basis(HK)),x->
ImagesRepresentative(HKxhomHP, ImagesRepresentative(HKhomHKx,x) )) );

for x in Elements(Basis(HK)) do
Append(HPrels, [ImagesRepresentative(HKhomHP,x)-ImagesRepresentative(HKhomHP2,x)]);
od;

od;
fi;
od;


HPP:=HP/VectorSpace(GF(prime),HPrels);


return Dimension(HPP);;
end;
#####################################################################


#####################################################################
PSeries:=
function(G,p,n)
local
        P,R,L,F;

P:=SylowSubgroup(G,p);
R:=ResolutionPrimePowerGroup(P,n+1);

####################################################################
F:=function(R);
return TensorWithIntegersModP(R,p);
end;
####################################################################

L:=List(Reversed([1..n]),i->PrimePart(G,R,F,i));
Add(L,1);
L:=Reversed(L);

return(PoincareSeries(L,n));

end;
#####################################################################

return PSeries(G,p,n);
end);
#####################################################################
#####################################################################

