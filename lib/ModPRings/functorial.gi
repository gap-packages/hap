
########################################################################
########################################################################
InstallGlobalFunction(HAP_FunctorialModPCohomologyRing,
function(f,prime,X)
local
deg, G,H,GisoG, HisoH, HhommG, RH, RG, RHmapRG, PG, RPG, PH, RPH, AG, AH, APG, APH,
PhomG, PhomH, RPGhomRG, RPHhomRH, mat, Ghom, Hhom, HhomG, BH, BG, BPH, BPG,
op, DimsPG, DimsPH, DimsPHstart, APGElementToVector, VectorToAPHElement,
homims, AlgHom, pos, b,x,dd,n;
#For a group homomorphism H-->G and for each degree n we'll construct a diagram
#
#                     B
#  y=xB   H^n(H)  <-------  H^n(G)    x
#
#           |                 |
#           | A            C  |
#           V                 V
#
#  t=yA  H^n(PH)             h^n(PG)  s=xC
#
#where A,B,C are matrices representing maps x |---> xA, x |---> xB, x |--->xC.
#When y =xA we determine x using x=SolutionMat(A,y).  Here x and y are row 
#vectors
#We write
#         A=Hhom[n]
#         B=HhomG[n]
#         C=Ghom[n]

H:=Source(f);
G:=Target(f);
HhommG:=f;

######CASE WHEN X IS AN INTEGER X=deg
if IsInt(X) then
deg:=X-2;  #THIS IS A BIT SLOPPY
GisoG:=IsomorphismGroups(SmallGroup(IdGroup(G)),G);;
HisoH:=IsomorphismGroups(SmallGroup(IdGroup(H)),H);;

if IsPGroup(H) then 
    RH:=ResolutionPrimePowerGroup(Source(HisoH),deg+1);
else
    RH:=ResolutionFiniteGroup(Source(HisoH),deg+1);;
    pos:=PositionProperty([1..Length(RH!.properties)],i->RH!.properties[i][1]="characteristic");
    RH!.properties[pos][2]:=prime;
fi;
RH!.group:=H;
RH!.elts:=List(RH!.elts,x->Image(HisoH,x));
if IsPGroup(G) then
    RG:=ResolutionPrimePowerGroup(Source(GisoG),deg+1);
else
    RG:=ResolutionFiniteGroup(Source(GisoG),deg+1);
    pos:=PositionProperty([1..Length(RG!.properties)],i->RG!.properties[i][1]="characteristic");
    RG!.properties[pos][2]:=prime;
fi;
RG!.group:=G;
RG!.elts:=List(RG!.elts,x->Image(GisoG,x));
fi;
######CASE WHEN X IS AN INTEGER X=deg DONE

######CASE WHEN X IS A PAIR X[1]=RG, X[2]=RH
if IsList(X) then
RG:=X[1];
RH:=X[2];
deg:=Minimum(Length(RG),Length(RH))-1;
fi;
######CASE WHEN X IS A PAIR X[1]=RG, X[2]=RH DONW


RHmapRG:=EquivariantChainMap(RH,RG,HhommG);; 

PG:=SylowSubgroup(G,prime);;
RPG:=ResolutionPrimePowerGroup(PG,deg+1);;
PH:=SylowSubgroup(H,prime);
RPH:=ResolutionPrimePowerGroup(PH,deg+1);;
AG:=ModPCohomologyRing_alt(G,RPG);;
AH:=ModPCohomologyRing_alt(H,RPH);;
APG:=Parent(AG);
APH:=Parent(AH);

PhomG:=GroupHomomorphismByFunction(PG,G,x->x);;
PhomH:=GroupHomomorphismByFunction(PH,H,x->x);;
RPGhomRG:=EquivariantChainMap(RPG,RG,PhomG);;
RPHhomRH:=EquivariantChainMap(RPH,RH,PhomH);;

mat:=function(x); return TransposedMat(HomomorphismAsMatrix(x)); end;

Ghom:=[];
Hhom:=[];
HhomG:=[];
for n in [1..deg] do
Ghom[n]:=mat(Homology(TensorWithIntegersModP(RPGhomRG,prime),n));;
Hhom[n]:=mat(Homology(TensorWithIntegersModP(RPHhomRH,prime),n));;
HhomG[n]:=mat(Homology(TensorWithIntegersModP(RHmapRG,prime),n));;
od;


BH:=CanonicalBasis(AH);;
BG:=CanonicalBasis(AG);;
BPH:=CanonicalBasis(APH);;
BPG:=CanonicalBasis(APG);;
op:= OperationAlgebraHomomorphism( APG, BPG, OnRight );
DimsPG:=[];
DimsPG[1]:=[[1]];
for n in [1..deg+1] do
DimsPG[n+1]:=Filtered([1..Length(BPG)], i->APG!.degree(BPG[i])=n);
od;

DimsPH:=[];
DimsPHstart:=[];
DimsPHstart[1]:=0;
DimsPH[1]:=[[1]];
for n in [1..deg+1] do
DimsPH[n+1]:=Filtered([1..Length(BPH)], i->APH!.degree(BPH[i])=n);
DimsPHstart[n+1]:=DimsPHstart[n]+Length(DimsPH[n]);
od;


####################################################
APGElementToVector:=function(w)
local v;
#inputs a homogeneous element w in APG.
v:=Image(op,w);
v:=v[1];
v:=v{DimsPG[APG!.degree(w)+1]};
return v;
end;
####################################################

####################################################
VectorToAPHElement:=function(v,d)
local w, i;
#inputs a vector v representing homogeneous element w in APH of degree d.
w:=Zero(APH);
for i in [1..Length(v)] do
w:=w+v[i]*BPH[DimsPHstart[d+1]+i];
od;
return w;
end;
####################################################


homims:=[BH[1]];
for b in Filtered(BG,a->AG!.degree(a)<=deg and AG!.degree(a)>0)   do
dd:=AG!.degree(b);
x:=APGElementToVector(b);
x:=SolutionMat(Ghom[dd]*One(GF(prime)),x);
x:=x*HhomG[dd];
x:=x*Hhom[dd];
x:=VectorToAPHElement(x,dd);
Add(homims,x);
od;

BG:=Filtered(BG,a->AG!.degree(a)<=deg);
AG:=Subalgebra(APG,BG);
AG!.degree:=APG!.degree;
BH:=Filtered(BH,a->AH!.degree(a)<=deg);
AH:=Subalgebra(APH,BH);
AH!.degree:=APH!.degree;

#Use    AlgebraHomomorphismByImages   to check this really is a homomorphism.
AlgHom:=AlgebraHomomorphismByImagesNC(AG,AH,BG,homims);
return AlgHom;
end);
########################################################################
########################################################################

