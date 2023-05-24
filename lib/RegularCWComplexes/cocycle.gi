#(C) Graham Ellis, 2005-2014
#RT:=0;
#####################################################################
#####################################################################
InstallMethod(CoboundaryMatrix,
"coboundary matrix of a cochain complex",
[IsHapCochainComplex,IsInt],
function(C,n)
local A,k;

A:=[];
for k in [1..C!.dimension(n)] do
Add(A, C!.boundary(n,k));
od;

return 1*TransposedMat(A);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallMethod(IdentityMap,
"Identity map on regular CW complex",
[IsHapRegularCWComplex],
function(X)
local map;

map:=function(n,k);
return k;
end;

return Objectify(HapRegularCWMap,
       rec(
           source:=X,
           target:=X,
           mapping:=map));

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ChainMapOfRegularCWMap,
function(f)
local S, T, map, sparsemap;

S:=ChainComplexOfRegularCWComplex(Source(f));
T:=ChainComplexOfRegularCWComplex(Target(f));

################################
map:=function(v,n)
local w,i,x;
w:=0*[1..T!.dimension(n)];
for i in Filtered([1..Length(v)],a-> not IsZero(v[a])) do
x:=f!.mapping(n,i);
if not x=fail then
#w[AbsInt(x)]:=w[AbsInt(x)]+SignInt(x)*v[i];
w[x]:=w[x]+1*v[i];
fi;
od;
return w;
end;
################################

################################
sparsemap:=function(v,n)
local w,i,x;
w:=[];
for i in Filtered([1..Length(v)],a-> not IsZero(v[a])) do
#if not x=fail then
x:=f!.mapping(n,i);
#Add(w,[x,v[i]]);
Add(w,[AbsInt(x),SignInt(x)*v[i]]);
#fi;
od;

return w;
end;
################################

return Objectify(HapChainMap,
       rec(
           source:=S,
           target:=T,
           mapping:=map,
           sparseMap:=sparsemap,
           properties:=[["characteristic", 0],["type","chainMap"]]));


end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexEquivalenceOfRegularCWComplex,
function(X)
local S,T,F,G, mapS, mapT,Deform,Crits,n, Htpy, bij;;

S:=ChainComplexOfRegularCWComplex(X);
T:=ChainComplexOfRegularCWComplexWithVectorField(X,"anything");
Deform:=T!.deform;
Htpy:=T!.htpy;
bij:=T!.bij;

Crits:=[];
for n in [0..Dimension(X)] do
Crits[n+1]:=Filtered(CriticalCells(X),x->x[1]=n);
Apply(Crits[n+1],x->x[2]);
od;

####################
mapT:=function(v,n)
local w, z, i, j, def, pos;
w:=0*[1..T!.dimension(n)];
z:=Filtered([1..Length(v)],i->not IsZero(v[i]));
for i in z do
def:=Deform(n,i);
for j in def do
pos:=Position(Crits[n+1],AbsInt(j));
w[pos]:=w[pos]+SignInt(j)*v[i];
od;
od;
return w;
end;
####################

F:=Objectify(HapChainMap,
       rec(
           source:=S,
           target:=T,
           mapping:=mapT,
           properties:=[["characteristic", 0],["type","chainMap"]]));

####################
mapS:=function(v,n)
local w, i, j, jj, cells,B,SB,Or;
w:=0*[1..S!.dimension(n)];
for i in Filtered( [1..Length(v)], a->not IsZero(v[a])  )  do

cells:=[Crits[n+1][i]];
if n>0 then
B:=X!.boundaries[n+1][Crits[n+1][i]];
B:=B{[2..Length(B)]};
Or:=X!.orientation[n+1][Crits[n+1][i]];
SB:=List([1..Length(B)],k->Or[k]*B[k]);
for j in SB do
Deform(n-1,j);
Append(cells,-SignInt(j)*Htpy[n][AbsInt(j)]);
od;
fi;

for jj in cells do
w[AbsInt(jj)]:=w[AbsInt(jj)]+SignInt(jj)*v[i];
od;


od;
return w;
end;
####################

G:=Objectify(HapChainMap,
       rec(
           source:=T,
           target:=S,
           mapping:=mapS,
           properties:=[["characteristic", 0],["type","chainMap"]]));


return [F,G];
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexHomeomorphismEquivalenceOfRegularCWComplex,
function(X)
local S,T,F,G, mapS, mapT;;
#THIS FUNCTION STILL NEEDS TO BE WRITTEN. THIS IS JUST A PLACEHOLDER!

S:=ChainComplexOfRegularCWComplex(X);
T:=S;

####################
mapT:=function(v,n);
return v;
end;
####################

F:=Objectify(HapChainMap,
       rec(
           source:=S,
           target:=T,
           mapping:=mapT,
           properties:=[["characteristic", 0],["type","chainMap"]]));

####################
mapS:=function(v,n);
return v;
end;
####################

G:=Objectify(HapChainMap,
       rec(
           source:=T,
           target:=S,
           mapping:=mapS,
           properties:=[["characteristic", 0],["type","chainMap"]]));


return [F,G];

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HAP_CocyclesAndCoboundaries,
function(arg)
local  
	C, n, toggle, Dimension,  
	M1, M2, row, 
	dim, BasisKerd2, BasisImaged1, Rels,
	Smith, SmithRecord, TorsionCoefficients,
	ColMat, InvColMat, 
	RemoveRowsMat, InsertRowsList, 
	CycleToClass, ClassToCycle,
	i, j, x, sum;

C:=arg[1];
n:=arg[2];
if Length(arg)>2 then toggle := arg[3]; else toggle := false; fi;

Dimension:=C!.dimension;

if n <0 then return false; fi;
if n=0 then return [0]; fi;

if Dimension(n)>0 then
M2:=CoboundaryMatrix(C,n);
else
M2:=[0*[1..Dimension(n-1)]]; 
fi;

BasisKerd2:=LLLReducedBasis(TransposedMat(M2),"linearcomb").relations;
Unbind(M2);
M1:=CoboundaryMatrix(C,n-1);
BasisImaged1:=LLLReducedBasis(TransposedMat(M1)).basis;
Unbind(M1);
dim:=Length(BasisImaged1);

Rels:=[];
for i in [1..dim] do
        Rels[i]:=SolutionMat(BasisKerd2,BasisImaged1[i]);
od;

if Length(Rels)=0 and Length(BasisKerd2)>0 then
Append(Rels,[List([1..Length(BasisKerd2)],x->0)]);
fi;				#CHECK THE MATHS HERE!

if toggle=false then 
return rec(
		cocyclesBasis:=BasisKerd2, 
		boundariesCoefficients:=Rels,
		torsionCoefficients:=fail,
		cocycleToClass:=fail,
		classToCocycle:=fail );
fi;
	################STOP HERE IF TOGGLE=FALSE####################


SmithRecord:= SmithNormalFormIntegerMatTransforms(Rels);
Smith:=SmithRecord.normal;
ColMat:=TransposedMat(SmithRecord.coltrans);
InvColMat:=Inverse(ColMat);	#Only valid for finite groups

TorsionCoefficients:=[];
for i in [1..Length(BasisKerd2)] do
if i<=Length(Smith) then
    TorsionCoefficients[i]:=Smith[i][i];
else
    TorsionCoefficients[i]:=0;
fi;
od;

InsertRowsList:=[];
RemoveRowsMat:=IdentityMat(Length(TorsionCoefficients));
for i in [1..Length(BasisKerd2)] do
	if TorsionCoefficients[i]=1 then 
	RemoveRowsMat[i]:=47;
	Append(InsertRowsList,[i]);
	fi;
od;
RemoveRowsMat:=Filtered(RemoveRowsMat,r->not (r=47));
if Length(RemoveRowsMat)=0 then
TorsionCoefficients:=[];
else
TorsionCoefficients:= RemoveRowsMat*TorsionCoefficients;
fi;

#####################################################################
CycleToClass:=function(v)
local u, i;

u:=SolutionMat(BasisKerd2,v);
u:=ColMat*u;
u:=RemoveRowsMat*u;

for i in [1..Length(u)] do
if TorsionCoefficients[i]>0 then
u[i]:=u[i]  mod TorsionCoefficients[i];
fi;
od;

return u;
end;
#####################################################################

#####################################################################
ClassToCycle:=function(u)
local v,w, i, temp;

for i in [1..Length(u)] do
if TorsionCoefficients[i]>0 then
u[i]:=u[i] mod TorsionCoefficients[i];     
fi;
od;

v :=[];
temp:=0;
for i in [1..Length(BasisKerd2)] do
if i in InsertRowsList then v[i]:=0; 
else 
temp:=temp+1;
v[i] := u[ temp ];
fi;
od;

v:=InvColMat*v;

w:=[];
for i in [1..Dimension(n)] do
w[i]:=0;
od;

for i in [1..Length(v)] do
w:=w + v[i]*BasisKerd2[i];
od;

return w;
end;
#####################################################################

return 	rec(
		cocyclesBasis:=BasisKerd2,
	 	boundariesCoefficients:=Rels,
	 	torsionCoefficients:=TorsionCoefficients,
	 	cocycleToClass:=CycleToClass,        
	 	classToCocycle:=ClassToCycle );

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HAP_IntegralCohomology, 
function(C,n)
local A, i, Smith, TorsionCoefficients;

A:=HAP_CocyclesAndCoboundaries(C,n);
Smith:= SmithNormalFormIntegerMat(A.boundariesCoefficients);

TorsionCoefficients:=[];

for i in [1..Length(A.cocyclesBasis)] do
        if i<=Length(Smith) then
        TorsionCoefficients[i]:=Smith[i][i];
        else
        TorsionCoefficients[i]:=0;
        fi;
od;

return Filtered(TorsionCoefficients, i-> not (i=1));

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HAP_IntegralClassToCocycle,  
function(arg)
local
	C, u, n, A,
	i;

C:=arg[1];
u:=arg[2];
n:=arg[3];
if Length(arg)>3 then 
A:=arg[4];
else
A:=HAP_CocyclesAndCoboundaries(C,n,true);
fi;

return A.classToCocycle(u);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(HAP_IntegralCocycleToClass, 
function(arg)
local
        C, v, n, A,
        i;

C:=arg[1];
v:=arg[2];
n:=arg[3];
if Length(arg)>3 then
A:=arg[4];
else
A:=HAP_CocyclesAndCoboundaries(C,n,true);
fi;

return A.cocycleToClass(v);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(CupProductOfRegularCWComplex_alt,
function(arg)
local X,Xd,XX,A, CX, HCX, CXd, CXX, XdmapXX, XdmapX, CXdmapCXX, 
      DIAG, tc,i,n, w, HCXXmapHCXd, tensor, quad2pair, HCXcab, HCXd, 
      HCXdsmallcab, CXdsmall, HXdsmall, CXdsmallmapCXd, HCXdmapHCXdsmall,
      HCXdsmall,range; 

X:=arg[1];
X:=ContractedComplex(X);
DIAG:=DiagonalApproximation(X);;
XdmapXX:=DIAG.inclusion;
XX:=Target(XdmapXX);
Xd:=Source(XdmapXX);
CXdmapCXX:=ChainMapOfRegularCWMap(XdmapXX);
HCXXmapHCXd:=HomToIntegers(CXdmapCXX);
CXd:=Source(CXdmapCXX);
CXX:=Target(CXdmapCXX);
XdmapX:=DIAG.projection;
CX:=ChainComplexOfRegularCWComplex(Target(XdmapX));
HCX:=HomToIntegers(CX);
#HCXd:=HomToIntegers(CXd);
CXdsmallmapCXd:=ChainComplexEquivalenceOfRegularCWComplex(Xd)[2];
CXdsmall:=Source(CXdsmallmapCXd);
HCXdmapHCXdsmall:=HomToIntegers(CXdsmallmapCXd);
HCXdsmall:=Target(HCXdmapHCXdsmall);
HCXd:=Source(HCXdmapHCXdsmall);

range:=[1..Length(CX)];
if Length(arg)=2 then range:=arg[2]; fi;

HCXcab:=
List(range,n->HAP_CocyclesAndCoboundaries(HCX,n,true));

HCXdsmallcab:=
List(range,n->HAP_CocyclesAndCoboundaries(HCXdsmall,n,true));

quad2pair:=XX!.quad2pair;

##############################
##############################
tensor:=function(p,q,vv,ww)
local  u, P, v, w, Q, i, j, cp, cq, cheat;
#inputs a  p-class vv and a q-class ww.
#outputs a p+q-cocycle u in HCXdsmall
#u:=0*[1..CXX!.dimension(p+q)];
u:=[];
v:=HCXcab[p]!.classToCocycle(vv);
w:=HCXcab[q]!.classToCocycle(ww);
P:=Filtered([1..Length(v)],i->not IsZero(v[i]));
Q:=Filtered([1..Length(w)],i->not IsZero(w[i]));

#cp:=Filtered(CriticalCells(X), x-> x[1]=p);
#cp:=List(cp,y->y[2]);
#Print(P,"   ", cp);
#P:=Filtered(P, i-> i in cp);
#cq:=Filtered(CriticalCells(X), x-> x[1]=q);
#cq:=List(cq,y->y[2]);
#Q:=Filtered(Q, i-> i in cq);
#For the above to work we need to choose the vector field on Xdsmall
#compatible with that on X. I'll do this when I write a lazy implementation
#of the direct product XxX.


for i in P do
for j in Q do
#u[quad2pair[p+1][q+1][i][j][2]] := v[i]*w[j];
Add(u,  [ quad2pair[p+1][q+1][i][j][2] , v[i]*w[j] ]);
od;od;


###############################
if p<q then cheat:=1;         #NEED TO THINK ABOUT THIS!
else cheat:=-1;
fi;
Apply(u,x->[x[1],cheat*x[2]]);
###############################

w:= HCXXmapHCXd!.sparseMap(u, p+q);
w:=HCXdmapHCXdsmall!.mapping(w,p+q);
return HCXdsmallcab[p+q]!.cocycleToClass(w);
end;
##############################
##############################

return tensor;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(CupProductOfRegularCWComplex,
function(arg)
local X,Xd,XX,A, CX, CXX, XdmapXX, XdmapX, CXmapCXX,
      DIAG, tc,i,n, w, HCXXmapHCX, tensor, quad2pair, HCXcab, 
      HCXsmallcab, CXsmall, HXsmall, CXsmallmapCX, HCXmapHCXsmall,
      HCXsmall,range, CXmapCXsmall, HCXsmallmapHCX, Equiv;

X:=arg[1];
X:=ContractedComplex(X);
CriticalCells(X);
XX:=DirectProductOfRegularCWComplexes(X,X,Dimension(X)+1);
CXmapCXX:=DiagonalChainMap(X);
HCXXmapHCX:=HomToIntegers(CXmapCXX);
Equiv:=ChainComplexEquivalenceOfRegularCWComplex(X);
CXsmallmapCX:=Equiv[2];
HCXmapHCXsmall:=HomToIntegers(CXsmallmapCX);
HCXsmall:=Target(HCXmapHCXsmall);   
CXmapCXsmall:=Equiv[1];
HCXsmallmapHCX:=HomToIntegers(CXmapCXsmall);

range:=[1..Length(HCXsmall)];
HCXsmallcab:=
List(range,n->HAP_CocyclesAndCoboundaries(HCXsmall,n,true));

quad2pair:=XX!.quad2pair;

##############################
##############################
tensor:=function(p,q,vv,ww)
local  vvv,www,u, P, v, w, Q, i, j, cp, cq, cheat;
#inputs a  p-class vv and a q-class ww.
#outputs a p+q-cocycle u in HCXdsmall
#u:=0*[1..CXX!.dimension(p+q)];
u:=[];
vvv:=HCXsmallcab[p]!.classToCocycle(vv);
www:=HCXsmallcab[q]!.classToCocycle(ww);
v:=HCXsmallmapHCX!.mapping(vvv,p);
w:=HCXsmallmapHCX!.mapping(www,q);

P:=Filtered([1..Length(v)],i->not IsZero(v[i]));
Q:=Filtered([1..Length(w)],j->not IsZero(w[j]));

for i in P do
for j in Q do
Add(u,  [ quad2pair[p+1][q+1][i][j][2] , v[i]*w[j] ]);
od;od;


###############################
if p<q then cheat:=1;         #NEED TO THINK ABOUT THIS!
else cheat:=-1;
fi;
Apply(u,x->[x[1],cheat*x[2]]);
###############################

#RT:=RT-Runtime();  #Takes all the time!
w:= HCXXmapHCX!.sparseMap(u, p+q);
#RT:=RT+Runtime();
w:=HCXmapHCXsmall!.mapping(w,p+q);
w:=HCXsmallcab[p+q]!.cocycleToClass(w);
return w;
end;
##############################
##############################

return tensor;
end);
#####################################################################
#####################################################################



################################################
################################################
InstallGlobalFunction(HAP_CupProductOfPresentation,
function(G)
local R, C, cup, one,two, diag, tensor , rel, sn,i,j,k,x,t,ii;

R:=ResolutionAsphericalPresentation(G,3);
#C:=TensorWithIntegers(R);
C:=HomToIntegers(R);

one:=HAP_CocyclesAndCoboundaries(C,1,true);
two:=HAP_CocyclesAndCoboundaries(C,2,true);

diag:=[];

for k in [1..R!.dimension(2)] do
rel:=List(R!.boundary(2,k),x->x[1]);
tensor:=[];

for i in [1..Length(rel)]  do
ii:=i+1;
if rel[i]<0 then ii:=i; fi;
for j in [ii..Length(rel)]  do
#if AbsInt(rel[i]) <= AbsInt(rel[j]) then sn:=1; fi;
#if AbsInt(rel[i]) > AbsInt(rel[j]) then sn:=-1; fi;
#sn:=sn*SignInt(rel[i])*SignInt(rel[j]);
sn:=SignInt(rel[i])*SignInt(rel[j]);  #Changed May 2017
Add(tensor,[sn*AbsInt(rel[i]), AbsInt(rel[j])] );
od;
od;

Add(diag,tensor);
od;

############################
cup:=function(a,b)
local A,B, AB, r,t;
A:=one.classToCocycle(a);
B:=one.classToCocycle(b);
AB:=0*[1..C!.dimension(2)];

for r in [1..Length(diag)] do
for t in diag[r] do
AB[r]:=AB[r]+SignInt(t[1])*A[AbsInt(t[1])]*B[t[2]];

od;
od;

return two.cocycleToClass(AB);
end;
###########################

###########################
if Cohomology(C,2)=[] then
cup:=function(a,b); return [0]; end;
fi;
###########################

return cup;
end);
################################################
################################################

