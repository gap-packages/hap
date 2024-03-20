RT:=0;
#####################################################################
#####################################################################
InstallGlobalFunction(HAP_CocyclesAndCoboundariesModP,
function(arg)
local  
	C, n,  Dimension,  
	M1, M2, row, 
	dim, dim2, BasisKerd2, BasisImaged1, 
        cohomology,
	CycleToClass, ClassToCycle,
	i, j, x, sum, ONE, char, V, v, lncoho;

C:=arg[1];
n:=arg[2];

char:=EvaluateProperty(C,"characteristic");
ONE:=One(GF(char));
Dimension:=C!.dimension;

if n <0 then return false; fi;
#if n=0 then return [0*ONE]; fi;

if Dimension(n)>0 then
M2:=CoboundaryMatrix(C,n);
else
M2:=[[ONE*0]]; 
fi;

#BasisKerd2:=SemiEchelonMatTransformation(TransposedMat(M2)).relations;
BasisKerd2:=NullspaceMat(TransposedMat(M2));
Unbind(M2);
M1:=CoboundaryMatrix(C,n-1);
if M1=[] then
M1:=[ONE*0*[1..Dimension(n)]];
if Dimension(n)=0 then M1:=[[ONE*0]]; fi;
fi;
#BasisImaged1:=SemiEchelonMatTransformation(TransposedMat(M1)).vectors;
BasisImaged1:=SemiEchelonMat(TransposedMat(M1)).vectors;
Unbind(M1);
dim:=Length(BasisImaged1);

cohomology:=[];
V:=MutableCopyMat(BasisImaged1);
for v in  BasisKerd2 do
Add( V, v*ONE);
V:=MutableCopyMat(SemiEchelonMatTransformation(V).vectors);
dim2:=Length(V);
if dim2>dim then Add(cohomology,v); dim:=dim2; fi;
od;

lncoho:=Length(cohomology);
if Length(cohomology)>0 then
cohomology:=ONE*SemiEchelonMatTransformation(cohomology).vectors;
#else
#cohomology:=[0*ONE];
fi;

cohomology:=Concatenation(cohomology,BasisImaged1);
Unbind(V); Unbind(BasisImaged1);

#####################################################################
CycleToClass:=function(v)
local cls;

cls:= SolutionMat(cohomology,v*ONE);
return cls{[1..lncoho]};
end;
#####################################################################

#####################################################################
ClassToCycle:=function(u)
local w, i;

w:=cohomology[1]*0;
for i in [1..Length(u)] do
if not IsZero(u[i]) then w:=w+u[i]*cohomology[i]; fi;
od;

return w;
end;
#####################################################################

return 	rec(
		cocyclesBasis:=BasisKerd2,
	 	cocycleToClass:=CycleToClass,        
	 	classToCocycle:=ClassToCycle );

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(CupProductOfRegularCWComplexModP,
function(arg)
local X,prime,Xd,XX,A, CX, CXX, XdmapXX, XdmapX, CXmapCXX,
      DIAG, tc,i,n, w, HCXXmapHCX, tensor, quad2pair, HCXcab, 
      HCXsmallcab, CXsmall, HXsmall, CXsmallmapCX, HCXmapHCXsmall,
      HCXsmall,range, CXmapCXsmall, HCXsmallmapHCX, Equiv;

X:=arg[1];
prime:=arg[2];
X:=ContractedComplex(X);
CriticalCells(X);
XX:=DirectProductOfRegularCWComplexes(X,X,Dimension(X)+1);
CXmapCXX:=DiagonalChainMap(X);
HCXXmapHCX:=HomToIntegersModP(CXmapCXX,prime);
Equiv:=ChainComplexEquivalenceOfRegularCWComplex(X);
CXsmallmapCX:=Equiv[2];
HCXmapHCXsmall:=HomToIntegersModP(CXsmallmapCX,prime);
HCXsmall:=Target(HCXmapHCXsmall);   
CXmapCXsmall:=Equiv[1];
HCXsmallmapHCX:=HomToIntegersModP(CXmapCXsmall,prime);

range:=[1..Length(HCXsmall)];
HCXsmallcab:=
List(range,n->HAP_CocyclesAndCoboundariesModP(HCXsmall,n));

quad2pair:=XX!.quad2pair;


##############################
##############################
tensor:=function(p,q,vv,ww)
local  vvv,www,u, P, v, w, Q, i, j, cp, cq, cheat;
#inputs a  p-class vv and a q-class ww.
#outputs a p+q-cocycle u in HCXdsmall
#u:=0*[1..CXX!.dimension(p+q)];
if p=0 then return ww; fi;
if q=0 then return vv; fi;
 
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



