
###########################################
###########################################
InstallGlobalFunction(QuotientChainMap,
function(X,ddata)
local data,CX, CT, B, Babs, SSBabs, DimensionT, BoundaryT, fnc, dimsT, ChnMap, mapping, n, v, x, tmp;

data:=1*ddata;

CX:=ChainComplexOfRegularCWComplex(X);

####################################
####################################
fnc:=function(rels,n)
local A, B,  x, i, j, k, ln;
# inputs the relations rels identifying cells of dimension n
# and outputs the list B where cell e_i^n gets sent to +/-e^n_{B[i]}
# by the quotient map. The sign is also recorded.
ln:=X!.nrCells(n);

A:=List([1..ln],i->[i]);

for x in rels do
i:=PositionProperty(A,a->x[1] in a or -x[1] in a);
j:=PositionProperty(A,a->x[2] in a or -x[2] in a);

if not i=j then
if (x[1] in A[i] and x[2] in A[j])
or
(-x[1] in A[i] and -x[2] in A[j])
then
Add(A, Concatenation(A[i],A[j]));
else
Add(A, Concatenation(A[i],-A[j]));
fi;
Remove(A,Minimum(i,j));
Remove(A,Maximum(i,j)-1);
fi;
od;

B:=[];
for i in [1..ln] do
j:=PositionProperty(A,a->i in a or -i in a);

k:=Minimum(List(A[j],AbsInt));
if (i in A[j] and k in A[j])
or
(not i in A[j] and not k in A[j])
then
B[i]:=k;
else B[i]:=-k;
fi;
od;
return B;
end;
####################################
####################################

B:=List([0..Dimension(X)],k->fnc(data[k+1],k));
Babs:=List(B,x->List(x,AbsInt));
SSBabs:=List(Babs,x->SSortedList(x));
#dimsT:=List(B,x->Length(SSortedList(List(x,AbsInt))) );   
dimsT:=List(SSBabs,x->Length(x));

#############################
DimensionT:=function(n);
if n<0 or n>Dimension(X) then return 0; fi;
return dimsT[n+1]; 
end;
#############################

#############################
mapping:=function(n,k)
local v,a;
#function: cells of X ----> cells of T   
#which outputs the signed number k' of the n-dimensional cell
#in the quotient corresponding to the k-th n-dimensional cell 
v:=0*[1..DimensionT(n)];
a:=B[n+1][k];
v[AbsInt(a)]:=SignInt(a);
return v;
end;
#############################

#############################
BoundaryT:=function(n,kk)
local m, bnd, bnd1,i,j,jj,k;

k:=SSBabs[n+1][kk];
m:=Position(Babs[n+1],k);
if k in B[n+1] then bnd:= 1*CX!.boundary(n,m);;
else
bnd:= -1*CX!.boundary(n,m);;
fi;

bnd1:=0*[1..DimensionT(n-1)];

for i in [1..CX!.dimension(n-1)] do
j:=B[n][i];
jj:=Position(SSBabs[n],AbsInt(j));
bnd1[AbsInt(jj)]:=bnd1[AbsInt(jj)]-SignInt(j)*bnd[i];
od;

return bnd1;
end;
#############################

CT:=Objectify(HapChainComplex,
           rec(
           dimension:=DimensionT,
           boundary:=BoundaryT,
           properties:=[
           ["length",Length(CX)],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

ChnMap:=Objectify(HapChainMap,
       rec(
           source:=CX,
           target:=CT,
           mapping:=mapping,
           properties:=[["characteristic", 0],["type","chainMap"]]));
ChnMap!.B:=B;
return ChnMap;
end);
###########################################
###########################################


