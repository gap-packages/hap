


##################################################
##################################################
InstallGlobalFunction(DirectProductOfRegularCWComplexes, 
function(X,Y)
local bndX, bndY, bnd, orientX, orientY, orient,
      quad2pair, pair2quad,
      XYmapX, XYmapY, XYmappingX, XYmappingY,
      i, j, ij, x, y, a, b, ia, ib, count, BND, ORIEN, M, bool;

bndX:=X!.boundaries;
bndY:=Y!.boundaries;

bnd:=List([0..Dimension(X)+Dimension(Y)],i->[]);

bool:=IsBound(Y!.orientation) and IsBound(X!.orientation);
if bool then
orientX:=X!.orientation;
orientY:=Y!.orientation;
orient:=List([0..Dimension(X)+Dimension(Y)],i->[]);
fi;

###############################
###############################
quad2pair:=[];
pair2quad:=List([1..1+Dimension(X)+Dimension(Y)],i->[]);;
count:=List([1..1+Dimension(X)+Dimension(Y)],i->0);

for i in [1..1+Dimension(X)] do
quad2pair[i]:=[];
for j in [1..1+Dimension(Y)] do
quad2pair[i][j]:=[];
ij:=i-1+j;
for x in [1..Length(bndX[i])] do
quad2pair[i][j][x]:=[];
for y in [1..Length(bndY[j])] do
count[ij]:=count[ij]+1;
quad2pair[i][j][x][y]:=[ij,count[ij]];
pair2quad[ij][count[ij]]:=[i,j,x,y];
od; od;
od; od;
###############################
###############################

###############################
###############################
for i in [1..1+Dimension(X)] do
for j in [1..1+Dimension(Y)] do
for x in [1..Length(bndX[i])] do
for y in [1..Length(bndY[j])] do
BND:=[0];
ORIEN:=[];

if i>1 then
a:=bndX[i][x];

for ia in [2..Length(a)] do   #a{[2..Length(a)]} do
Add(BND,quad2pair[i-1][j][a[ia]][y][2]);
BND[1]:=BND[1]+1;
##############
if bool then
Add(ORIEN,(-1)^(i)*orientX[i][x][ia-1]);
fi;
##############
od;
fi;

if j>1 then 
b:=bndY[j][y];
for ib in [2..Length(b)] do #b{[2..Length(b)]} do
Add(BND,quad2pair[i][j-1][x][b[ib]][2]);
BND[1]:=BND[1]+1;
##############
if bool then
Add(ORIEN,(-1)^(i)*orientY[j][y][ib-1]);
fi;
##############
od;
fi;
#Add(bnd[i-1+j],BND);
#Add(orient[i-1+j],ORIEN);
bnd[i-1+j][quad2pair[i][j][x][y][2]]:=BND;
orient[i-1+j][quad2pair[i][j][x][y][2]]:=ORIEN;
od;od;
od;od;

bnd[1]:=List(bnd[1],i->[1,0]);
Add(bnd,[]);
##############################
##############################

if bool then
M:= RegularCWComplex(bnd,orient);
else
M:= RegularCWComplex(bnd);
OrientRegularCWComplex(M);
fi;

#######################################################
#### NEED TO THINK ABOUT PROJECTION MAPS BEFORE COMPLETING
#### THIS SECTION!

########################
XYmappingX:=function(n,k)
local pq;
pq:=pair2quad[n+1][k];
if pq[2]=1 then  return (-1)^(pq[1])*pq[3]; #WRONG! MUST CORRECT
else return fail; fi;
end;
########################

########################
XYmappingY:=function(n,k)
local pq;
pq:=pair2quad[n+1][k];
if pq[1]=1 then return (-1)^pq[2]*pq[4];
else return fail; fi;
end;
########################

####
#### END OF INCOMPLETE SECTION
######################################################

XYmapX:=Objectify(HapRegularCWMap,
       rec(
           source:=M,
           target:=X,
           mapping:=XYmappingX));

XYmapY:=Objectify(HapRegularCWMap,
       rec(
           source:=M,
           target:=Y,
           mapping:=XYmappingY));



M!.firstProjection:=XYmapX;
M!.secondProjection:=XYmapY;
M!.quad2pair:=quad2pair;

return M;
end);
##################################################
##################################################

##################################################
##################################################
InstallGlobalFunction(DiagonalApproximation,
function(X)
local W, M, R, RM, D, DD, DDinv, bound, orient, fp,
      mapW, mapX, maprec, x, y, i, j, n, quad2pair, b, MmapW, MmapX; 

OrientRegularCWComplex(X);

W:=DirectProductOfRegularCWComplexes(X,X);
quad2pair:=W!.quad2pair;

D:=[];
DD:=List(W!.boundaries,i->[]);;
DDinv:=List(W!.boundaries,i->[]);;
maprec:=List(W!.boundaries,i->[]);;

for n in [1..Length(X!.boundaries)] do
for i in [1..Length(X!.boundaries[n])] do
Add(D,quad2pair[n][n][i][i]);
od;
od;

for x in D do
Add(DD[x[1]],x[2]);
b:=BoundaryOfRegularCWCell(W,x[1]-1,x[2]);
for y in b do
Add(DD[y[1]],y[2]);
od;
od;

Apply(DD,x->SSortedList(x));
bound:=[];
orient:=[];
for n in [1..Length(DD)] do
bound[n]:=List(DD[n],x->StructuralCopy(W!.boundaries[n][x]));
orient[n]:=List(DD[n],x->StructuralCopy(W!.orientation[n][x]));
for i in [1..Length(DD[n])] do
DDinv[n][DD[n][i]]:=i;
maprec[n][i]:=DD[n][i];
od;
od;

for n in [2..Length(bound)] do
for x in bound[n] do
for j in [2..Length(x)] do
x[j]:=DDinv[n-1][x[j]];
od;
od;
od;

M:= RegularCWComplex(bound);
M!.orientation:=orient;

R:=DeformationRetract(M);
RM:=Source(R);

#########################
mapW:=function(n,i)
local ii;
ii:=R!.mapping(n,i);
return maprec[n+1][ii];
end;
#########################


MmapW:= Objectify(HapRegularCWMap,
       rec(
           source:=RM,
           target:=W,
           mapping:=mapW));
#fp:=W!.firstProjection;
fp:=W!.secondProjection;
fp:=fp!.mapping;

#########################
mapX:=function(n,i);
return fp(n,mapW(n,i));
end;
#########################

MmapX:= Objectify(HapRegularCWMap,
       rec(
           source:=RM,
           target:=X,
           mapping:=mapX));


return rec(projection:=MmapX, inclusion:=MmapW);



end);
##################################################
##################################################

##################################################
##################################################
InstallGlobalFunction(CWMap2ChainMap,
function(F)
local C,D, X,Y, map, critX, critY, L, n, i;

if not IsHapRegularCWMap(F) then 
Print("Input must be a map of regular CW-complexes.\n");
fi;

X:=Source(F);
Y:=Target(F);
L:=CriticalCellsOfRegularCWComplex(X);
critX:=[];
for n in [0..Dimension(X)] do
critX[n+1]:=Filtered(L,x->x[1]=n); #I hope the order is preserved!!
Apply(critX[n+1],x->x[2]);
od;
L:=CriticalCellsOfRegularCWComplex(Y);
critY:=[];
for n in [0..Dimension(Y)] do
critY[n+1]:=Filtered(L,x->x[1]=n); #I hope the order is preserved!!
Apply(critY[n+1],x->x[2]);
od;
C:=ChainComplexOfRegularCWComplexWithVectorField(X,"anything");
#for n in [0..Length(C)] do
#for i in [1..X!.nrCells(n)] do
#C!.deform(n,i);
#od;od;
D:=ChainComplex(Y);


####################
####################
map:=function(v,n)
local w,i,j,jj,k,kk, x, cells,B,Or;

w:=List([1..D!.dimension(n)],i->0);
#for i in [1..Length(v)] do
for i in Filtered( [1..Length(v)], i->not IsZero(v[i])  )  do

cells:=[critX[n+1][i]];
if n>0 then
B:=X!.boundaries[n+1][critX[n+1][i]];
Or:=X!.orientation[n+1][critX[n+1][i]];;
B:=List([1..Length(Or)],s->Or[s]*B[s+1]);
#for j in X!.boundaries[n+1][critX[n+1][i]] do
for j in B do
C!.deform(n-1,j);
Append(cells,-SignInt(j)*C!.htpy[n][AbsInt(j)]);
od;
fi;

for jj in cells do
j:=F!.mapping(n,AbsInt(jj));
if not j=fail then
j:=SignInt(jj)*j;
x:=D!.deform(n,j);
  for k in x do
  kk:=AbsInt(k);
  kk:=Position(critY[n+1],kk); #Could easily speed this line up!!
  w[kk]:=w[kk]+SignInt(k)*v[i];
  od;
fi;
od;
od;
return w;
end;
####################
####################

return Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=map,
           properties:=[ ["type","chainMap"],
           ["characteristic", 0]
           ]));

end);
##################################################
##################################################



