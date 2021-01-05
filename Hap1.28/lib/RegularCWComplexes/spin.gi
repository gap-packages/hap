
#DeclareGlobalFunction("Spin");
####################################################################
################### Spinning CW-complexes ##########################
####################################################################
############ Input: inclusion map from a subcomplex U to ###########
################### its parent (n-dimensional) regular #############
################### CW-complex X ###################################
########### Output: (n+1)-dimensional CW-complex S(X) ##############
################### corresponding to the `spinning' ################
################### of the original complex about U ################
####################################################################
InstallGlobalFunction(Spin, 
function(inc)
local X,Y,U, map, bndX, bndY, bnd, orient, quad2pair, pair2quad,
      indx,bndd,p,pr,INDX,pos,i, j, ij, n, x, y, a, b, ia, ib, 
      count, BND, ORIEN, M;

X:=Target(inc);
U:=Source(inc);
indx:=U!.boundaries;;
map:=inc!.mapping;
Y:=[ [[1,0],[1,0]], [[2,1,2],[2,1,2]], [[2,1,2]], [] ]; 
Y:=RegularCWComplex(Y);
#Y is the unit disk

bndX:=X!.boundaries;
bndY:=Y!.boundaries;
bnd:=List([0..Dimension(X)+Dimension(Y)],i->[]);
orient:=List([0..Dimension(X)+Dimension(Y)],i->[]);

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
##############
od;
fi;

if j>1 then 
b:=bndY[j][y];
for ib in [2..Length(b)] do #b{[2..Length(b)]} do
Add(BND,quad2pair[i][j-1][x][b[ib]][2]);
BND[1]:=BND[1]+1;
##############
##############
od;
fi;
bnd[i-1+j][quad2pair[i][j][x][y][2]]:=BND;
orient[i-1+j][quad2pair[i][j][x][y][2]]:=ORIEN;
od;od;
od;od;

bnd[1]:=List(bnd[1],i->[1,0]);
Add(bnd,[]);
##############################
##############################


#Let's now work on removing those cells of the
#form f x e with e the unique 2-cell in Y and f NOT in U.

indx:=List([1..Length(X!.boundaries)], a-> [1..Length(X!.boundaries[a])]);
for n in [1..Length(indx)-1] do
for j in [1..Length(U!.boundaries[n])] do
indx[n][map(n-1,j)]:=0;
od;
od;

#INDX[n] is a list of n-1-cells       
#bnd[n] is a list of boundaries of n-1-cells
#map(n,i) is the image of the ith n-cell in U
INDX:=List([1..Length(bnd)], a-> [1..Length(bnd[a])]);
for p in [1..1+Dimension(X)] do
for x in [1..X!.nrCells(p-1)] do
if indx[p][x]>0 then
pr:=quad2pair[p][3][x][1];
INDX[pr[1]][pr[2]]:=0;
fi;
od;
od;

for n in [1..Length(INDX)] do
INDX[n]:=Filtered(INDX[n],a->a>0);
od;

for n in [1..Length(bnd)] do
bnd[n]:=bnd[n]{INDX[n]};
od;

###########################
pos:=function(n,i);
return Position(INDX[n],i);
end;
###########################

#For safety, we'll be inefficient and create a copy of bnd rather 
#than modify bnd. This gets around silly mistakes that could arise 
#with pointers.

bndd:=List([1..Length(bnd)],i->[]);

bndd[1]:=bnd[1];
for n in [1..Length(bnd)-1] do
for x in bnd[n+1] do
Add(bndd[n+1], Concatenation([x[1]], List(x{[2..Length(x)]}, i-> pos(n,i))))  ;
od;
od;

M:= RegularCWComplex(1*bndd);
OrientRegularCWComplex(M);

return M;
end);
##################################################
##################################################



