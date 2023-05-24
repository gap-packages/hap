

##################################################
##################################################
InstallGlobalFunction(DirectProductOfRegularCWComplexesLazy, 
function(X,Y)
local bndX, bndY, boundary, 
      quad2pair, pair2quad,orient,
      NrCells,
      i, j, ij, x, y, a, b, ia, ib, count, BND, M, 
      m, n, p, q;

#THIS FUNCTION IS NOT FULLY IMPLEMENTED AND NEEDS TO BE USED WITH CARE

bndX:=X!.boundaries;
bndY:=Y!.boundaries;

###############################
###############################
NrCells:=function(n)
local p,q, dim;

dim:=0;
for p in [0..n] do
q:=n-p;
dim:=dim + X!.nrCells(p)*Y!.nrCells(q);
od;
return dim;  

end;
###############################
###############################

boundary:=[];
boundary[1]:=function(k) return [1,0]; end;
for n in [2..1+Dimension(X)+Dimension(Y)] do
###############################
###############################
boundary[n]:=function(k)
local bnd;
#Boundary of k-th cell in dimension n+1
bnd:=[];

return bnd;
end;
###############################
###############################
od;

###############################
###############################
BND:=List([1..1+Dimension(X)+Dimension(Y)],
          n->[boundary[n], NrCells(n-1)]);
BND:=List(BND,f->LazyList(f[1],[["length",f[2]]]));
###############################
###############################

M:=Objectify(HapRegularCWComplex,
       rec(
           nrCells:=NrCells,
           boundaries:=BND,
           coboundaries:=fail,
           orientation:=fail,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           quad2pair:=fail,
           properties:=[["dimension",Dimension(X)+Dimension(Y)] ]));


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


M!.quad2pair:=quad2pair;
M!.pair2quad:=pair2quad;

return M;
end);
##################################################
##################################################

#Y:=RegularCWComplex(ClosedSurface(1));
#YY:=DirectProductOfRegularCWComplexesLazy(Y,Y);
#Print(Homology(YY,2),"\n");
#Print(YY,"\n");

