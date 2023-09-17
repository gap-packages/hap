
#################################################################
#################################################################
InstallGlobalFunction(FilteredCubicalComplexToFilteredRegularCWComplex,
function(arg)
local M, C, BinLst, LstBin, dim, dimm, flen,
      Boundaries, Coboundaries, filtration, Orientation, Properties,
      Boundary, ArrayValueDim, n, i, j, v, bnd, k, b;

M:=arg[1];
if Length(arg)>1 then dim:=arg[2];
else
dim:=EvaluateProperty(M,"dimension");
fi;

flen:=Maximum(Flat(M!.filtration));

C:=SparseFilteredChainComplexOfFilteredCubicalComplex(M);
BinLst:=C!.coordinateToPosition;
LstBin:=C!.positionToCoordinate;

dimm:=Position(List([0..dim],i->C!.dimension(i)),0);
if dimm=fail then dimm:=dim; else dimm:=dimm-2; fi;
Properties:=[["dimension", dimm ],["filtration_length",flen]];


ArrayValueDim:=ArrayValueFunctions(EvaluateProperty(M,"dimension"));


#######################################
Boundary:=function(n,j)
local pos, neg, bnd;
bnd:=C!.boundary(n,j);
pos:=Filtered(bnd,x->x[2]>0);
Apply(pos,x->x[1]);
neg:=Filtered(bnd,x->x[2]<0);
Apply(neg,x->x[1]);

return [pos,neg];
end;
################################


##############################
Boundaries:=[];
Boundaries[1]:=List([1..C!.dimension(0)],x->[1,0]);
Orientation:=[];
Orientation[1]:=List([1..C!.dimension(0)],x->[1]);

for n in [1..dim] do
Boundaries[n+1]:=[];
Orientation[n+1]:=[];
for i in [1..C!.dimension(n)] do
v:=StructuralCopy(Boundary(n,i));
bnd:=Concatenation([Length(v[1])+Length(v[2])],Flat(v));
Add(Boundaries[n+1],bnd);
Add(Orientation[n+1],
Concatenation(List([1..Length(v[1])],a->1),List([1..Length(v[2])],a->-1)));
od;

od;
Boundaries[dim+2]:=[];
##############################

### COBOUNDARIES BEGIN ######################
Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
for n in [0..dim] do
  k:=2*(n+1)+1;#k:=1+2^(n+1);
  Coboundaries[n+1]:=List(Boundaries[n+1],i->[0]);
  for j in [1..Length(Boundaries[n+2])] do
    b:=Boundaries[n+2][j];
    #k:=Length(b);
    for i in b{[2..k]} do
      Coboundaries[n+1][i][1]:=Coboundaries[n+1][i][1]+1;
      Add(Coboundaries[n+1][i],j);
    od;
  od;
#  for b in Coboundaries[n+1] do
#  Append(b,List([1..Length(b)-1],a->1));
#  od;
od;
Coboundaries[dim+1]:=List(Boundaries[dim+1],a->[0]);
### COBOUNDARIES END ###############################



return Objectify(HapFilteredRegularCWComplex,
       rec(
           nrCells:=C!.dimension,
           filteredDimension:=C!.filteredDimension,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=Orientation,
           properties:=Properties));

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(StructuralCopyOfFilteredRegularCWComplex,
function(Y)
local nrcells, Boundaries, filtdim, filtdimrec, F, n, i;

Boundaries:=1*Y!.boundaries;
F:=EvaluateProperty(Y,"filtration_length");

nrcells:=function(n);
return Length(Boundaries[n+1]);
end;

filtdimrec:=[];
for n in [0..Dimension(Y)] do
filtdimrec[n+1]:=[];
for i in [1..F] do
filtdimrec[n+1][i]:=Y!.filteredDimension(i,n);
od;
od;

filtdim:=function(i,n);
return filtdimrec[n+1][i];
end;

return Objectify(HapFilteredRegularCWComplex,
       rec(
           nrCells:=nrcells,
           filteredDimension:=filtdim,
           boundaries:=Boundaries,
           coboundaries:=1*Y!.coboundaries,
           vectorField:=StructuralCopy(Y!.vectorField),
           inverseVectorField:=StructuralCopy(Y!.inverseVectorField),
           criticalCells:=StructuralCopy(Y!.criticalCells),
           orientation:=StructuralCopy(Y!.orientation),
           properties:=StructuralCopy(Y!.properties)));

end);
#################################################################
#################################################################


##########################################################
##########################################################
InstallGlobalFunction(SparseChainComplexOfFilteredRegularCWComplex,
function(Y)
local
        C;

C:=SparseChainComplexOfRegularCWComplex(Y);


return
Objectify(HapFilteredSparseChainComplex,
           rec(
           dimension:=C!.dimension,
           boundary:=C!.boundary,
           filteredDimension:=Y!.filteredDimension,
           properties:=[
           ["length",EvaluateProperty(Y,"dimension")],
           ["filtration_length",EvaluateProperty(Y,"filtration_length")],
           ["type","FilteredChainComplex"],
           ["characteristic",0]]
           ));


end);
##########################################################
##########################################################

#############################################
#############################################
InstallOtherMethod(Dimension,
"Dimension of  filtered regular CW-complex",
[IsHapFilteredRegularCWComplex],
function(f) return EvaluateProperty(f,"dimension");
return EvaluateProperty(f,"dimension");
end);
#############################################
#############################################

#############################################
#############################################
InstallOtherMethod(Dimension,
"Dimension of  filtered simplicial complex",
[IsHapFilteredSimplicialComplex],
function(f) return EvaluateProperty(f,"dimension");
return EvaluateProperty(f,"dimension");
end);
#############################################
#############################################


#############################################
#############################################
InstallOtherMethod(Size,
"Size of  filtered regular CW-complex",
[IsHapFilteredRegularCWComplex],
function(Y)
local sz, i; 
sz:=0;
for i in [0..Dimension(Y)] do
sz:=sz+Y!.nrCells(i);
od;
return sz;
end);
#############################################
#############################################


#############################################
#############################################
InstallGlobalFunction(HAPContractFilteredRegularCWComplex,
function(arg)
local Y;
#TO BE WRITTEN
end);
############################################
############################################

############################################
############################################
InstallGlobalFunction(ContractedFilteredRegularCWComplex,
function(Y)
local  Nbounds, Ncobounds, fl, dim, n, t, k, FiltDim, T, S, x, y, s, bnd, cob,
       pos, bool, perm, cnt, d, d1, F, i, W, dims, allowables, retained, bas,
       Tboundaries,Tcoboundaries,red;

Nbounds:=1*Y!.boundaries;
Ncobounds:=1*Y!.coboundaries;
fl:=EvaluateProperty(Y,"filtration_length");
dim:=Dimension(Y);

allowables:=[];         #allowables[t+1][n+1] is the list of the n-cell numbers 
for t in [1..fl+2] do   #removed during collapse of the t-term of filtration.
allowables[t]:=List([0..dim],n->[]);
od;
#####################################
FiltDim:=function(t,d);
if t<1 then return 0;
else return Y!.filteredDimension(t,d);
fi;
end;
#####################################

retained:=[];;  #retained[t+1][n+1] is the list of the n-cell numbers remaining
for t in [1..fl+2] do  #in t-th term of the filtration after collapse
retained[t]:=List([0..dim],i->[]);
od;

for t in [1..fl] do
 T:=FiltrationTerm(Y,t);
 Tboundaries:=T!.boundaries;
 Tcoboundaries:=T!.coboundaries;
 bool:=true;
 while bool do
  bool:=false;
  for n in Reversed([0..dim-1]) do
   for k in Concatenation(allowables[t][n+1],[FiltDim(t-1,n)+1..FiltDim(t,n)]) do
    if Tcoboundaries[n+1][k][1]=1 then #coboundary of cell k in dimension n 
     x:=Tcoboundaries[n+1][k][2];       #consists only of cell x in dimension n+1.
bool:=true;
     bnd:=Tboundaries[n+2][x];   #if cell k not in t-1 st term then all cells in 
     bnd:=bnd{[2..Length(bnd)]}; #its coboundary also not in t-1 st term
     for y in bnd do
      cob:=Tcoboundaries[n+1][y];
      s:=cob[1]-1;
      cob:=cob{[2..Length(cob)]};
      pos:=Position(cob,x);
      Remove(cob,pos);
      Tcoboundaries[n+1][y]:=Concatenation([s],cob);
     od;
     Add( allowables[t+1][n+1] , k);
     Add( allowables[t+1][n+2] , x);
     Tcoboundaries[n+1][k][1]:=0;
    fi;
   od;
  od;
 od;
od;

retained:=[List([0..dim+1],i->[])];;
for t in [1..fl] do
retained[t+1]:=[];
for n in [0..dim] do
retained[t+1][n+1]:=Difference([1..FiltDim(t,n)],allowables[t+1][n+1]);
od;
od;

perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

red:=List([1..dim+1],i->[]);

for d in [1..dim+1] do
   perm[d]:=[];
   cnt:=0;
      for t in [1..fl] do
         bas:=Difference(retained[t+1][d],retained[t][d]);
Append(red[d],bas);
         for i in [1..Length(bas)] do
            perm[d][bas[i]]:=i+cnt;
         od;
         cnt:=cnt+Length(bas);
      od;
od;

for d in [1..dim+1] do
Nbounds[d]:=Nbounds[d]{red[d]};
od;

dims:=[];
for t in [1..fl] do
dims[t]:=[];
for d in [0..dim] do
dims[t][d+1]:= Length(retained[t+1][d+1]);
od;
od;


for d in [1..dim+1] do
if d>1 then
d1:=d-1;
for x in Nbounds[d] do
for i in [2..Length(x)] do
x[i]:=perm[d1][x[i]];
od;
od;
fi;
od;

W:=RegularCWComplex(Nbounds);;      

############################################
FiltDim:=function(t,d)
local m;
if t<1 then return 0; fi;

return dims[t][d+1];
end;
############################################



return Objectify(HapFilteredRegularCWComplex,
rec(
           nrCells:=W!.nrCells,
           filteredDimension:=FiltDim,
           boundaries:=W!.boundaries,
           coboundaries:=W!.coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=W!.orientation,
           properties:=[ ["dimension",Dimension(W)], ["filtration_length",fl]]));


end);
############################################
############################################


#####################################################################
#####################################################################
FilteredSimplicialComplexToFilteredCWComplex:=
function(M)
local
         Y, fd;

Y:=SimplicialComplexToRegularCWComplex(M);
Add(Y!.properties,["filtration_length",M!.filtrationLength]);

return Objectify(HapFilteredRegularCWComplex,
       rec(
           nrCells:=Y!.nrCells,
           filteredDimension:=M!.filteredDimension,
           boundaries:=Y!.boundaries,
           coboundaries:=Y!.coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=Y!.orientation,
           properties:=Y!.properties));


end;
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(FiltrationTermOfRegularCWComplex,
function(F,t)
local Y, Boundaries, Orientation, n, bool;

Boundaries:=1*F!.boundaries;

if not F!.orientation=fail then bool:=true; else bool:=false; fi;

if bool then
Orientation:=1*F!.orientation;
fi;

for n in [0..Dimension(F)] do
Boundaries[n+1]:=Boundaries[n+1]{[1..F!.filteredDimension(t,n)]};
od;

if bool then
for n in [0..Dimension(F)] do
Orientation[n+1]:=Orientation[n+1]{[1..F!.filteredDimension(t,n)]};
od;
fi;

if bool then
Y:=RegularCWComplex(Boundaries,Orientation);
else
Y:=RegularCWComplex(Boundaries);
fi;

return Y;

end);
#################################################################
#################################################################

