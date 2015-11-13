
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
local
      Y, Start, Contract, nn, dim, bool, BOOL, FREE, degrees,F, CNT;

Y:=arg[1];
if Length(arg)=2 then Start:=arg[2]; else Start:=1; fi;

CNT:=0;
#############################################
##### The work-horse function.###############
Contract:=function(n,F)
local

      b, C, i, j, t, cob, pos, bool, cnt, CoboundCondition,
      Free, mybool, 
      MCoboundaries, U;

#This function removes pairs of n- and (n+1)-cells if possible.
#U=Upper, M=Middle and L=Lower dimensional cells.

####################
####################
if CNT=0 then
Y!.degrees:=[];
for i in [1..Dimension(Y)+1] do
Y!.degrees[i]:=List([1..Y!.filteredDimension(1,i-1)],a->1);
t:=1;
  while t<EvaluateProperty(Y,"filtration_length") do
  t:=t+1;
    for j in [1+Y!.filteredDimension(t-1,i-1)..Y!.filteredDimension(t,i-1)] do
     Y!.degrees[i][j]:=t;
    od;
  od;
od;
fi;
CNT:=1;
####################
####################

MCoboundaries:=Y!.coboundaries[n+1];
C:=Length(MCoboundaries);

#############
#############
CoboundCondition:=function(i)
local cnt;
if Y!.degrees[n+1][i]>F then return [false]; fi;
if F>1 and Y!.degrees[n+1][i]<F then return [false]; fi;
cnt:=MCoboundaries[i]{[2..Length(MCoboundaries[i])]};
cnt:=Filtered(cnt,x->Y!.degrees[n+2][x]<=F);
if Length(cnt)=1 then 
Y!.degrees[n+1][i]:=Y!.degrees[n+1][i]+1;
Y!.degrees[n+2][cnt[1]]:=Y!.degrees[n+2][cnt[1]]+1;
return [true,cnt[1]]; fi;
return [false]; 
end;
#############
#############

#######################
#######################THIS TAKES ALL THE TIME
FREE:=Filtered([Start..Y!.filteredDimension(F,n)], i->Y!.degrees[n+1][i]=F); 

mybool:=false;
for i in FREE do
if  CoboundCondition(i)[1] then mybool:=true;fi;
od;


return mybool;
#######################
#######################

end;
####End of work-horse function.#############
############################################
dim:=EvaluateProperty(Y,"dimension");

for F in [Start..EvaluateProperty(Y,"filtration_length")] do

bool:=true;
BOOL:=true;
nn:=dim-1;

while BOOL or nn>0 do
BOOL:=false;
  for nn in Reversed([0..dim-1]) do
    while bool do
      bool:=Contract(nn,F);
      if bool=true then BOOL:=true; fi;
    od;
    bool:=true;
  od;
od;

od;

end);
############################################
############################################

############################################
############################################
InstallGlobalFunction(ContractedFilteredRegularCWComplex,
function(YY)
local Y, W, bnd,  perm,F, cnt, d, i,j, newboundaries, adjust, f, properties,
filtdim, orien, nrcells;

Y:=StructuralCopyOfFilteredRegularCWComplex(YY);
F:=EvaluateProperty(Y,"filtration_length");
HAPContractFilteredRegularCWComplex(Y);


perm:=[];
cnt:=[];

####
for d in [1..Length(Y!.boundaries)] do
cnt[d]:=List([1..F+1],i->0);
for i in [1..Length(Y!.boundaries[d])] do
cnt[d][Y!.degrees[d][i]]:=1*cnt[d][Y!.degrees[d][i]]+1;
od;
od;
####

####
adjust:=[];
for d in [1..Length(Y!.boundaries)] do
adjust[d]:=[0];
for i in [2..F+1] do
adjust[d][i]:=adjust[d][i-1]+cnt[d][i-1];
od;
od;
####

####
for d in [1..Length(Y!.boundaries)] do
perm[d]:=[];
cnt[d]:=List([1..F+1],i->0);
for i in [1..Length(Y!.boundaries[d])] do
cnt[d][Y!.degrees[d][i]]:=1*cnt[d][Y!.degrees[d][i]]+1;
perm[d][i]:=1*cnt[d][Y!.degrees[d][i]]+adjust[d][Y!.degrees[d][i]];

od;
od;
####





newboundaries:=List([1..Length(Y!.boundaries)],i->[]);
orien:=List([1..Length(Y!.boundaries)],i->[]);
for d in [1..Length(Y!.boundaries)] do
newboundaries[d]:=[];
orien[d]:=[];
for i in [1..Length(Y!.boundaries[d])] do
newboundaries[d][perm[d][i]]:=StructuralCopy(Y!.boundaries[d][i]);
orien[d][perm[d][i]]:=StructuralCopy(Y!.orientation[d][i]);
od;
od;

for d in [2..Length(Y!.boundaries)] do
for i in [1..Length(Y!.boundaries[d])] do
bnd:=newboundaries[d][i];
for j in [2..Length(bnd)] do
bnd[j]:=perm[d-1][bnd[j]];
od;
od;
od;

####################
filtdim:=function(i,d);

return adjust[d+1][i];

end;
####################

for d in [1..Length(newboundaries)] do
newboundaries[d]:=newboundaries[d]{[1..filtdim(F,d-1)]};
od;



W:=RegularCWComplex(newboundaries);
properties:=W!.properties;
Add(properties,["filtration_length",F]);

return Objectify(HapFilteredRegularCWComplex,
       rec(
           nrCells:=W!.nrCells,
           filteredDimension:=filtdim,
           boundaries:=W!.boundaries,
           coboundaries:=W!.coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=orien,
           properties:=properties));


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

