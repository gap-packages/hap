#(C) Graham Ellis

##################################################################
##################################################################
InstallGlobalFunction(SimplicialComplexToRegularCWSpace,
function(K)
local
	NrCells,Boundaries,tmp,TMP,Coboundaries,Properties,
        cnt,b,bb,k,n,s,x,i,j,dim ;

####################
NrCells:=function(n);
if n>dim then return 0; fi;
return Length(Filtered(Boundaries[n+1],x->not x[1]=0));
end;
####################

dim:=Dimension(K);
Properties:=[["dimension",dim]];


### BOUNDARIES BEGIN ######################
Boundaries:=[]; #Boundaries[n+1] contains the info on n-cells
Boundaries[1]:=List([1..K!.nrSimplices(0)],x->[1,0]);

                 ##We denote by 0 the unique vertex in dimension -1.

for n in [1..dim] do

  Boundaries[n+1]:=[];
  tmp:=List(Boundaries[1],x->[]);
  TMP:=List(Boundaries[1],x->[]);
  cnt:=0;
  for s in K!.simplicesLst[n] do
    cnt:=cnt+1;
    Add(tmp[s[1]],s);
    Add(TMP[s[1]],cnt);
  od;
 for k in [1..K!.nrSimplices(n)] do
 bb:=K!.simplices(n,k);

 b:=List(bb,x->  Difference(bb,[x]) );
 Apply(b,x->   TMP[x[1]][Position(tmp[x[1]],x)] );

 Boundaries[n+1][k]:=Concatenation([Length(b)],b);
 od;
od;
Boundaries[dim+2]:=[];
### BOUNDARIES END ###############################

### COBOUNDARIES BEGIN ######################
Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
for n in [0..dim] do
  k:=n+3;
  Coboundaries[n+1]:=List(Boundaries[n+1],i->[0]);
  for j in [1..Length(Boundaries[n+2])] do
    b:=Boundaries[n+2][j]; 
    for i in b{[2..k]} do
      Coboundaries[n+1][i][1]:=Coboundaries[n+1][i][1]+1;
      Add(Coboundaries[n+1][i],j);
    od;
  od;
  #for b in Coboundaries[n+1] do
  #Append(b,List([1..Length(b)-1],a->1)); 
  #od;
od;
Coboundaries[dim+1]:=List(Boundaries[dim+1],a->[0]);
### COBOUNDARIES END ###############################


return Objectify(HapRegularCWSpace,
       rec(
           nrCells:=NrCells,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           properties:=Properties));

end);
##################################################################
##################################################################

#############################################
#############################################
InstallOtherMethod(Dimension,
"Dimension of  regular CW-space",
[IsHapRegularCWSpace],
function(f) return EvaluateProperty(f,"dimension");
return EvaluateProperty(f,"dimension");
end);
#############################################
#############################################

#############################################
#############################################
InstallGlobalFunction(HAPContractRegularCWSpace,
function(Y)
local
      Contract, nn, dim, bool, BOOL;

#############################################
##### The work-horse function.###############
Contract:=function(n)
local

      b, C, i, j, t, cob, pos, bool,
      Free, UBoundaries, UCoboundaries,
      MBoundaries, MCoboundaries, LCoboundaries, U;

#This function removes pairs of n- and (n+1)-cells if possible.
#U=Upper, M=Middle and L=Lower dimensional cells.

####################
####################
if Y!.vectorField=fail then
Y!.vectorField:=List([1..Dimension(Y)],i->[]);
Y!.inverseVectorField:=List([1..Dimension(Y)],i->[]);
Y!.bnd:=StructuralCopy(Y!.boundaries);
Y!.cobnd:=StructuralCopy(Y!.coboundaries);
fi;
####################
####################

MCoboundaries:=Y!.cobnd[n+1];
MBoundaries:=Y!.bnd[n+1];
UCoboundaries:=Y!.cobnd[n+2];
UBoundaries:=Y!.bnd[n+2];
if n>0 then
  LCoboundaries:=Y!.cobnd[n];
fi;
C:=Length(MCoboundaries);

#######################
#######################THIS TAKES ALL THE TIME
Free:=Filtered([1..C],i->MCoboundaries[i][1]=1);
if Length(Free)=0 then return false;fi;
#######################
#######################


for i in Free do
if MCoboundaries[i][1]=1 then
Y!.vectorField[n+1][MCoboundaries[i][2]]:=i;
Y!.inverseVectorField[n+1][i]:=MCoboundaries[i][2];
      ###
  if n>0 then
    b:=MBoundaries[i];
    for j in StructuralCopy(b{[2..1+b[1]]}) do
     t:=LCoboundaries[j][1];
     LCoboundaries[j][1]:=LCoboundaries[j][1]-1;
     cob:=LCoboundaries[j];
     pos:=Position(cob{[2..t+1]},i);
     LCoboundaries[j]:=Concatenation(cob{[1..pos]},cob{[2+pos..Length(cob)]});
    od;
  fi;
      ###
    U:=MCoboundaries[i][2];
    b:=UBoundaries[U];
    for j in StructuralCopy(b{[2..1+b[1]]}) do
     t:=MCoboundaries[j][1];
     MCoboundaries[j][1]:=MCoboundaries[j][1]-1;
if t=2 then Add(Free,j);fi;############################ADDED
     cob:=MCoboundaries[j];
     pos:=Position(cob{[2..t+1]},U);
     MCoboundaries[j]:=Concatenation(cob{[1..pos]},cob{[2+pos..Length(cob)]});
    od;
      ###
  MBoundaries[i]:=[0];
  UBoundaries[U]:=[0];
  UCoboundaries[U]:=[0];
  MCoboundaries[i]:=[0];
fi;
od;

Y!.bnd[n+2]:=UBoundaries;
Y!.cobnd[n+2]:=UCoboundaries;
Y!.bnd[n+1]:=MBoundaries;
Y!.cobnd[n+1]:=MCoboundaries;
if n>0 then
  Y!.cobnd[n]:=LCoboundaries;
fi;

Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.bnd[k+1],x->not x[1]=0));
            end;

if Length(Free)>0 then return true;
else return false; fi;

end;
####End of work-horse function.#############
############################################

dim:=EvaluateProperty(Y,"dimension");
bool:=true;
BOOL:=true;
nn:=dim-1;

while BOOL or nn>0 do
BOOL:=false;
  for nn in Reversed([0..dim-1]) do
    while bool do
      bool:=Contract(nn);
      if bool=true then BOOL:=true; fi;
    od;
    bool:=true;
  od;
od;

end);
############################################
############################################

#####################################################################
#####################################################################
InstallOtherMethod(Size,
"Volume of a regular CW-space",
[IsHapRegularCWSpace],
function(Y) return Sum(List( [1..Length(Y!.boundaries)],i->Y!.nrCells(i-1)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(HAPRemoveCellFromRegularCWSpace,
function(Y,dim,n)
local  bnd, x,tmp;

#Remove the n-th cell in dimension dim

####################
####################
if Y!.vectorField=fail then
Y!.vectorField:=List([1..Dimension(Y)],i->[]);
Y!.inverseVectorField:=List([1..Dimension(Y)],i->[]);
Y!.bnd:=StructuralCopy(Y!.boundaries);
Y!.cobnd:=StructuralCopy(Y!.coboundaries);
fi;
####################
####################


dim:=dim+1;

if dim=1 then Y!.bnd[1][n]:=[0]; return [dim-1,n]; fi;

bnd:=Y!.bnd[dim][n];
bnd:=bnd{[2..Length(bnd)]};
Y!.bnd[dim][n]:=[0];

if dim>1 then
for x in bnd do
tmp:=Y!.cobnd[dim-1][x];
tmp[1]:=tmp[1]-1;
#tmp[Length(tmp)]:=-42;
tmp[Position(tmp{[2..Length(tmp)]},n)+1]:=-42;

tmp:=Filtered(tmp,i->not i = -42);

Y!.cobnd[dim-1][x]:=tmp;
if IsBound(Y!.free) then
if IsBound(Y!.free[dim-1]) then
if tmp[1]=1 then AddSet(Y!.free[dim-1],x); fi;
fi;
fi;
od;
fi;

return [dim-1,n];

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(CriticalCellsOfRegularCWSpace,
function(arg)
local Y,ContractSpace,cells,dim,c,pos,ppos;

Y:=arg[1];
if not Y!.criticalCells=fail then
return Y!.criticalCells;
fi;

if Length(arg)>1 then
   ContractSpace:=HAPContractRegularCWSpace_Alt;
else
   ContractSpace:=HAPContractRegularCWSpace;
fi;


#######
dim:=0;
while true do
if Y!.nrCells(dim)=0 then break; fi;
dim:=dim+1;
od;
dim:=dim-1;
#######

cells:=[];
ContractSpace(Y);

while true do

  if 
#Size(Y)=0  
Sum(List( [1..Length(Y!.bnd)],i->Y!.nrCells(i-1)))=0
then  
Y!.criticalCells:=cells; 
Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
Unbind(Y!.bnd);
Unbind(Y!.cobnd);
return cells; fi;

  pos:=0;

  while true do
    pos:=pos+1;
    ppos:=PositionProperty(Y!.bnd[dim+1]{[pos..Length(Y!.bnd[dim+1])]},                          x->x[1]>0);    if ppos=fail then dim:=dim-1; break; fi;

  pos:=pos+ppos-1;

    c:=HAPRemoveCellFromRegularCWSpace(Y,dim,pos);

    Add(cells,c);

    ContractSpace(Y);

  od;
od;

Y!.criticalCells:=cells;
Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
Unbind(Y!.bnd);
Unbind(Y!.cobnd);

return cells;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(CubicalComplexToRegularCWSpace,
function(M)
local C, Properties, Boundaries, Coboundaries, BinLst, 
      LstBin, dim,  bnd,  Boundary, ArrayValueDim, Orientation, 
      Dimension,   n, i, j, k, b, v;

dim:=EvaluateProperty(M,"dimension");
ArrayValueDim:=ArrayValueFunctions(dim);

C:=ChainComplex(M);
BinLst:=C!.coordinateToPosition;
LstBin:=C!.positionToCoordinate;

Properties:=[["dimension",dim]];

Dimension:=C!.dimension;


#######################################
Boundary:=function(n,j)
local x,poscells,negcells,nn,a,b,cnt;
poscells:=[];
negcells:=[];

cnt:=0;
nn:=LstBin[n+1][j];
for x in [1..Length(nn)] do
  if IsEvenInt(nn[x]) then
    cnt:=cnt+1;
    a:=StructuralCopy(nn);
    a[x]:=a[x]+1;
    b:=StructuralCopy(nn);
    b[x]:=b[x]-1;
    if IsOddInt(cnt) then
        Add(poscells,a);
        Add(negcells,b);
    else
        Add(poscells,b);
        Add(negcells,a);
    fi;
  fi;
od;

Apply(poscells,x->ArrayValueDim(BinLst,x));
Apply(negcells,x->ArrayValueDim(BinLst,x));
return [poscells,negcells];
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


return Objectify(HapRegularCWSpace,
       rec(
           nrCells:=Dimension,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=Orientation,
           properties:=Properties));


end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(ChainComplexOfRegularCWSpace,
function(Y)
local
	C, Dimension, Boundary, one, zero, n, dim, characteristic;

Dimension:=Y!.nrCells;
dim:=EvaluateProperty(Y,"dimension");

zero:=[];
for n in [1..dim+1] do
zero[n]:=List([1..Dimension(n-1)],i->0);
od;

if not IsBound(Y!.orientation) then
characteristic:=2;
one:=One(GF(2));
######################
Boundary:=function(n,k)
local b,i,j,B;

b:=StructuralCopy(zero[n]);
B:=Y!.boundaries[n+1][k];

for i in [2..Length(B)] do
b[B[i]]:=1;
od;

return one*b;
end;
######################
else 
characteristic:=0;
######################
Boundary:=function(n,k)
local b,i,j,B,sn;

b:=StructuralCopy(zero[n]);
B:=Y!.boundaries[n+1][k];
sn:=Y!.orientation[n+1][k];

for i in [2..Length(B)] do
b[B[i]]:=sn[i-1];
od;

return b;
end;
######################
fi;

return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimension,
           boundary:=Boundary,
           properties:=[
           ["length",dim],
           ["type","chainComplex"],
           ["characteristic",characteristic]]
           ));


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(ChainComplexOfRegularCWSpaceWithVectorField,
function(Y)
local
        basis, bij,Dimension, Boundary, one, zero, b, n, dim, characteristic, DeformCell,DeformCellSgn;

dim:=EvaluateProperty(Y,"dimension");

basis:=[];
bij:=[];
for n in [0..dim] do
basis[n+1]:=Filtered(CriticalCellsOfRegularCWSpace(Y),x->x[1]=n);
Apply(basis[n+1],x->x[2]);
bij[n+1]:=[];
for b in [1..Length(basis[n+1])] do
bij[n+1][basis[n+1][b]]:=b;
od;
od;

###############################
Dimension:=function(n);
return Length(basis[n+1]);
end;
###############################


zero:=[];
for n in [1..dim+1] do
zero[n]:=List([1..Dimension(n-1)],i->0);
od;


###############################
###############################
DeformCell:=function(n,k)
local x,f,bnd,def;            #This will return a list of n-cells 
                              #into which the k-th n-cell is deformed.

if [n,k] in Y!.criticalCells then
return [k];
fi;
if n>0 then
if IsBound(Y!.vectorField[n][k]) then return []; fi;
fi;
 
f:=Y!.inverseVectorField[n+1][k];
bnd:=Y!.boundaries[n+2][f];
def:=[];
for x in [2..Length(bnd)] do
if not bnd[x]=k then
Append(def,DeformCell(n,bnd[x]));
fi;
od;

return def;
end;
###############################
###############################

###############################
###############################
DeformCellSgn:=function(n,kk)
local sgnn,x,f,k,sgnk,cnt,bnd,def,sn,tog,def1,def2;            
		        	#This will return a list of signed n-cells
                                #into which the k-th n-cell is deformed.

k:=AbsInt(kk);
sgnk:=SignInt(kk);
if [n,k] in Y!.criticalCells then
return [kk];
fi;

if n>0 then
if IsBound(Y!.vectorField[n][k]) then return []; fi;
fi;

f:=Y!.inverseVectorField[n+1][k];
bnd:=Y!.boundaries[n+2][f];
sn:=Y!.orientation[n+2][f];

def:=[]; def1:=[];def2:=[];

for x in [2..Length(bnd)] do
if not bnd[x]=k then
Add(def1,sn[x-1]*bnd[x]);
else
sgnn:=sn[x-1];
break;
fi;
od;
cnt:=x+1;

for x in [cnt..Length(bnd)] do
Add(def2,sn[x-1]*bnd[x]);
od;

if sgnn=1 then 
def:=-Concatenation(Reversed(def1),Reversed(def2));
else
def:=Concatenation(def2,def1);
fi;

Apply(def,x->DeformCellSgn(n,x));

if sgnk=1 then return Flat(def);
else
return -Reversed(Flat(def));
fi;
end;
###############################
###############################


if not IsBound(Y!.orientation) then
characteristic:=2;
one:=One(GF(2));
######################
Boundary:=function(n,k)
local b,i,j,B;

b:=StructuralCopy(zero[n]);
B:=Y!.boundaries[n+1][basis[n+1][k]];
B:=B{[2..Length(B)]};
Apply(B,x->DeformCell(n-1,x));
B:=Concatenation(B);
Apply(B,i->bij[n][i]);

for i in B do
b[i]:=b[i]+1;
od;

return one*b;
end;
######################
else
characteristic:=0;
######################
Boundary:=function(n,k)
local b,i,j,B,sn;

b:=StructuralCopy(zero[n]);
B:=Y!.boundaries[n+1][basis[n+1][k]];
B:=B{[2..Length(B)]};
sn:=Y!.orientation[n+1][basis[n+1][k]];
B:=List([1..Length(B)],i->sn[i]*B[i]);

Apply(B,x->DeformCellSgn(n-1,x));
B:=Concatenation(B);
Apply(B,i->SignInt(i)*bij[n][AbsInt(i)]);


for i in B do
#b[AbsInt(i)]:=b[AbsInt(i)]+SignInt(i)*sn[AbsInt(i)];
b[AbsInt(i)]:=b[AbsInt(i)]+SignInt(i);
od;

return b;
end;
######################
fi;

if IsBound(Y!.orientation) then DeformCell:=DeformCellSgn; fi;
return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimension,
           boundary:=Boundary,
           deform:=DeformCell,
           basis:=basis,
           bij:=bij,
           properties:=[
           ["length",dim],
           ["type","chainComplex"],
           ["characteristic",characteristic]]
           ));


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod( ChainComplex,
"for regular CW spaces, using discrete vector fields",
 [IsHapRegularCWSpace],
 function(Y)
 CriticalCellsOfRegularCWSpace(Y);
 return ChainComplexOfRegularCWSpaceWithVectorField(Y);
 end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod( Homology,
"Homology of a regular CW spaces, using discrete vector fields",
 [IsHapRegularCWSpace,IsInt],
 function(Y,n) local C;
 if not IsBound(Y!.orientation) then
 Print("Can only compute the mod 2 homology as no orientation is available.\n");
 fi;
 C:=ChainComplex(Y);
 return Homology(C,n);
 end);
##########################################################
##########################################################


