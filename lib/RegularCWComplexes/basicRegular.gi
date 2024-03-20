#(C) Graham Ellis

######################################
######################################
InstallGlobalFunction(HomotopyTruncation,
function(W,N)
local Y, V, A, n, x, i;

if Dimension(W)=N then
return ContractedComplex(W);
fi;

Y:=1*W!.boundaries{[1..N+2]};
Add(Y,[]);
Y:=RegularCWComplex(Y);
CocriticalCellsOfRegularCWComplex(Y,N+1);
V:=SSortedList(Flat(Y!.vectorField[N+1]));
V:=Filtered([1..Length(Y!.boundaries[N+1])],i-> not i in V);

Y!.boundaries[N+1]:=Y!.boundaries[N+1]{V};
Y!.boundaries[N+2]:=[];

return ContractedComplex(RegularCWComplex(Y!.boundaries));

end);
######################################
######################################

##################################################################
##################################################################
InstallGlobalFunction(SimplicialComplexToRegularCWComplex_alt,
function(arg)
local
	K,DM,NrCells,Boundaries,tmp,TMP,Coboundaries,Properties,
        Orientation, 
        cnt,ln,b,bb,k,n,s,x,i,j,dim,bbb ;

K:=IntegerSimplicialComplex(arg[1]);
if Length(arg)>1 then dim:=arg[2]; else dim:=Dimension(K); fi; 

####################
NrCells:=function(n);
if n>dim then return 0; fi;
return Length(Filtered(Boundaries[n+1],x->not x[1]=0));
end;
####################

#dim:=Dimension(K);
Properties:=[["dimension",dim]];

#############################
Orientation:=[];
#Orientation[1]:=ListWithIdenticalEntries(K!.nrSimplices(0),[1]);
Orientation[1]:=List([1..K!.nrSimplices(0)],i->[1]); #more memory but safer!!
for n in [1..dim] do
  tmp:=[];  
  for i in [1..n+1] do
  Add(tmp,(-1)^(i+1));
  od;
#Orientation[n+1]:=ListWithIdenticalEntries(K!.nrSimplices(n),tmp);
Orientation[n+1]:=List([1..K!.nrSimplices(n)],i->StructuralCopy(tmp));
od;

#############################

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
 bb:=SSortedList(bb);

 #b:=List(bb,x->  Difference(bb,[x]) )#;
 b:=[]; ln:=n+1; #ln:=Length(bb);
 
 for i in [1..ln] do
 bbb:=1*bb;
 Remove(bbb,i);
 b[i]:=bbb;
 od;
 Apply(b,x->   TMP[x[1]][Position(tmp[x[1]],x)] );  

 #Boundaries[n+1][k]:=Concatenation([Length(b)],b);
 Add(b,ln,1);
 Boundaries[n+1][k]:=b;

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


return Objectify(HapRegularCWComplex,
       rec(
           nrCells:=NrCells,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           orientation:=Orientation,
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
"Dimension of  regular CW-complex",
[IsHapRegularCWComplex],
function(f) return EvaluateProperty(f,"dimension");
return EvaluateProperty(f,"dimension");
end);
#############################################
#############################################

#############################################
#############################################
InstallGlobalFunction(HAPContractRegularCWComplex,
function(arg)
local
      Y,d,Contract, nn, dim, bool, BOOL, FREE;

#NB: to obtain the algorithm described in the Alzobydi-Ellis paper
#under Option 1 you need to comment out lines A and C and uncomment line B.
#The default setting is for Option 2 of the paper.

Y:=arg[1];
if Length(arg)>1 then d:=arg[2];  #Added 12/03/2018
else d:=0; fi;

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
if not IsBound(FREE) then FREE:=[1..C]; fi;

Free:=[];
for i in FREE do
if MCoboundaries[i][1]=1 then Add(Free,i);fi;           #Line A
#if MCoboundaries[i][1]=1 then Add(Free,i); break; fi;  #Line B
od;


if Length(Free)=0 then Unbind(FREE); return false;fi;
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
if t=2 then Add(Free,j);fi;############################ADDED  #Line C
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

if Length(Free)>0 then FREE:=Free; return true;
else Unbind(FREE); return false; fi;

end;
####End of work-horse function.#############
############################################

dim:=EvaluateProperty(Y,"dimension");
bool:=true;
BOOL:=true;
nn:=dim-1;

while BOOL or nn>d do
BOOL:=false;
  for nn in Reversed([d..dim-1]) do
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
"Volume of a regular CW-complex",
[IsHapRegularCWComplex],
function(Y) return Sum(List( [1..Length(Y!.boundaries)],i->Y!.nrCells(i-1)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(HAPRemoveCellFromRegularCWComplex,
function(Y,dim,n)
local  bnd, x,tmp, cobnd;

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

bnd:=Y!.bnd[dim][n];
bnd:=bnd{[2..Length(bnd)]};
Y!.bnd[dim][n]:=[0];

cobnd:=Y!.cobnd[dim][n];          ####Added this loop July 2012
cobnd:=cobnd{[2..Length(cobnd)]};                             #
for x in cobnd do                                             # 
tmp:=Y!.bnd[dim+1][x];                                        #
tmp[1]:=tmp[1]-1;                                             #
tmp[Position(tmp{[2..Length(tmp)]},n)+1]:=-42;                #
tmp:=Filtered(tmp,i->not i = -42);                            #
Y!.bnd[dim+1][x]:=tmp;                                        #
od;                                              ############## 


if dim=1 then  return [dim-1,n]; fi;

for x in bnd do
tmp:=Y!.cobnd[dim-1][x];
tmp[1]:=tmp[1]-1;
tmp[Position(tmp{[2..Length(tmp)]},n)+1]:=-42;

tmp:=Filtered(tmp,i->not i = -42);

Y!.cobnd[dim-1][x]:=tmp;
if IsBound(Y!.free) then
if IsBound(Y!.free[dim-1]) then
if tmp[1]=1 then AddSet(Y!.free[dim-1],x); fi;
fi;
fi;
od;

return [dim-1,n];

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(CriticalCellsOfRegularCWComplex,
function(arg)
local Y,ContractSpace,cells,dim,c,pos,ppos,  b,x, bbooll;

Y:=arg[1];
if not Y!.criticalCells=fail then
return Y!.criticalCells;
fi;

if Dimension(Y)=0 then 
Y!.criticalCells:=List([1..Y!.nrCells(0)],i->[0,i]);
return Y!.criticalCells;
fi;

##############################
if Length(arg)>1 then 
cells:=CocriticalCellsOfRegularCWComplex(Y,arg[2]);
if arg[2]<EvaluateProperty(Y,"dimension") then 
Y!.criticalCells:=fail;  fi;
return cells;
fi;
##############################


   ContractSpace:=HAPContractRegularCWComplex;


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


#######
#######
if dim=0 then bbooll:=true; else
bbooll:=false;
for b in Y!.bnd[dim+1][pos]{[2..Length(Y!.bnd[dim+1][pos])]} do
if bbooll then break; fi;
if Y!.cobnd[dim][b][1]=2 then bbooll:=true; break; fi;
od;
fi;
#######
#######
 
    c:=HAPRemoveCellFromRegularCWComplex(Y,dim,pos);

    Add(cells,c);

if bbooll then
    ContractSpace(Y);
fi;

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
InstallGlobalFunction(CubicalComplexToRegularCWComplex,
function(arg)
local M,dim,C, Properties, Boundaries, Coboundaries, BinLst, 
      LstBin,   bnd,  Boundary, ArrayValueDim, Orientation, 
      Dimension,  dimm,  n, i, j, k, b, v;

M:=arg[1];
if Length(arg)>1 then dim:=arg[2];
else
dim:=EvaluateProperty(M,"dimension");
fi;
ArrayValueDim:=ArrayValueFunctions(EvaluateProperty(M,"dimension"));

C:=ChainComplex(M);
BinLst:=C!.coordinateToPosition;
LstBin:=C!.positionToCoordinate;

dimm:=Position(List([0..dim],i->C!.dimension(i)),0);
if dimm=fail then dimm:=dim; else dimm:=dimm-2; fi;
Properties:=[["dimension", dimm ]];

if Length(arg)=1 then
Dimension:=C!.dimension;
else
###
Dimension:=function(n);
  if n<dim then return C!.dimension(n);
  else return 0; fi;
end;
#######
fi;


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


return Objectify(HapRegularCWComplex,
       rec(
           nrCells:=Dimension,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           vectorField:=fail,
           inverseVectorField:=fail,
           criticalCells:=fail,
           orientation:=Orientation,
           coordinateToPosition:=BinLst,
           positionToCoordinate:=LstBin,
           properties:=Properties));


end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallGlobalFunction(ChainComplexWithChainHomotopy,
function(Y)
local C, n, k, chainHomotopy;

#Inputs a regular CW complex Y. It creates a discrete vector field on Y.
#Let Y' be the homotopy equivalent space with cells correponding to the
#critical cells of Y. It outputs the chain complex C=C(Y') with a component
#C!.chainHomotopy which is a chain homotopy on the bigger chain complex C(Y)
#corresponding to the homotopy equivalence C(Y) -->C(Y').

CriticalCells(Y);
C:=ChainComplexOfRegularCWComplexWithVectorField(Y,true);
for n in [0..Dimension(Y)] do
for k in [1..Y!.nrCells(n)] do
C!.deform(n,k);
od;
od;

##############################
chainHomotopy:=function(w,n)
local u, k, h, a, abs;
u:=[1..Y!.nrCells(n+1)]*0;

for k in [1..Length(w)] do
if not w[k]=0 then
h:=C!.htpy[n+1][k];
for a in h do
abs:=AbsInt(a);
u[abs]:=u[abs]+w[k]*SignInt(a);
od;
fi;
od;
return u;
end;
##############################

C!.chainHomotopy:=chainHomotopy;
return C;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallGlobalFunction(ChainComplexOfRegularCWComplex,
function(Y)
local
	C, Dimension, Boundary, one, zero, n, dim, characteristic;


dim:=EvaluateProperty(Y,"dimension");
##########################
Dimension:=function(n);
if n<0 or n>dim then return 0; fi;
return Length(Y!.boundaries[n+1]);
end;
##########################


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

if n>dim then return [one]; fi;

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

#if Dimension(n)=0  then return [1..Dimension(n-1)]*0; fi;
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
InstallGlobalFunction(ChainComplexOfRegularCWComplexWithVectorField,
function(arg)
local
        Y,basis, bool, bij,Dimension, Boundary, one, zero, b, n, dim, 
        characteristic, DeformCellRec, DeformCell,DeformCellSgn, 
        DeformCellSgnHtpy, HomotopicalDeformCell, 
        HDCrec, DCSrec, DCSHrec, AlgRed, HtpyRed, BoundaryRec;


Y:=arg[1];

#############################################
HtpyRed:=function(w)
local cnt, tog;

cnt:=1;
tog:=true;

while tog do
tog:=false;
while cnt<Length(w) do
if w[cnt]=-w[cnt+1] then
w[cnt]:=0; w[cnt+1]:=0;
cnt:=cnt+2; tog:=true;
else cnt:=cnt+1;
fi;
od;
w:=Filtered(w,i->not i=0);
od;
return w;
end;
#############################################

AlgRed:=AlgebraicReduction;


dim:=EvaluateProperty(Y,"dimension");

basis:=[];
bij:=[];
for n in [0..dim] do
basis[n+1]:=Filtered(CriticalCellsOfRegularCWComplex(Y),x->x[1]=n);
Apply(basis[n+1],x->x[2]);
bij[n+1]:=[];
for b in [1..Length(basis[n+1])] do
bij[n+1][basis[n+1][b]]:=b;
od;
od;

###############################
Dimension:=function(n);
if IsBound(basis[n+1]) then
return Length(basis[n+1]);
else
return 0; fi;
end;
###############################


zero:=[];
for n in [1..dim+1] do
zero[n]:=List([1..Dimension(n-1)],i->0);
od;

DeformCellRec:=List([0..EvaluateProperty(Y,"dimension")],i->[]); #Added 28 June 2021
###############################
###############################
DeformCell:=function(n,k)
local x,f,bnd,def;            #This will return a list of n-cells 
                              #into which the k-th n-cell is deformed.

if IsBound(DeformCellRec[n+1][k]) then
return DeformCellRec[n+1][k];
fi;


if [n,k] in Y!.criticalCells then
DeformCellRec[n+1][k]:=[k];
return [k];
fi;

if n>0 then
if IsBound(Y!.vectorField[n][k]) then 
DeformCellRec[n+1][k]:=[];
return []; fi;
fi;
 
f:=Y!.inverseVectorField[n+1][k];
bnd:=Y!.boundaries[n+2][f];
def:=[];
for x in [2..Length(bnd)] do
if not bnd[x]=k then
Append(def,DeformCell(n,bnd[x]));
fi;
od;

DeformCellRec[n+1][k]:=def; 
return def;
end;
###############################
###############################

DCSrec:=List([1..dim+1],i->[]);;
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

if IsBound(DCSrec[n+1][k]) then
if sgnk=1 then return DCSrec[n+1][k];
else
return -DCSrec[n+1][k];
fi;
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
def:=-Concatenation(def1,def2);
else
def:=Concatenation(def2,def1);
fi;

Apply(def,x->DeformCellSgn(n,x));

def:=Flat(def);
Apply(def,x->[x,0]);

def:=AlgRed(def);

Apply(def,x->x[1]);

DCSrec[n+1][k]:=def;

if sgnk=1 then return def;
else
return -def;
fi;
end;
###############################
###############################

DCSHrec:=List([1..dim+1],i->[]);;
###############################
###############################
DeformCellSgnHtpy:=function(n,kk)
local sgnn,x,f,k,sgnk,cnt,bnd,def,defcp,sn,tog,def1,def2;
                                #This will return a list of signed n-cells
                                #into which the k-th n-cell is deformed.
k:=AbsInt(kk);
sgnk:=SignInt(kk);
if [n,k] in Y!.criticalCells then
DCSHrec[n+1][k]:=[];
return [kk];
fi;

if n>0 then
if IsBound(Y!.vectorField[n][k]) then 
DCSHrec[n+1][k]:=[];
return []; fi;
fi;

if IsBound(DCSrec[n+1][k]) and IsBound(DCSHrec[n+1][k]) then
if sgnk=1 then return DCSrec[n+1][k];
else
return -DCSrec[n+1][k];
fi;
fi;

f:=Y!.inverseVectorField[n+1][k];

bnd:=Y!.boundaries[n+2][f];
sn:=Y!.orientation[n+2][f];          #CHECK THIS!!

def:=[]; 

for x in [2..Length(bnd)] do
if not bnd[x]=k then
Add(def,sn[x-1]*bnd[x]);
else
sgnn:=sn[x-1];
break;
fi;
od;
cnt:=x+1;

DCSHrec[n+1][k]:=[sgnn*f];

for x in [cnt..Length(bnd)] do
Add(def,sn[x-1]*bnd[x]);
od;

def:=sgnn*def;

defcp:=StructuralCopy(def);
Apply(def,xx->-DeformCellSgnHtpy(n,xx));
for x in defcp do
Append(DCSHrec[n+1][k], -SignInt(x)*DCSHrec[n+1][AbsInt(x)]);
od;

def:=Flat(def);
Apply(def,x->[x,0]);

def:=AlgRed(def);

Apply(def,x->x[1]);

DCSrec[n+1][k]:=def;

if sgnk=1 then return def;
else
return -def;
fi;
end;
###############################
###############################


HDCrec:=List([1..dim+1],i->[]);;
###############################
###############################
HomotopicalDeformCell:=function(n,kk)
local sgnn,x,y,f,k,sgnk,cnt,bnd,def,sn,tog,def1,def2,B,S,i,ii,jj,indx,A,C;
                                #This will return an ordered list of signed 
				#n-cells into which the k-th n-cell is 
				#deformed.

k:=AbsInt(kk);
sgnk:=SignInt(kk);
if [n,k] in Y!.criticalCells then
return [kk];
fi;

if n>0 then
if IsBound(Y!.vectorField[n][k]) then return []; fi;
fi;

if IsBound(HDCrec[n+1][k]) then
if sgnk=1 then return HDCrec[n+1][k];
else
return -Reversed(HDCrec[n+1][k]);
fi;
fi;

f:=Y!.inverseVectorField[n+1][k];
bnd:=Y!.boundaries[n+2][f];
sn:=1*Y!.homotopyOrientation[n+2][f];  ##

if n=1 then
########################################
# Let's order the elements of bnd.
B:=1*bnd{[2..Length(bnd)]};
S:=[1]; 
indx:=[2..Length(B)];
for ii in [1..Length(B)-1] do
   A:=Y!.boundaries[2][B[S[ii]]];
   A:=A{[2,3]};
   for jj in indx do
      C:=Y!.boundaries[2][B[jj]];
      C:=C{[2,3]};
      if Length(Intersection(A,C))>0 then
         Add(S,jj);
         RemoveSet(indx,jj);
         break;
      fi;
   od;
od;
bnd:=Concatenation([Length(B)],B{S});
sn:=sn{S};;


if Length(bnd)>3 then
###### This is a really thoughtless way to get the signs right!!!
###### And it is also wasteful of time.
A:=[Y!.boundaries[2][bnd[2]]{[2,3]}];
for i in [3..Length(bnd)]  do
x:=Y!.boundaries[2][bnd[i]]{[2,3]}; y:=A[Length(A)];;
if x[1] in y then sn[i-1]:=Y!.homotopyOrientation[2][bnd[i]][1]; Add(A,x{[1,2]});  
else  sn[i-1]:=-Y!.homotopyOrientation[2][bnd[i]][1]; Add(A,x{[2,1]}); fi; 

od;
if A[1][1] in A[2] then sn[1]:=-Y!.homotopyOrientation[2][bnd[1]][1]; 
else sn[1]:=Y!.homotopyOrientation[2][bnd[1]][1]; fi;
######
######
########################################
fi;
fi;


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

Apply(def,x->HomotopicalDeformCell(n,x));


def:=Flat(def);

def:=HtpyRed(def);

HDCrec[n+1][k]:=def;

if sgnk=1 then return def;
else
return -Reversed(def);
fi;
end;
###############################
###############################

BoundaryRec:=List([1..dim+1],i->[]);



if not IsBound(Y!.orientation) then
characteristic:=2;
one:=One(GF(2));
######################
Boundary:=function(n,k)
local b,i,j,B;

if IsBound(BoundaryRec[n+1][AbsInt(k)]) then return 
SignInt(k)*BoundaryRec[n+1][AbsInt(k)]; fi;

b:=StructuralCopy(zero[n]);
B:=Y!.boundaries[n+1][basis[n+1][k]];
B:=B{[2..Length(B)]};
Apply(B,x->DeformCell(n-1,x));
B:=Concatenation(B);
Apply(B,i->bij[n][i]);

for i in B do
b[i]:=b[i]+1;
od;

BoundaryRec[n+1][k]:=one*b;
return 1*BoundaryRec[n+1][k];
end;
######################
else
characteristic:=0;
######################
Boundary:=function(n,k)
local b,i,j,B,sn;

if IsBound(BoundaryRec[n+1][AbsInt(k)]) then return
SignInt(k)*BoundaryRec[n+1][AbsInt(k)]; fi;

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

BoundaryRec[n+1][k]:=b;
return 1*BoundaryRec[n+1][k];
end;
######################
fi;

if IsBound(Y!.orientation) then DeformCell:=DeformCellSgn; fi;
if Length(arg)=2 then DeformCell:=DeformCellSgnHtpy; fi;
if Length(arg)=1 then DCSHrec:=fail; fi;


return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimension,
           boundary:=Boundary,
           deform:=DeformCell,
           deformCellSgn:=DeformCellSgn,
           htpy:=DCSHrec,
           homotopicalDeform:=HomotopicalDeformCell,
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

###########################################################
##########################################################
InstallOtherMethod( RegularCWComplex,
"for regular CW spaces, do nothing",
 [IsHapRegularCWComplex],
 function(Y)
 return Y;
 end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod( ChainComplex,
"for regular CW spaces, using discrete vector fields",
 [IsHapRegularCWComplex],
 function(Y)
 CriticalCellsOfRegularCWComplex(Y);
 return ChainComplexOfRegularCWComplexWithVectorField(Y);
 end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod( CochainComplex,
"for regular CW spaces, using discrete vector fields",
 [IsHapRegularCWComplex],
 function(Y) local C;
 CriticalCellsOfRegularCWComplex(Y);
 C:=ChainComplexOfRegularCWComplexWithVectorField(Y);
 return HomToIntegers(C);
 end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod( SparseChainComplex,
"for regular CW spaces, using discrete vector fields",
 [IsHapRegularCWComplex],
 function(Y)
 CriticalCellsOfRegularCWComplex(Y);
 return SparseChainComplexOfRegularCWComplexWithVectorField(Y);
 end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod( Homology,
"Homology of a regular CW spaces, using discrete vector fields",
 [IsHapRegularCWComplex,IsInt],
 function(Y,n) local C, H, m, bool;
 if not IsBound(Y!.orientation) then
 OrientRegularCWComplex(Y);
 fi;
 m:=Minimum(n+1,Dimension(Y));
 bool:=Y!.vectorField=fail or Y!.criticalCells=fail;
 if bool then 
    if m=Dimension(Y) then CriticalCellsOfRegularCWComplex(Y); 
    else
    CocriticalCellsOfRegularCWComplex(Y,m); fi;
 fi;
 C:=ChainComplex(Y);
 H:=Homology(C,n);
 if m<Dimension(Y) and bool then 
 Y!.vectorField:=fail; 
 Y!.criticalCells:=fail; 
 Y!.properties:=Filtered(Y!.properties,x->not x[1]="codim");
 fi;
return H;
 end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod( Cohomology,
"Coomology of a regular CW spaces, using discrete vector fields",
 [IsHapRegularCWComplex,IsInt],
 function(Y,n) local C, H, m, bool;
 if not IsBound(Y!.orientation) then
 OrientRegularCWComplex(Y);
 fi;
 m:=Minimum(n+1,Dimension(Y));
 bool:=Y!.vectorField=fail or Y!.criticalCells=fail;
 if bool then
    if m=Dimension(Y) then CriticalCellsOfRegularCWComplex(Y);
    else
    CocriticalCellsOfRegularCWComplex(Y,m); fi;
 fi;
 C:=ChainComplex(Y);
 H:=Cohomology(HomToIntegers(C),n);
 if m<Dimension(Y) and bool then
 Y!.vectorField:=fail;
 Y!.criticalCells:=fail;
 Y!.properties:=Filtered(Y!.properties,x->not x[1]="codim");
 fi;
return H;
 end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(SparseChainComplexOfRegularCWComplex,
function(Y)
local
        C, Dimension, Boundary, one,  n, dim, characteristic;

#Dimension:=Y!.nrCells;
##########################
Dimension:=function(n);
if n<0 then return 0; fi;
return Length(Y!.boundaries[n+1]);
end;
##########################

dim:=EvaluateProperty(Y,"dimension");


if not IsBound(Y!.orientation) then
characteristic:=2;
one:=One(GF(2));
######################
Boundary:=function(n,k)
local b,i,j,B;

b:=[];
B:=Y!.boundaries[n+1][k];

for i in [2..Length(B)] do
Add(b,[B[i],one]);
od;

return b;
end;
######################
else
characteristic:=0;
######################
Boundary:=function(n,k)
local b,i,j,B,sn;

b:=[];
B:=Y!.boundaries[n+1][k];
sn:=Y!.orientation[n+1][k];

for i in [2..Length(B)] do
Add(b,[B[i],sn[i-1]]);
od;

return b;
end;
######################
fi;

return
Objectify(HapSparseChainComplex,
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
InstallGlobalFunction(SparseChainComplexOfRegularCWComplexWithVectorField,
function(Y)
local
        basis, bij,Dimension, Boundary, one, zero, b, n, dim, characteristic, DeformCell,DeformCellSgn;

dim:=EvaluateProperty(Y,"dimension");

basis:=[];
bij:=[];
for n in [0..dim] do
basis[n+1]:=Filtered(CriticalCellsOfRegularCWComplex(Y),x->x[1]=n);
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

b:=[];
B:=Y!.boundaries[n+1][basis[n+1][k]];
B:=B{[2..Length(B)]};
Apply(B,x->DeformCell(n-1,x));
B:=Concatenation(B);
Apply(B,i->bij[n][i]);

for i in B do
Add(b,[i,one]);
#b[i]:=b[i]+1;
od;

b:=SortedList(b);
for i in [1..Length(b)-1] do
if b[i][1]=b[i+1][1] then
b[i+1][2]:=b[i+1][2]+b[i][2];
b[i][2]:=0;
fi;
od;
b:=Filtered(b,x->not IsZero(x[2]));

return b;
end;
######################
else
characteristic:=0;
######################
Boundary:=function(n,k)
local b,i,j,B,sn;

b:=[];
B:=Y!.boundaries[n+1][basis[n+1][k]];
B:=B{[2..Length(B)]};
sn:=Y!.orientation[n+1][basis[n+1][k]];
B:=List([1..Length(B)],i->sn[i]*B[i]);

Apply(B,x->DeformCellSgn(n-1,x));
B:=Concatenation(B);
Apply(B,i->SignInt(i)*bij[n][AbsInt(i)]);



for i in B do
#b[AbsInt(i)]:=b[AbsInt(i)]+SignInt(i);
Add(b,[AbsInt(i),SignInt(i)]);
od;

b:=SortedList(b);
for i in [1..Length(b)-1] do
if b[i][1]=b[i+1][1] then 
b[i+1][2]:=b[i+1][2]+b[i][2];
b[i][2]:=0;
fi;
od;
b:=Filtered(b,x->not IsZero(x[2]));

return b;
end;
######################
fi;

if IsBound(Y!.orientation) then DeformCell:=DeformCellSgn; fi;
return
Objectify(HapSparseChainComplex,
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


#######################################################
#######################################################
InstallGlobalFunction(HAPRegularCWComplex,
function(arg)
local bnd, Y, dim, Coboundaries, k, i, j, n, b;

bnd:=arg[1];
Y:=Objectify(HapRegularCWComplex, rec());
dim:=PositionProperty(bnd,IsEmpty)-2;
Y!.properties:=[["dimension",dim]];


### COBOUNDARIES BEGIN ######################
Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
for n in [0..dim] do
  #k:=2*(n+1)+1;#k:=1+2^(n+1);
  Coboundaries[n+1]:=List(bnd[n+1],i->[0]);
  for j in [1..Length(bnd[n+2])] do
    b:=bnd[n+2][j];
    k:=Length(b);
    for i in b{[2..k]} do
      Coboundaries[n+1][i][1]:=Coboundaries[n+1][i][1]+1;
      Add(Coboundaries[n+1][i],j);
    od;
  od;
od;
Coboundaries[dim+1]:=List(bnd[dim+1],a->[0]);
### COBOUNDARIES END ###############################


Y!.boundaries:=bnd;
Y!.coboundaries:=Coboundaries;
Y!.vectorField:=fail;
Y!.inverseVectorField:=fail;
Y!.criticalCells:=fail;
if Length(arg)=2 then Y!.orientation:=arg[2]; 
else OrientRegularCWComplex(Y); fi;

####################
Y!.nrCells:=function(n);
if n>dim then return 0; fi;
return Length(Filtered(Y!.boundaries[n+1],x->not x[1]=0));
end;
####################

return Y;
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(ContractedRegularCWComplex,
function(arg)
local W,dd,Y, V,perm, d, d1, n, x, i, b, cnt, bnd,  dim , F, bool, orien;
W:=arg[1];

if Length(arg)>1 then dd:=arg[2];  ##Added 12/03/2018
else dd:=0; fi;

if IsBound(W!.orientation) then
Y:=RegularCWComplex(W!.boundaries,W!.orientation);
else
Y:=RegularCWComplex(W!.boundaries);
fi;

if IsBound(Y!.orientation) then bool:=true;
orien:=StructuralCopy(Y!.orientation);
else bool:=false; fi;

HAPContractRegularCWComplex(Y,dd);

bnd:=StructuralCopy(Y!.bnd);

perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;

for d in [1..Length(bnd)] do
F:=Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
#bnd[d]:=Filtered(bnd[d],x->not x[1]=0);
bnd[d]:=bnd[d]{F};
if bool and IsBound(orien[d]) then orien[d]:=orien[d]{F}; fi;
if d>1 then
d1:=d-1;
for x in bnd[d] do
for i in [2..Length(x)] do
x[i]:=perm[d1][x[i]];
od;
od;
fi;
od;

if bool then V:=RegularCWComplex(bnd, orien);
else V:= RegularCWComplex(bnd);fi;

V!.perm:=perm;

return V;
end);
#######################################################
#######################################################

################################################
################################################
InstallGlobalFunction(VerticesOfRegularCWCell,
function(Y,n,k)
local V, U, N, tmp, v ;

if n=0 then return [k]; fi;

N:=n+1;
V:=StructuralCopy(Y!.boundaries[N][k]);
V:=V{[2..Length(V)]};
V:=SSortedList(V);
N:=N-1;

while N>1 do
tmp:=[];
for v in V do
U:=StructuralCopy(Y!.boundaries[N][v]);
U:=U{[2..Length(U)]};
Append(tmp,U);
od;
tmp:=SSortedList(tmp);
V:=tmp;
N:=N-1;
od;

return V;
end);
################################################
################################################

################################################
################################################
InstallGlobalFunction(BoundaryOfRegularCWCell,
function(Y,n,k)
local V, U, N, tmp, v , cells;

if n=0 then return []; fi;

N:=n+1;
V:=StructuralCopy(Y!.boundaries[N][k]);
V:=V{[2..Length(V)]};
V:=SSortedList(V);
cells:=List(V,i->[N-2,i]);
N:=N-1;

while N>1 do
tmp:=[];
for v in V do
U:=StructuralCopy(Y!.boundaries[N][v]);
U:=U{[2..Length(U)]};
Append(tmp,U);
Append(cells,List(U,i->[N-2,i]));   #Changed N-1 to N-2  (May 2022)
od;
tmp:=SSortedList(tmp);
cells:=SSortedList(cells);
V:=tmp;
N:=N-1;
od;

return cells;
end);
################################################
################################################



#######################################################
#######################################################
InstallGlobalFunction(SimplifiedRegularCWComplex,
function(Y)
local W , a, b, OnceSimplifiedRegularCWComplex;


#######################################################
#######################################################
OnceSimplifiedRegularCWComplex:=function(W)
local Y, perm, cnt, JoinCells, d, d1, n, x, i, b,  cobnd, bnd,  dim , F, bool, orien,  pos ;

if IsBound(W!.orientation) then
Y:=RegularCWComplex(1*W!.boundaries,1*W!.orientation);
else
Y:=RegularCWComplex(1*W!.boundaries);

fi;

if IsBound(Y!.orientation) then bool:=true;
orien:=1*Y!.orientation;

else bool:=false; fi;

bnd:=Y!.boundaries;
cobnd:=Y!.coboundaries;


###################################################
###################################################
JoinCells:=function(d1,n)
                                 #The n-th cell in dimension d=d1-1 is removed
                                 #assuming it has a coboundary of size 2.
local V1,V2,V3,cob, d, a, b, d2,d3, m, s, t, pos, poss,u;

##
##CHECK IF REMOVAL SHOULD TAKE PLACE
if not (cobnd[d1][n][1] =2 and bnd[d1][n][1]>0) then return false; fi;
d2:=d1+1;
d3:=d2+1;
cob:=1*cobnd[d1][n];
if not SortedList(cobnd[d2][cob[2]])= SortedList(cobnd[d2][cob[3]]) then return false; fi;
V1:=BoundaryOfRegularCWCell(Y,d1-1,n);
V2:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][2]);
V3:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][3]);
if not 
 1+Size(V1)=Size(Intersection(V2,V3))  
then return false; fi;
##
##CHECK DONE

##
##REDUCE COBOUNDARIES OF BOUNDARIES OF nTH CELL
if d1>1 then
d:=d1-1;
for m in bnd[d1][n]{[2..Length(bnd[d1][n])]} do
t:=1*cobnd[d][m]{[2..Length(cobnd[d][m])]};
poss:=Position(t,n);
Remove(t,poss);
cobnd[d][m]:=Concatenation([Length(t)],t);
od;
fi;
##
##COBOUNDARIES OF BOUNDARIES REDUCED

##
##REMOVE nTH CELL, ITS COBOUNDARY, AND ADJUST ITS PRESENCE IN BOUNDARIES
##OF ITS COBOUNDARIES
bnd[d1][n]:=[0];
cobnd[d1][n]:=[0];

s:=bnd[d2][cob[2]];
s:=s{[2..Length(s)]};
pos:=Position(s,n);
Remove(s,pos);

if bool then a:=orien[d2][cob[2]][pos]; Remove(orien[d2][cob[2]],pos); fi;

t:=bnd[d2][cob[3]];
t:=t{[2..Length(t)]};
pos:=Position(t,n);
Remove(t,pos);

if bool then b:=orien[d2][cob[3]][pos]; Remove(orien[d2][cob[3]],pos); fi;

u:=SortedList(Concatenation(s,t));

bnd[d2][cob[2]]:=Concatenation([Length(u)],u);

if bool then orien[d2][cob[2]]:=Concatenation(orien[d2][cob[2]],-a*b*orien[d2][cob[3]]); fi;
##
##nTH CELL AND ITS PRESENCE REMOVED


##
##FOR SECOND CELL OF DIMENSION n+1  REDUCE THE BOUNDARIES OF ITS 
##COBOUNDARIES

for m in cobnd[d2][cob[3]]{[2..Length(cobnd[d2][cob[3]])]} do
t:=bnd[d3][m]{[2..Length(bnd[d3][m])]};

pos:=Position(t,cob[3]);
Remove(t,pos);
bnd[d3][m]:=Concatenation([Length(t)],t);
if bool then Remove(orien[d3][m],pos); fi;
od;


##
##SECOND CELL  BOUNDARIES OF ITS COBOUNDARIES REDUCED

##
##REMOVE PRESENCE OF SECOND CELL IN COBOUNDARIES OF ITS BOUNDARIES
for m in bnd[d2][cob[3]]{[2..Length(bnd[d2][cob[3]])]} do
if cobnd[d1][m][1]>0 then

t:=cobnd[d1][m]{[2..Length(cobnd[d1][m])]};

pos:=Position(t,cob[3]);
t[pos]:=cob[2];
t:=SSortedList(t);
cobnd[d1][m]:=Concatenation([Length(t)],t);

fi;
od;


bnd[d2][cob[3]]:=[0];

cobnd[d2][cob[3]]:=[0];
return true;
end;


###################################################
###################################################



###################################################
######SIMPLIFICATION STARTS########################

for d in [0..Dimension(Y)-1] do
d1:=d+1;

for n in [1..Length(bnd[d1])] do
if cobnd[d1][n][1] =2 and bnd[d1][n][1]>0 then 
JoinCells(d1,n);   
fi;
od;

od;

######SIMPLIFICATION DONE##########################
###################################################

perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;

for d in [1..Length(bnd)] do
F:=1*Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
bnd[d]:=1*((bnd[d]){F});
if bool and IsBound(orien[d]) then orien[d]:=orien[d]{F}; fi;

if d>1 then
d1:=d-1;
for x in bnd[d] do
for i in [2..Length(x)] do
x[i]:=1*perm[d1][x[i]];

od;
od;
fi;
od;

#if bool then return RegularCWComplex(bnd, orien);
if false then return RegularCWComplex(bnd, orien);
else
return RegularCWComplex(bnd);
fi;
end;
#######################################################
#######################################################


W:=OnceSimplifiedRegularCWComplex(Y);


a:=Size(Y);
b:=Size(W);

while a>b do
W:=OnceSimplifiedRegularCWComplex(W);
a:=b;
b:=Size(W);
od;

return W;
end);
#######################################################
#######################################################


#####################################################
#####################################################
InstallGlobalFunction(IsPureRegularCWComplex,
function(Y)
local n, x, dim, bool;

dim:=Dimension(Y);
bool:=true;

for n in [1..dim] do
if not bool then break; fi;
for x in Y!.coboundaries[n] do
if x[1]=0 then bool:=false; break; fi;
od;
od;

return bool;
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(BoundaryOfPureRegularCWComplex,
function(Y)
local F, n, m, d, t, i, x, bool, perm, cnt, d1,  dim, bnd, cobnd, orien, B;

if not IsPureRegularCWComplex(Y) then return fail; fi;


dim:=Dimension(Y);

bnd:=Y!.boundaries*1;
cobnd:=Y!.coboundaries*1;
if IsBound(Y!.orientation) then
orien:=Y!.orientation*1; bool:=true;
else bool:=false;
fi;

bnd[dim+1]:=[];

d:=dim;
  for n in [1..Length(bnd[d])] do
    if cobnd[d][n][1]>1 then
      if d>1 then
        t:=bnd[d][n];
        for m in t{[2..Length(t)]} do
          cobnd[d-1][m][1]:=cobnd[d-1][m][1]-1;
        od;
      fi;
      bnd[d][n]:=[0];
    fi;
  od;


for d in Reversed([1..dim-1]) do
  for n in [1..Length(bnd[d])] do
    if cobnd[d][n][1]=0 then
      if d>1 then
        t:=bnd[d][n];
        for m in t{[2..Length(t)]} do
          cobnd[d-1][m][1]:=cobnd[d-1][m][1]-1;
        od;
      fi;
      bnd[d][n]:=[0];
    fi;
  od;
od;


perm:=[];  ##perm[d][i] will be the new position of
           ##the old i-th cell of dimension d;

for d in [1..Length(bnd)] do
perm[d]:=[];
cnt:=0;
for n in [1..Length(bnd[d])] do
if bnd[d][n][1]=0 then cnt:=cnt+1;
else
perm[d][n]:=n-cnt;
fi;
od;
od;

for d in [1..Length(bnd)] do
F:=Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
#bnd[d]:=Filtered(bnd[d],x->not x[1]=0);
bnd[d]:=bnd[d]{F};
if bool and IsBound(orien[d]) then orien[d]:=orien[d]{F}; fi;
if d>1 then
d1:=d-1;
for x in bnd[d] do
for i in [2..Length(x)] do
x[i]:=perm[d1][x[i]];
od;
od;
fi;
od;


if IsBound(orien) then
B:=RegularCWComplex(bnd,orien);
B!.perm:=perm;
return B;
else
B:=RegularCWComplex(bnd);
B!.perm:=perm;
return B;
fi;

end);
#####################################################
#####################################################

##################################################
##################################################
InstallGlobalFunction(OrientRegularCWComplex,
function(Y)
local bnd, cobnd, orien, dim, d, d1, d2, x, i, j, b, bb, sn, m, S, T, s, t, 
       L, bool ;

if IsBound(Y!.orientation) then
if not Y!.orientation=fail then
return ;
fi;fi;

dim:=Dimension(Y);
bnd:=Y!.boundaries;
cobnd:=Y!.coboundaries;
orien:=bnd*0;
for x in orien do
Apply(x,y->y{[2..Length(y)]});
od;

Apply(orien[1], x->x+1);
Apply(orien[2],x->[1,-1]);


#######################
for d in [3..dim+1] do
  d1:=d-1;
  d2:=d-2;
  for i in [1..Length(bnd[d])] do
    b:=bnd[d][i]{[2..Length(bnd[d][i])]};
    bb:=[];
    for j in [1..Length(b)] do
      Add(bb, bnd[d1][b[j]]{[2..Length(bnd[d1][b[j]])]}   );
    od;
    orien[d][i][1]:=1;
    S:=[1..Length(b)];
    T:=[1..Length(b)];
    for s in [2..Length(b)] do
    Unbind(S[s]);
    od;
    Unbind(T[1]);

    while 0 in orien[d][i] do
    ###############################
    bool:=false;
    for s in S do
    for t in T do
      L:=Intersection(bb[s], bb[t]);
      if Length(L)>0 then
        S[t]:=t;
        Unbind(T[t]);   
        bool:=true;
        if orien[d][i][s]*orien[d1][b[s]][Position(bb[s],L[1])]=
           orien[d1][b[t]][Position(bb[t],L[1])]
           then orien[d][i][t]:=-1;
           else
           orien[d][i][t]:=1;
        fi; 
        break;
      fi;
    od;
    od;
    ###############################
    od;
  od;
od;
#######################


Y!.orientation:=orien;
end);
##################################################
##################################################

#############################################
#############################################
InstallGlobalFunction(HAP_Sequence2Boundaries,
function(Y)

local orien, b, i, bb, newb, neworien, pos, s, t;

for i in [1..Y!.nrCells(2)] do
  b:=Y!.boundaries[3][i];
  orien:=Y!.orientation[3][i];
  bb:=List(b{[2..Length(b)]},  j->Y!.boundaries[2][j]{[2,3]});

  s:=bb[1];
  Unbind(bb[1]);
  newb:=[b[1],b[2]];
  neworien:=[orien[1]];
  while Length(newb)<Length(b) do
    for t in bb do
      if Length(Intersection(s,t))>0 then pos:=Position(bb,t); break; fi;
    od;
    s:=bb[pos];
    Unbind(bb[pos]);
    Add(newb,b[pos+1]);
    Add(neworien,orien[pos]);
  od;
  Y!.boundaries[3][i]:=newb;
  Y!.orientation[3][i]:=neworien;


od;
end);
#############################################
#############################################

############################################
InstallGlobalFunction(CubicalComplex,
function(A)
local dim, dims;

dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);


return Objectify(HapCubicalComplex,
           rec(
           binaryArray:=A,
           properties:=[
           ["dimension",dim],
           ["arraySize",dims]]
           ));

end);
############################################
############################################
############################################
InstallGlobalFunction(ReadImageAsWeightFunction,
function(file,f)
local MM,F, M, A, B, C, Y, k,k1, i,j,W, weight, coord, x1;

MM:=ReadImageAsFilteredPureCubicalComplex(file,3);
A:=PureCubicalComplexToCubicalComplex(MM);;
A:=0*A!.binaryArray;

for i in [1..f] do
F:=FiltrationTerm(MM,i);
F:=PureCubicalComplexToCubicalComplex(F);
A:=A+F!.binaryArray;
od;


B:=A*0;;
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
if A[i][j]>0 then B[i][j]:=1; fi;
od;od;

M:=CubicalComplex(B);
C:=ChainComplex(M);

W:=[];
for k in [1..1+Dimension(M)] do
W[k]:=[];k1:=k-1;
for i in [1..C!.dimension(k1)] do
coord:=C!.positionToCoordinate[k][i];
W[k][i]:=A[coord[2]][coord[1]];
od;
od;

Unbind(C);
Y:=CubicalComplexToRegularCWComplex(M);
Unbind(M);

#####################
weight:=function(k,i);
return W[k+1][i];
end;
#####################


return [Y, weight];

end);
#######################################################
#######################################################


######################################################
#######################################################
InstallGlobalFunction(EulerIntegral,
function(Y,weight)
local k, i, eulint, alpha, sn;

eulint:=0;
for k in [0..Dimension(Y)] do
sn:=(-1)^k;
alpha:=0;
for i in [1..Y!.nrCells(k)] do
alpha:=alpha+weight(k,i);
od;
eulint:=eulint+sn*alpha;
od;

return eulint;

end);
######################################################
######################################################


############################################
############################################
InstallGlobalFunction(HAP_SimplicialPairToCWMap,
function(T,S)
local TT,SS, VT, V, SimpT, SimpS, i, fn, F;

VT:=T!.vertices;
###################
fn:=function(i);
#return i;
return Position(VT,i);
end;
###################

SimpT:=[];
for i in [1..Length(T!.simplicesLst)] do
SimpT[i]:=List(T!.simplicesLst[i],x->List(x,fn));
od;

VT:=S!.vertices;
###################
fn:=function(i);
#return i;
return Position(VT,i);
end;
###################

SimpS:=[];
for i in [1..Length(S!.simplicesLst)] do
SimpS[i]:=List(S!.simplicesLst[i],x->List(x,fn));
od;

TT:=SimplicialComplex(Concatenation(SimpT));
SS:=SimplicialComplex(Concatenation(SimpS));

TT:=RegularCWComplex(TT);

SS:=RegularCWComplex(SS);

################
################
F:=function(n,k);
return Position(T!.simplicesLst[n+1],S!.simplicesLst[n+1][k]);
end;
################
################


return Objectify(HapRegularCWMap,
       rec(
           source:=SS,
           target:=TT,
           mapping:=F,
           SimpT:=SimpT,
           SimpS:=SimpS));

end);
############################################
############################################


#####################################
#####################################
InstallGlobalFunction(HAP_PureCubicalPairToCWMap,
function(M,A)
local F, YA, YM;

YA:=RegularCWComplex(A);
YM:=RegularCWComplex(M);

################
################
F:=function(n,k)
local x,i, v;
x:=YA!.positionToCoordinate[n+1][k];
v:=1*YM!.coordinateToPosition;
for i in Reversed(x) do
v:=v[i];
od;
return v;
end;
################
################

return Objectify(HapRegularCWMap,
       rec(
       source:=YA,
       target:=YM,
       mapping:=F
       ));
end);
######################################
######################################

######################################
######################################
InstallGlobalFunction(GraphOfRegularCWComplex,
function(Y)
local A, n, x;

n:=Y!.nrCells(0);
A:=NullMat(n,n);

for x in Y!.boundaries[2] do
A[x[2]][x[3]]:=1;
A[x[3]][x[2]]:=1;
od;

return IncidenceMatrixToGraph(A);

end);
######################################
######################################

######################################
######################################
InstallGlobalFunction(HomotopyGraph,
function(W)
local Y, V, A, n, x, i;

if Dimension(W)=1 then
return ContractedComplex(Graph(W));
fi;

Y:=1*W!.boundaries{[1..3]};
Add(Y,[]);
Y:=RegularCWComplex(Y);
CocriticalCellsOfRegularCWComplex(Y,2);
V:=SSortedList(Flat(Y!.vectorField[2]));

n:=Y!.nrCells(0);
A:=NullMat(n,n);

#for x in Y!.boundaries[2] do
for i in [1..Y!.nrCells(1)] do
if not i in V then
x:=Y!.boundaries[2][i];
A[x[2]][x[3]]:=1;
A[x[3]][x[2]]:=1;
fi;
od;

return ContractedComplex(IncidenceMatrixToGraph(A));

end);
######################################
######################################


######################################
######################################
InstallGlobalFunction(DeformationRetract,
function(Y)
local R,map, perm, invperm, P, IP, i;

R:=ContractedComplex(Y);

perm:=R!.perm;
invperm:=[];

for P in perm do
IP:=[];
for i in [1..Length(P)] do
if IsBound(P[i]) then
IP[P[i]]:=i;
fi;
od;
Add(invperm,IP);
od;


##################
map:=function(n,k);
return invperm[n+1][k];
end;
##################

return Objectify(HapRegularCWMap,
       rec(
           source:=R,
           target:=Y,
           mapping:=map));

end);
######################################
######################################

###############################
###############################
InstallGlobalFunction(NonRegularCWBoundary,
function(Y,n,k)
local DeformCell , bnd;

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

bnd:=Y!.boundaries[n+1][k];
bnd:=bnd{[2..Length(bnd)]};
Apply(bnd,i->DeformCell(n-1,i));
bnd:=Concatenation(bnd);
bnd:=SortedList(bnd);
return bnd;

end);
###############################
###############################



#######################################################
#######################################################
InstallGlobalFunction(HAP_BaryCentricSubdivisionRegularCWComplex,
function(Y)
local  B, chains, newchains, maxchains, shifts, s, n, b, c, i;

B:=1*Y!.coboundaries;
for n in [1..Length(B)] do
for b in B[n] do
Remove(b,1);
od;
od;

shifts:=[];
s:=0;
for n in [0..Dimension(Y)] do
shifts[n+1]:=s;
s:=s+Y!.nrCells(n);;
od;

maxchains:=[];
chains:=List([1..Y!.nrCells(0)],i->[i]);
newchains:=[];

for n in [1..Dimension(Y)+1] do
for c in chains do
   if Length(B[n][c[Length(c)]])=0 then Add(maxchains,c);
else
   for i in B[n][c[Length(c)]] do
   Add(newchains, Concatenation(c,[i]));
   od;
fi;
od;
chains:=newchains;
newchains:=[];
od;

for c in maxchains do
for n in [1..Length(c)] do
c[n]:=c[n]+shifts[n];
od;
od;

return MaximalSimplicesToSimplicialComplex(maxchains);
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallOtherMethod(BarycentricSubdivision,
"for regular CW complexes",
[IsHapRegularCWComplex],
function(Y);
return HAP_BaryCentricSubdivisionRegularCWComplex(Y);
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(HAP_Triangulation, 
function(Y)
local K, B, chains, newchains, maxchains, reduce, vsets, s, n, b, c, i;

B:=1*Y!.coboundaries;
for n in [1..Length(B)] do
for b in B[n] do
Remove(b,1);
od;
od;

maxchains:=[];
chains:=List([1..Y!.nrCells(0)],i->[i]);
newchains:=[];

###################################
reduce:=function(x);
if Length(x)>1 then return [x[1],x[Length(x)]];
else return x;
fi;
end;
###################################

for n in [1..Dimension(Y)+1] do
for c in chains do
   if Length(B[n][c[Length(c)]])=0 then Add(maxchains,c);
else
   for i in B[n][c[Length(c)]] do
   Add(newchains, Concatenation(c,[i]));
   od;
fi;
od;
Apply(newchains,reduce);
newchains:=SSortedList(newchains);
chains:=newchains;
newchains:=[];
od;

Apply(maxchains, x->[x[1],x[Length(x)]]);

vsets:=List([1..Y!.nrCells(0)],a->[]);
for c in maxchains do
AddSet(vsets[c[1]],c[2]);
od;
vsets:=SSortedList(vsets);
return MaximalSimplicesToSimplicialComplex(vsets);
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallOtherMethod(Nerve,
"For pure Cw complexes",
[IsHapRegularCWComplex],
function(Y)
return HAP_Triangulation(Y);
end);
#######################################################
#######################################################

#######################################################
#######################################################
InstallGlobalFunction(ComposeCWMaps,
function(F,G)
local X,Y,map;       #THIS FUNTION IS NOT NEEDED

#  FG: X ---G---> W  ---F---> Y
X:=Source(G);
Y:=Target(F);
#######################
map:=function(n,k);
return F!.mapping(n,G!.mapping(n,k));
end;
#######################

return Objectify(HapRegularCWMap,
       rec(
           source:=X,
           target:=Y,
           mapping:=map  ));
end);
#######################################################
#######################################################

###############################################################
###############################################################
InstallGlobalFunction(SimplicialComplexToRegularCWComplex,
function(KK)
local K, posn, popP, P, S, B, n, b, s, i, bb, N;

K:=IntegerSimplicialComplex(KK);
S:=K!.simplicesLst;
P:=[];
###################################
popP:=function()
local cnt,s;
if n=1 then P[n]:=SSortedList(S[n]); P[n]:=List(P[n],x->x[1]);fi;

if n=2 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do
cnt:=cnt+1;
if IsInt(P[n][s[1]]) then P[n][s[1]]:=[]; fi;
P[n][s[1]][s[2]]:= cnt;
od;
fi;

if n=3 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do
cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]]) then P[n][s[1]][s[2]]:=[]; fi;
P[n][s[1]][s[2]][s[3]]:= cnt;
od;
fi;

if n=4 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do
cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]]) then P[n][s[1]][s[2]][s[3]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]]:= cnt;
od;
fi;

if n=5 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do
cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]]) then P[n][s[1]][s[2]][s[3]][s[4]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]]:= cnt;
od;
fi;

if n=6 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]]:= cnt;
od;
fi;

if n=7 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]]:= cnt;
od;
fi;

if n=8 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]]:= cnt;
od;
fi;

if n=9 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]]:= cnt;
od;
fi;

if n=10 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]]:= cnt;
od;
fi;

if n=11 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]]:= cnt;
od;
fi;

if n=12 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]]:= cnt;
od;
fi;

if n=13 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]]:= cnt;
od;
fi;

if n=14 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]]:= cnt;
od;
fi;

if n=15 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]]:= cnt;
od;
fi;

if n=16 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]][s[16]]:= cnt;
od;
fi;

if n=17 then P[n]:=P[n-1];cnt:=0;
for s in S[n] do cnt:=cnt+1;
if IsInt(P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]][s[16]]) then P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]][s[16]]:=[]; fi;
P[n][s[1]][s[2]][s[3]][s[4]][s[5]][s[6]][s[7]][s[8]][s[9]][s[10]][s[11]][s[12]][s[13]][s[14]][s[15]][s[16]][s[17]]:= cnt;
od;
fi;



if n>17 then P[n]:=SSortedList(S[n]); fi;
end;
###################################


###################################
posn:=function(bb)
local p;
if n=1 then return P[n][bb[1]]; fi;
if n=2 then return P[n][bb[1]][bb[2]]; fi;
if n=3 then return P[n][bb[1]][bb[2]][bb[3]]; fi;
if n=4 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]]; fi;
if n=5 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]]; fi;
if n=6 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]]; fi;
if n=7 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]]; fi;
if n=8 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]]; fi;
if n=9 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]]; fi;
if n=10 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]]; fi;
if n=11 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]]; fi;
if n=12 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]]; fi;
if n=13 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]][bb[13]]; fi;
if n=14 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]][bb[13]][bb[14]]; fi;
if n=15 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]][bb[13]][bb[14]][bb[15]]; fi;
if n=16 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]][bb[13]][bb[14]][bb[15]][bb[16]]; fi;
if n=17 then return P[n][bb[1]][bb[2]][bb[3]][bb[4]][bb[5]][bb[6]][bb[7]][bb[8]][bb[9]][bb[10]][bb[11]][bb[12]][bb[13]][bb[14]][bb[15]][bb[16]][bb[17]]; fi;

p:=Position(P[n],bb);
return p;
end;
###################################

B:=[List(S[1],x->[1,0])];
for n in [1..Length(S)-1] do
   popP();
   B[n+1]:=[];
   for s in S[n+1] do
      b:=[];
      for i in [1..n+1] do
         bb:=Concatenation(s{[1..i-1]},s{[i+1..n+1]});

         Add(b, posn(bb));

      od;
   b:=Concatenation( [Length(b)], b);
   Add(B[n+1],b);
   od;
od;

return RegularCWComplex(B);
end);
###############################################################
###############################################################

###################################################
###################################################
InstallGlobalFunction("ClassifyingSpaceFiniteGroup",
function(G,n)
local R,Y,T,W,i;

R:=ResolutionFiniteGroup(G,n);
T:=Group(One(G));
R!.stabilizer:=function(n,i); return T; end;
R:=HAP_BaryCentricSubdivisionGComplex(R);
R:=ResolutionToEquivariantCWComplex(R);
Y:=EquivariantCWComplexToRegularCWComplex(R,G);

#In general we can't be sure that R corresponds to a *regular* CW-complex.
#The paper "computing group resolutions" explains that R corresponds to a
##CW-complex X. This X can certainly be represented, up to homotopy,
#by the data type for a HAPRegularCWComplex but (despite the false claim in 
#the book "An invitation to computational Homotopy") we need to check if it is
#a regular CW-complex. In the case of cyclic groups with standard minimal 
#resolution it clearly is.

if IsCyclic(G) then
return SimplifiedComplex(Y);
fi;

#Let's check to see if X is regular
#If X fails the test then  we could add more cells to X. But this is not yet implemented.
W:=EquivariantCWComplexToRegularCWComplex(R,Group(One(G)));
CriticalCells(W);
for i in [1..n-1] do
if not Homology(W,i)=[] then
Print("This function is currently unable to construct a regular CW classyfying space for that group. You could try entering an isomorphic copy of the group whose elements are listed in a different order.\n");
return fail; fi;
od;
return SimplifiedComplex(Y);
end);
###################################################
###################################################

##############################################################
##############################################################
InstallGlobalFunction(RegularCWComplexReordered,
function(Y)
local n, k, b, m, sigma, sigmai, bnds, newbnds;

bnds:=1*Y!.boundaries;
for n in [0..Dimension(Y)] do
sigma:=Random(SymmetricGroup(Y!.nrCells(n)));
sigmai:=sigma^-1;
newbnds:=List([1..Y!.nrCells(n)],i->bnds[n+1][i^sigma]);;
bnds[n+1]:=newbnds;
newbnds:=[];
for k in [1..Y!.nrCells(n+1)] do
b:=bnds[n+2][k];
m:=b[1];
b:=b{[2..b[1]+1]};
Apply(b,i->i^sigmai);
b:=SortedList(b);
b:=Concatenation([m],b);
Add(newbnds,b);
od;
bnds[n+2]:=newbnds;
od;

return RegularCWComplex(bnds);;
end);
##############################################################
##############################################################
