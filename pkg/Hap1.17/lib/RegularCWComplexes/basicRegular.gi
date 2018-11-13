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
InstallGlobalFunction(SimplicialComplexToRegularCWComplex,
function(arg)
local
	K,DM,NrCells,Boundaries,tmp,TMP,Coboundaries,Properties,
        Orientation, 
        cnt,b,bb,k,n,s,x,i,j,dim ;

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
if MCoboundaries[i][1]=1 then Add(Free,i);fi;
od;

#Print([Length(FREE),Length(Free)],"  ");

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
        Y,basis, bool, bij,Dimension, Boundary, one, zero, b, n, dim, characteristic, DeformCell,DeformCellSgn, DeformCellSgnHtpy, HomotopicalDeformCell, 
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
local sgnn,x,f,k,sgnk,cnt,bnd,def,sn,tog,def1,def2;
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
sn:=Y!.homotopyOrientation[n+2][f];  ##

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
Apply(def,x->[x,0]);

def:=HtpyRed(def);


Apply(def,x->x[1]);

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
cells:=List(V,i->[N-1,i]);
N:=N-1;

while N>1 do
tmp:=[];
for v in V do
U:=StructuralCopy(Y!.boundaries[N][v]);
U:=U{[2..Length(U)]};
Append(tmp,U);
Append(cells,List(U,i->[N-1,i]));
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
local Y, perm, cnt, JoinCells, d, d1, n, x, i, b,  cobnd, bnd,  dim , F, bool, orien,  pos;

if IsBound(W!.orientation) then
Y:=RegularCWComplex(StructuralCopy(W!.boundaries),StructuralCopy(W!.orientation));
else
Y:=RegularCWComplex(StructuralCopy(W!.boundaries));
fi;

if IsBound(Y!.orientation) then bool:=true;
orien:=StructuralCopy(Y!.orientation);
else bool:=false; fi;

bnd:=Y!.boundaries;
cobnd:=Y!.coboundaries;


###################################################
###################################################
JoinCells:=function(d1,n)
                                 #The n-th cell in dimension d=d1-1 is removed
                                 #assuming it has a coboundary of size 2.
local V1,V2,V3,cob, d, a, b, d2,d3, m, s, t, pos, poss ;

##
##CHECK IF REMOVAL SHOULD TAKE PLACE

V1:=BoundaryOfRegularCWCell(Y,d1-1,n);
V2:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][2]);
V3:=BoundaryOfRegularCWCell(Y,d1,cobnd[d1][n][3]);
#if not Size(V1)=Size(Intersection(V1,V2,V3)) then Print([d1-1,n],V1,V2,V3,"\n\n");
# fi;
if 
not (
cobnd[d1][n][1] =2 
and bnd[d1][n][1]>0  
and 1+Size(V1)=Size(Intersection(V2,V3)) )
then return false; fi;
##
##CHECK DONE

d2:=d1+1;
d3:=d2+1;
cob:=StructuralCopy(cobnd[d1][n]);
if not SortedList(cobnd[d2][cob[2]])= SortedList(cobnd[d2][cob[3]]) then return false; fi;

##
##REDUCE COBOUNDARIES OF BOUNDARIES OF nTH CELL
if d1>1 then
d:=d1-1;
for m in bnd[d1][n]{[2..Length(bnd[d1][n])]} do
t:=cobnd[d][m]{[2..Length(cobnd[d][m])]};
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

bnd[d2][cob[2]]:=Concatenation([Length(s)+Length(t)],s,t);

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
F:=Filtered([1..Length(bnd[d])],i->not bnd[d][i][1]=0);
bnd[d]:=bnd[d]{F};
if bool and IsBound(orien[d]) then orien[d]:=orien[d]{F}; fi;

if d>1 then
d1:=d-1;
for x in bnd[d] do
for i in [2..Length(x)] do
x[i]:=StructuralCopy(perm[d1][x[i]]);
od;
od;
fi;
od;

if bool then return RegularCWComplex(bnd, orien);
else
return RegularCWComplex(bnd);
fi;
end;
#######################################################
#######################################################


W:=OnceSimplifiedRegularCWComplex(Y);
#return W;

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

MM:=ReadImageAsFilteredPureCubicalComplex(file,f);
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
