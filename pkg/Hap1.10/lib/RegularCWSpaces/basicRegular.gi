#C Graham Ellis

##################################################################
##################################################################
InstallGlobalFunction(SimplicialComplexToRegularCWSpace,
function(K)
local
	NrCells,Boundaries,tmp,TMP,Coboundaries,VectorField,Properties,
        cnt,b,bb,k,n,s,x,i,j,dim;

####################
NrCells:=function(n);
if n>dim then return 0; fi;
return Length(Filtered(Boundaries[n+1],x->not x[1]=0));
end;
####################

VectorField:=fail;
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

 Boundaries[n+1][k]:=Concatenation([Length(b)],b,List(b,j->1));
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
  for b in Coboundaries[n+1] do
  Append(b,List([1..Length(b)-1],a->1)); 
  od;
od;
Coboundaries[dim+1]:=List(Boundaries[dim+1],a->[0]);
### COBOUNDARIES END ###############################


return Objectify(HapRegularCWSpace,
       rec(
           nrCells:=NrCells,
           boundaries:=Boundaries,
           coboundaries:=Coboundaries,
           vectorField:=VectorField,
           properties:=Properties));

end);
##################################################################
##################################################################


#############################################
#############################################
InstallGlobalFunction(ContractRegularCWSpace,
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

MCoboundaries:=Y!.coboundaries[n+1];
MBoundaries:=Y!.boundaries[n+1];
UCoboundaries:=Y!.coboundaries[n+2];
UBoundaries:=Y!.boundaries[n+2];
if n>0 then
  LCoboundaries:=Y!.coboundaries[n];
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
      ###
  if n>0 then
    b:=MBoundaries[i];
    for j in StructuralCopy(b{[2..1+b[1]]}) do
     t:=LCoboundaries[j][1];
     LCoboundaries[j][1]:=LCoboundaries[j][1]-1;
     cob:=LCoboundaries[j];
     pos:=Position(cob{[2..t+1]},i);
     LCoboundaries[j]:=Concatenation(cob{[1..pos]},cob{[2+pos..t+pos]},
                                          cob{[t+pos+2..Length(cob)]});
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
     MCoboundaries[j]:=Concatenation(cob{[1..pos]},cob{[2+pos..t+pos]},
                                          cob{[t+pos+2..Length(cob)]});
    od;
      ###
  MBoundaries[i]:=[0];
  UBoundaries[U]:=[0];
  UCoboundaries[U]:=[0];
  MCoboundaries[i]:=[0];
fi;
od;

Y!.boundaries[n+2]:=UBoundaries;
Y!.coboundaries[n+2]:=UCoboundaries;
Y!.boundaries[n+1]:=MBoundaries;
Y!.coboundaries[n+1]:=MCoboundaries;
if n>0 then
  Y!.coboundaries[n]:=LCoboundaries;
fi;

Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
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
InstallGlobalFunction(RemoveCellFromRegularCWSpace,
function(Y,dim,n)
local  bnd, x,tmp;

#Remove the n-th cell in dimension dim

dim:=dim+1;

if dim=1 then Y!.boundaries[1][n]:=[0]; return [dim-1,n]; fi;

bnd:=Y!.boundaries[dim][n];
bnd:=bnd{[2..(Length(bnd)+1)/2]};
Y!.boundaries[dim][n]:=[0];

if dim>1 then
for x in bnd do
tmp:=Y!.coboundaries[dim-1][x];
tmp[1]:=tmp[1]-1;
tmp[Length(tmp)]:=-42;
tmp[Position(tmp{[2..Length(tmp)]},n)+1]:=-42;

tmp:=Filtered(tmp,i->not i = -42);

Y!.coboundaries[dim-1][x]:=tmp;
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
if Length(arg)>1 then
   ContractSpace:=ContractRegularCWSpace_Alt;
else
   ContractSpace:=ContractRegularCWSpace;
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

  if Size(Y)=0  then  return cells; fi;

  pos:=0;

  while true do
    pos:=pos+1;
    ppos:=PositionProperty(Y!.boundaries[dim+1]{[pos..Length(Y!.boundaries[dim+1])]},                          x->x[1]>0);    if ppos=fail then dim:=dim-1; break; fi;

  pos:=pos+ppos-1;

    c:=RemoveCellFromRegularCWSpace(Y,dim,pos);

    Add(cells,c);

    ContractSpace(Y);

  od;
od;

return cells;
end);
##########################################################
##########################################################



