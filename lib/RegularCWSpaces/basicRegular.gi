#C Graham Ellis

##################################################################
##################################################################
InstallGlobalFunction(SimplicialComplexToRegularCWSpace,
function(K)
local
	NrCells,Boundaries,tmp,TMP,Coboundaries,VectorField,Properties,
        cnt,b,bb,k,n,s,x,i,dim;

####################
NrCells:=function(n);
if n>dim then return 0; fi;
return Length(List(Boundaries[n+1],x->not x[1]=0));
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

 b:=List(bb,x-> Difference(bb,[x])    );
 b:=List(b,x->   TMP[x[1]][Position(tmp[x[1]],x)] );

 ####
 #The above should have the same effect as
 #b:=List(bb,x->K!.enumeratedSimplex(Difference(bb,[x]))); 
 ####

 Boundaries[n+1][k]:=Concatenation([Length(b)],b,List(b,j->1));
 od;
od;
Boundaries[dim+2]:=[];
### BOUNDARIES END ###############################

### COBOUNDARIES BEGIN ######################
Coboundaries:=[];; #Coboundaries[n+1] contains the info on n-cells.
for n in [0..dim] do

  Coboundaries[n+1]:=[];
  tmp:=List(Boundaries[1],x->[]);
  TMP:=List(Boundaries[1],x->[]);
  cnt:=0;
  for s in K!.simplicesLst[n+2] do
    cnt:=cnt+1;
    for i in [1..Length(s)] do
       Add(tmp[s[i]],s);
       Add(TMP[s[i]],cnt);
    od; 
  od;
 for k in [1..K!.nrSimplices(n)] do
 bb:=K!.simplices(n,k);
 b:=List(bb,i->tmp[i]);
 b:=Intersection(b);
 b:=List(b,x->TMP[x[1]][Position(tmp[x[1]],x)]);
 Coboundaries[n+1][k]:=Concatenation([Length(b)],b,List(b,j->1));
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
Free:=Filtered([1..C],i->MCoboundaries[i][1]=1);

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



