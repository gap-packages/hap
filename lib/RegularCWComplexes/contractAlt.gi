#C Graham Ellis

#############################################
#############################################
InstallGlobalFunction(HAPContractRegularCWComplex_Alt,
function(Y)
local
      Contract, ii, nn, dim, bool, BOOL;


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
if not IsBound(Y!.free) then Y!.free:=[]; fi;
if not IsBound(Y!.free[n+1]) then 
Y!.free[n+1]:=Filtered([1..C],i->MCoboundaries[i][1]=1);
fi;
#Free:=Filtered([1..C],i->MCoboundaries[i][1]=1);
#if Length(Free)=0 then return false;fi;
if Length(Y!.free[n+1])=0 then return false;fi;

#######################
#######################



#for i in Free do
for i in Y!.free[n+1] do
RemoveSet(Y!.free[n+1],i);

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
#if t=2 then Add(Free,j);fi;############################ADDED
if t=2 then AddSet(Y!.free[n+1],j);fi;
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

#if Length(Free)>0 then return true;
if Length(Y!.free[n+1])>0 then return true;
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





