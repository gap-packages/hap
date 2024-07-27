#(C) Graham Ellis

#############################################
#############################################
InstallGlobalFunction(HAPCocontractRegularCWComplex,
function(Y)
local
      Contract, nn, dim, bool, BOOL;

#############################################
##### The work-horse function.###############
Contract:=function(n)
local

      b, C, i, j, t, cob, pos, bool,
      Free, UBoundaries, UCoboundaries,
      MBoundaries, MCoboundaries, LCoboundaries, LBoundaries, U;

#This function removes pairs of n- and (n-1)-cells if possible.
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
if n<dim then
UCoboundaries:=Y!.cobnd[n+2];
fi;
UBoundaries:=Y!.bnd[n+2];
if n>0 then
  LCoboundaries:=Y!.cobnd[n];
  LBoundaries:=Y!.bnd[n];
fi;
C:=Length(MBoundaries);

#######################
#######################THIS TAKES ALL THE TIME
Free:=Filtered([1..C],i->MBoundaries[i][1]=1);
if Length(Free)=0 then return false;fi; 
#######################
#######################


for i in Free do
if MBoundaries[i][1]=1 then
Y!.inverseVectorField[n][MBoundaries[i][2]]:=i;
Y!.vectorField[n][i]:=MBoundaries[i][2];
      ###

     j:=MBoundaries[i][2];
        if n>0 then
     for t in LCoboundaries[j]{[2..Length( LCoboundaries[j])]} do
     b:=MBoundaries[t];
     pos:=Position(b{[2..Length(b)]},j);
     b[pos+1]:=-42;
     b:=Filtered(b,a-> not a=-42);
     b[1]:=b[1]-1; 
     if b[1]=0 then b[1]:=-1; fi;
     if b[1]=1 then Add(Free,t); fi;
     MBoundaries[t]:=b; 
     od;
     MBoundaries[i]:=[0];
        fi;
        if n>0 then
     LBoundaries[j]:=[0];
     LCoboundaries[j]:=[0];
        fi;
      ###
        if n<dim then
    for t in MCoboundaries[i]{[2..Length(MCoboundaries[i])]} do 
    b:=UBoundaries[t];
     pos:=Position(b{[2..Length(b)]},i);
     b[pos+1]:=-42;
     b:=Filtered(b,a-> not a=-42);
     b[1]:=b[1]-1;  
     if b[1]=0 then b[1]:=-1; fi;
     UBoundaries[t]:=b;
    od;
        fi;
     MCoboundaries[i]:=[0];
      ###
fi;
od;

Y!.bnd[n+2]:=UBoundaries;
if n<dim then
Y!.cobnd[n+2]:=UCoboundaries;
fi;
Y!.bnd[n+1]:=MBoundaries;
Y!.cobnd[n+1]:=MCoboundaries;
if n>0 then
  Y!.cobnd[n]:=LCoboundaries;
fi;

Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") or k<0 then return 0; fi;
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

while BOOL  do
BOOL:=false;
  for nn in [1..dim] do
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

##########################################################
##########################################################
InstallGlobalFunction(CocriticalCellsOfRegularCWComplex,
function(arg)
local Y,SK,ContractSpace,cells,dim,c,cpos,pos,ppos;

Y:=arg[1];

if Dimension(Y)=0 then                        #
return List([1..Y!.nrCells(0)],i->[0,i]); fi; #ADDED APRIL 2024

SK:=1+Minimum(arg[2],EvaluateProperty(Y,"dimension")); 
if (not Y!.criticalCells=fail) and
EvaluateProperty(Y,"codim")>=SK-1 then
return Filtered(Y!.criticalCells,x->x[1]<SK);
fi;

   ContractSpace:=HAPCocontractRegularCWComplex;


#######
if EvaluateProperty(Y,"codim")=fail then dim:=0;
else dim:=EvaluateProperty(Y,"codim");
fi;
#while true do
#if not Y!.nrCells(dim)=0 then break; fi;
#if dim=EvaluateProperty(Y,"dimension") then break; fi;
#dim:=dim+1;
#od;
#######

if not Y!.criticalCells=fail then
cells:=Y!.criticalCells;
else
cells:=[];
fi;
ContractSpace(Y);

#dim:=0;
while true do

  if
#Sum(List( [1..Length(Y!.bnd)],i->Y!.nrCells(i-1)))=0
Sum(List( [1..SK],i->Y!.nrCells(i-1)))=0

then 

###################
#Y!.vectorField:=fail;
cpos:=PositionProperty(Y!.properties,x->x[1]="codim");
if cpos=fail then
Add(Y!.properties,["codim",dim]);
else
Y!.properties[cpos][2]:=dim;
fi;
Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
#Unbind(Y!.bnd);
#Unbind(Y!.cobnd);
####################

Y!.criticalCells:=cells; return Filtered(cells,x->x[1]<SK); fi;

  pos:=0;



  while true do
    pos:=pos+1;
    ppos:=PositionProperty(Y!.bnd[dim+1]{[pos..Length(Y!.bnd[dim+1])]},                          x->not x[1]=0);    if ppos=fail then dim:=dim+1;  break; fi;

  pos:=pos+ppos-1;

if dim>SK-1 then  

######################
#Y!.vectorField:=fail;
cpos:=PositionProperty(Y!.properties,x->x[1]="codim");
if cpos=fail then
Add(Y!.properties,["codim",dim]);
else
Y!.properties[cpos][2]:=dim;
fi;

Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
#Unbind(Y!.bnd);
#Unbind(Y!.cobnd);
######################

Y!.criticalCells:=cells; return Filtered(cells,x->x[1]<SK);  fi;

    c:=HAPRemoveCellFromRegularCWComplex(Y,dim,pos);

    Add(cells,c);

    ContractSpace(Y);

  od;
od;

#########################
#Y!.vectorField:=fail;
cpos:=PositionProperty(Y!.properties,x->x[1]="codim");
if cpos=fail then
Add(Y!.properties,["codim",dim]);
else
Y!.properties[cpos][2]:=dim;
fi;

Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
#Unbind(Y!.bnd);
#Unbind(Y!.cobnd);
#########################

Y!.criticalCells:=cells; return Filtered(cells,x->x[1]<SK);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(HAPRemoveVectorField,
function(Y);
if Y!.vectorField=fail then return true; fi;

Y!.vectorField:=fail;
Y!.inverseVectorField:=fail;
Y!.criticalCells:=fail;
Y!.nrCells:=function(k);
            if k>EvaluateProperty(Y,"dimension") then return 0; fi;
            return Length(Filtered(Y!.boundaries[k+1],x->not x[1]=0));
            end;
Unbind(Y!.bnd);
Unbind(Y!.cobnd);

end);
##########################################################
##########################################################


