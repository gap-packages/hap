
################################################################
################################################################
InstallGlobalFunction(SimplifiedSparseChainComplex,
function(arg)
local C,bounds, cobounds, lb, n, k, i, j, b,c,B,x,
      Dimension,Boundary,first,bnd,Replace,NormForm,
      NewGens,ZeroCells,BoundaryRec,F;

C:=arg[1];

####################Changed this function, July 2026
####################
NormForm:=function(b)
local S,a,pos,ls;
Sort(b);;
S:=SSortedList(List(b,x->x[1]));
ls:=Length(S);
if ls=Length(b) then return b; fi;
S:=List(S,x->[x,0]);
pos:=1;
for i in [1..Length(b)] do
a:=b[i];
S[pos][2]:=S[pos][2]+a[2];
if pos<ls then
if a[1]<b[i+1][1]  then pos:=pos+1; fi;
fi;
od;
return Filtered(S,x->not x[2]=0);
end;
####################
####################

####################
####################
bounds:=List([1..Length(C)],i->[]);
cobounds:=List([1..Length(C)],i->[]);
ZeroCells:=[1..C!.dimension(0)];

for n in [1..Length(C)] do
for k in [1..C!.dimension(n)] do
bounds[n][k]:=NormForm(1*C!.boundary(n,k));
od;
od;

for n in [1..Length(C)] do
cobounds[n]:=List([1..C!.dimension(n-1)],i->[]);
for k in [1..C!.dimension(n)] do
b:=C!.boundary(n,k);
for x in b do
Add(cobounds[n][x[1]],k);
od;
od;
od;

for n in [1..Length(cobounds)] do
for k in [1..Length(cobounds[n])] do
cobounds[n][k]:=SSortedList(cobounds[n][k]);
od;od;


####################
####################

if Length(arg)=1 then
####################
####################
first:=function(x) #find first cell with coefficient equal to +/-1
return PositionProperty(x,y->AbsInt(y[2])=1);
#local i;
#for i in [1..Length(x)] do
#if AbsInt(x[i][2])=1 then return i; fi;
#od;
#return fail;
end;
####################
####################
fi;
if Length(arg)=2 then
####################
####################
first:=function(x) #find first cell with coefficient equal to +/-1
if Length(x)>arg[2] then return fail; fi;
return PositionProperty(x,y->AbsInt(y[2])=1);
#local i;
#for i in [1..Length(x)] do
#if AbsInt(x[i][2])=1 then return i; fi;
#od;
#return fail;
end;
####################
####################
fi;


####################
####################
Replace:=function(n,b,bnd)
local cbnd, i, B,BB,x,y,z;
cbnd:=1*cobounds[n][b];
for i in cbnd do
B:=bounds[n][i];
if not B=0 then 
BB:=[];
  for x in B do
  if not x[1]=b then
  Add(BB,x); 
  else
  y:=1*bnd;
  Apply(y,a->[a[1],x[2]*a[2]]);
  Append(BB,y);
   for z in y do
   AddSet(cobounds[n][z[1]],i);
   od;
  fi;
  od;
bounds[n][i]:=NormForm(BB);
fi;
od;
return true;
end;
####################
####################






##################Workhorse###########################
for n in [1..Length(bounds)] do
#for k in [1..Length(bounds[n])] do
F:=[1..Length(bounds[n])];
SortBy(F , k->Length(bounds[n][k]));
for k in F do
i:=first(bounds[n][k]);
if not i=fail then
##################
bnd:=1*Concatenation(bounds[n][k]{[1..i-1]},bounds[n][k]{[i+1..Length(bounds[n][k])]});
b:=1*bounds[n][k][i];
if n>1 then bounds[n-1][b[1]]:=0; fi;
if n=1 then RemoveSet(ZeroCells,b[1]); fi;
if b[2]=1 then Apply(bnd,x->[x[1],-x[2]]); fi;
Replace(n,b[1],bnd);#replace the n-1-cell [b[1],1] by bnd;
bounds[n][k]:=0;


if n<Length(cobounds) then
for j in cobounds[n+1][k] do
bounds[n+1][j]:=Filtered(bounds[n+1][j],x->not x[1]=k);
od;
fi;

if n<Length(cobounds) then
cobounds[n+1][k]:=0;
fi;
cobounds[n][b[1]]:=0;

##################
fi;
od;
od;
######################################################

NewGens:=[];  #NewGens[n+1] will be the n-gens that remain in simplified complex
NewGens[1]:=ZeroCells;
for n in [1..Length(bounds)] do
NewGens[n+1]:=Filtered([1..Length(bounds[n])],k-> not bounds[n][k]=0);
od;

for n in [1..Length(bounds)] do
bounds[n]:=Filtered(bounds[n],i->not i=0);
od;


###################################
Dimension:=function(n)
if n<0 or n>=Length(NewGens) then return 0; fi;
return Length(NewGens[n+1]);
end;
###################################

###################################
Boundary:=function(n,k)
if n>Length(bounds) then return []; fi;
return 
List(bounds[n][k], x->[Position(NewGens[n],x[1]),x[2]]);
end;
###################################

BoundaryRec:=[];
for n in [1..Length(bounds)] do
BoundaryRec[n]:=[];
for k in [1..Dimension(n)] do
BoundaryRec[n][k]:=1*Boundary(n,k);
od;
od;

lb:=Length(bounds);
Unbind(bounds);
Unbind(cobounds);
Unbind(ZeroCells);
###################################
Boundary:=function(n,k)
if n>lb then return []; fi;
return BoundaryRec[n][k];
end;
###################################

return  Objectify(HapSparseChainComplex,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=
                [["length",EvaluateProperty(C,"length")],
                ["type", "chainComplex"],
                ["characteristic",
                EvaluateProperty(C,"characteristic")] ]));


end);
################################################################
################################################################


###########################################################
###########################################################
InstallGlobalFunction(ContractedComplexViaChild,
function(arg)
local  C,r,n,L,bool,file,tmpdir,t,cmd,D;

C:=arg[1];
if Length(arg)=2 then r:=arg[2]; else r:=10^10; fi;  #SLOPPY!

##First check to see if any boundaries have length <=r ####
for n in [1..Length(C)] do
L:=List([1..C!.dimension(1)],k->Length(C!.boundary(1,k)) );;
L:=Filtered(L,x->not x=0);;
if Minimum(L) <=r then bool:=true; fi;
od;
if not bool then return C; fi;
###########################################################

tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "complex.txt" );

t:=ChildProcess();;

HAPPrintTo(file,C);
NextAvailableChild([t]);
cmd:=Concatenation("C:=HAPRead(\"",file,"\");");
#ChildCommand("C:=HAPRead(\"/tmp/sparse.txt\");",t);
ChildCommand(cmd,t);
cmd:=Concatenation("D:=ContractedComplex(C,",String(r),");");
NextAvailableChild([t]);
ChildCommand(cmd,t);
NextAvailableChild([t]);
cmd:=Concatenation("HAPPrintTo(\"",file,"\",D);");
#ChildCommand("HAPPrintTo(\"/tmp/sparse.txt\",D);",t);
ChildCommand(cmd,t);
NextAvailableChild([t]);
D:=HAPRead(file);

ChildClose(t);
return D;
end);
##########################################################
##########################################################

