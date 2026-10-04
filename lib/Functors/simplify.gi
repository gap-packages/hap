
################################################################
################################################################
InstallGlobalFunction(SimplifiedSparseChainComplex,
function(arg)
local C,bounds, cobounds, lb, n, k, i, j, b,c,B,x,
      Dimension,Boundary,first,bnd,Replace,NormForm,
      NewGens,ZeroCells,BoundaryRec,PNG,merge,D,CmapD,DmapC,
      CmappingD,DmappingC,ReplacementRec,rr,chnmp;

C:=arg[1];
if Length(arg)=1 then rr:=infinity; chnmp:=false; fi;
if Length(arg)=2 then
          if IsInt(arg[2]) then rr:=arg[2]; chnmp:=false; fi;
          if IsBool(arg[2]) then rr:=infinity; chnmp:=true; fi;
fi;
if Length(arg)=3 then
          if IsInt(arg[2]) then rr:=arg[2]; chnmp:=true; fi;
          if IsInt(arg[3]) then rr:=arg[3]; chnmp:=true; fi;
fi;

####################Changed this function, July 2026.  Then changed remaining
####################code so that NormForm is rarely called!
NormForm:=function(b)
local S,a,pos,ls,L,bool,i;
#cnt:=cnt+1;
bool:=true;
for i in [1..Length(b)-1] do
   if not b[i][1]<b[i+1][1] then bool:=false; break; fi;
od;
if bool then return b; fi;
Sort(b);;
S:=SSortedList(List(b,x->x[1]));
ls:=Length(S);
if ls=Length(b) then return b; fi;
Apply(S,x->[x,0]);
pos:=1;
for i in [1..Length(b)] do
a:=b[i];
S[pos][2]:=S[pos][2]+a[2];
if pos<ls then
if a[1]<b[i+1][1]  then pos:=pos+1; fi;
fi;
od;
S:=Filtered(S,x->not x[2]=0);
return S;
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
#bounds[n][k]:=NormForm(1*C!.boundary(n,k));
bounds[n][k]:=C!.boundary(n,k);
if not IsSortedList(bounds[n][k]) then bounds[n][k]:=NormForm(bounds[n][k]); fi;
#WARNING: I suppose bounds[n][k] could be sorted but yet not in normal form!!!
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
cobounds[n+1]:=List([1..Length(bounds[n])],i->[]);
####################
####################


if rr=infinity then
####################
####################
first:=function(x) #find first cell with coefficient equal to +/-1
return PositionProperty(x,y->AbsInt(y[2])=1);
end;
####################
####################
fi;
if rr<infinity then
####################
####################
first:=function(x) #find first cell with coefficient equal to +/-1
if Length(x)>rr then return fail; fi;
return PositionProperty(x,y->AbsInt(y[2])=1);
end;
####################
####################
fi;

if chnmp then
ReplacementRec:=List([0..Length(C)],n->[]);
# ReplacementRec[i] is a record of the ordered sequence of replacements of 
# n-cells, where bnd occurs in the boundary of the k-th (n+1)-cell.  
# The (n+1)-cell will be deleted by other code.
####################
####################
##This function replaces the (n-1)-cell b throughout by bnd
Replace:=function(n,b,bnd)
local cbnd, pos, B,BB,x,Y,z,i,c;
Add(ReplacementRec[n-1],[b,bnd,k]);   #k is defined whenever Replace is applied
cbnd:=cobounds[n][b];  #removed 1*
for i in cbnd do
   B:=bounds[n][i];  #removed 1*
   if not (B=0 or B=[]) then
      pos:=PositionProperty(B,a->a[1]=b);
      #pos:=PositionSet(List(B,a->a[1]),b);
      if IsInt(pos) then
          c:=B[pos][2];
          Y:=1*bnd;
          Apply(Y,a->[a[1],c*a[2]]);
          Remove(B,pos);
          bounds[n][i]:=merge(B,Y);
          for z in Y do
              AddSet(cobounds[n][z[1]],i);
          od;
      fi;
   fi;
od;
return true;
end;
####################
####################
else
####################
####################
Replace:=function(n,b,bnd)
local cbnd, pos, B,BB,x,Y,z,i,c;
cbnd:=cobounds[n][b];  #removed 1*
for i in cbnd do
   B:=bounds[n][i];  #removed 1*
   if not (B=0 or B=[]) then
      pos:=PositionProperty(B,a->a[1]=b);
      #pos:=PositionSet(List(B,a->a[1]),b);
      if IsInt(pos) then
          c:=B[pos][2];
          Y:=1*bnd;
          Apply(Y,a->[a[1],c*a[2]]);
          Remove(B,pos);
          bounds[n][i]:=merge(B,Y);
          for z in Y do
              AddSet(cobounds[n][z[1]],i);
          od;
      fi;
   fi;
od;
return true;
end;
####################
####################
fi;

####################
####################
#We want to merge sort B and Y assuming that B and Y are sorted
#returns NormForm(Concatenation(B,Y));
merge:=function(B,Y)
local U, pB, pY, lB, lY;

lB:=Length(B); lY:=Length(Y);
U:=[];
pB:=1; pY:=1;
while pB<=lB and pY<=lY do
   if B[pB][1]<Y[pY][1] then
      Add(U,B[pB]); pB:=pB+1;
   elif B[pB][1]>Y[pY][1] then
      Add(U,Y[pY]); pY:=pY+1;
   else
      if B[pB][2]<>-Y[pY][2] then
         Add(U,[B[pB][1], B[pB][2]+Y[pY][2]]);
      fi;
      pB:=pB+1; pY:=pY+1;
   fi;
od;
while pB<=lB do Add(U,B[pB]); pB:=pB+1; od;
while pY<=lY do Add(U,Y[pY]); pY:=pY+1; od;

return U;
end;
####################
####################


##################Workhorse###########################
for n in [1..Length(bounds)] do
for k in [1..Length(bounds[n])] do
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
cobounds[n+1][k]:=[];
fi;
cobounds[n][b[1]]:=[];

##################
fi;
od;
od;
######################################################


NewGens:=[];  #NewGens[n+1] will be the n-gens remaining in simplified complex
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

PNG:=[];
for n in [1..Length(bounds)] do
PNG[n]:=List([1..C!.dimension(n-1)],i->Position(NewGens[n],i));
od;

###################################
Boundary:=function(n,k)
if n>Length(bounds) then return []; fi;
return 
#List(bounds[n][k], x->[Position(NewGens[n],x[1]),x[2]]);
List(bounds[n][k], x->[PNG[n][x[1]],x[2]]);
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
#Unbind(bounds);
#Unbind(cobounds);
Unbind(ZeroCells);
###################################
Boundary:=function(n,k)
if n>lb then return []; fi;
return BoundaryRec[n][k];
end;
###################################

D:=  Objectify(HapSparseChainComplex,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                bounds:=bounds,
                cobounds:=cobounds,
                properties:=
                [["length",EvaluateProperty(C,"length")],
                ["type", "chainComplex"],
                ["characteristic",
                EvaluateProperty(C,"characteristic")] ]));

if chnmp then
CmappingD:=fail;  #To do
DmappingC:=fail;

CmapD:=Objectify(HapChainMap,
        rec(
           source:=C,
           target:=D,
           mapping:=CmappingD,
           properties:=[ ["type","chainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(C,"characteristic"),
           EvaluateProperty(C,"characteristic"))],
           ]
           ));
DmapC:=Objectify(HapChainMap,
        rec(
           source:=D,
           target:=C,
           mapping:=DmappingC,
           properties:=[ ["type","chainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(C,"characteristic"),
           EvaluateProperty(C,"characteristic"))],
           ]
           ));

D!.ReplacementRec:=ReplacementRec;
D!.CmapD:=CmapD;
D!.DmapC:=DmapC;
fi;

return D;
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

