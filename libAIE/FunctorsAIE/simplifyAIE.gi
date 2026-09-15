################################################################
################################################################
# simplify_optimized2.gi
#
# Optimized implementation of SimplifiedSparseChainComplex.
# The external calling convention and returned HapSparseChainComplex
# are kept compatible with the original implementation.
################################################################

InstallGlobalFunction(SimplifiedSparseChainComplex,
function(arg)
local C,bounds,cobounds,lb,n,k,i,j,b,c,B,x,bnd,Replace,NormForm,
      NewGens,ZeroCells,BoundaryRec,PNG,merge,FindCellPos,
      FindUnitPos,NegBnd,Dimension,Boundary,first,unitpos;

C:=arg[1];

################################################################
# Normalize a sparse boundary.
#
# IMPORTANT: retain the original external semantics: NormForm is
# only called when IsSortedList(b) is false.  Thus a sorted list
# with repeated first coordinates is left alone, as before.
################################################################
NormForm:=function(b)
local S,a,pos,ls,L,bool,i;
  bool:=true;
  for i in [1..Length(b)-1] do
    if not b[i][1]<b[i+1][1] then
      bool:=false;
      break;
    fi;
  od;
  if bool then return b; fi;

  Sort(b);
  S:=SSortedList(List(b,x->x[1]));
  ls:=Length(S);
  if ls=Length(b) then return b; fi;
  Apply(S,x->[x,0]);
  pos:=1;
  for i in [1..Length(b)] do
    a:=b[i];
    S[pos][2]:=S[pos][2]+a[2];
    if pos<ls then
      if a[1]<b[i+1][1] then pos:=pos+1; fi;
    fi;
  od;
  S:=Filtered(S,x->not x[2]=0);
  return S;
end;

################################################################
# Store sparse boundaries.
################################################################
bounds:=List([1..Length(C)],i->[]);
cobounds:=List([1..Length(C)],i->[]);
ZeroCells:=[1..C!.dimension(0)];

for n in [1..Length(C)] do
  for k in [1..C!.dimension(n)] do
    bounds[n][k]:=C!.boundary(n,k);
    if not IsSortedList(bounds[n][k]) then
      bounds[n][k]:=NormForm(bounds[n][k]);
    fi;
  od;
od;

################################################################
# Build reverse incidence lists.
#
# Since k is traversed increasingly, each coboundary list is
# already sorted; consequently SSortedList is unnecessary.
################################################################
for n in [1..Length(C)] do
  cobounds[n]:=List([1..C!.dimension(n-1)],i->[]);
  for k in [1..C!.dimension(n)] do
    for x in bounds[n][k] do
      Add(cobounds[n][x[1]],k);
    od;
  od;
  cobounds[n+1]:=List([1..Length(bounds[n])],i->[]);
od;

################################################################
# Find the first coefficient of absolute value 1.
#
# Use a simple while loop: constructing [1..Length(x)] for every
# boundary is surprisingly expensive in large GAP computations.
################################################################
FindUnitPos:=function(x)
local q,l;
  l:=Length(x);
  q:=1;
  while q<=l do
    if AbsInt(x[q][2])=1 then return q; fi;
    q:=q+1;
  od;
  return fail;
end;

if Length(arg)=1 then
  first:=FindUnitPos;
elif Length(arg)=2 then
  first:=function(x)
    if Length(x)>arg[2] then return fail; fi;
    return FindUnitPos(x);
  end;
fi;

################################################################
# Cache the first +/-1 position.  A boundary only changes when it
# is touched by Replace, or when one entry is deleted from it in
# the next dimension.  Therefore we can update the cache locally
# instead of scanning every boundary in the workhorse loop.
################################################################
unitpos:=List([1..Length(bounds)],n->[]);
for n in [1..Length(bounds)] do
  unitpos[n]:=List([1..Length(bounds[n])],i->fail);
  for k in [1..Length(bounds[n])] do
    if Length(arg)=1 then
      unitpos[n][k]:=FindUnitPos(bounds[n][k]);
    elif Length(arg)=2 then
      if Length(bounds[n][k])<=arg[2] then
        unitpos[n][k]:=FindUnitPos(bounds[n][k]);
      fi;
    fi;
  od;
od;

################################################################
################################################################
# Merge two sorted sparse boundaries.
################################################################
merge:=function(B,Y,n,i)
local U,pB,pY,lB,lY,v;
  lB:=Length(B);
  lY:=Length(Y);
  U:=[];
  pB:=1;
  pY:=1;

  # Update cobounds only when a Y-cell was not already present in B.
  # In that case i is already in the corresponding coboundary list,
  # so the old unconditional AddSet loop was redundant.
  while pB<=lB and pY<=lY do
    if B[pB][1]<Y[pY][1] then
      Add(U,B[pB]);
      pB:=pB+1;
    elif B[pB][1]>Y[pY][1] then
      Add(U,Y[pY]);
      AddSet(cobounds[n][Y[pY][1]],i);
      pY:=pY+1;
    else
      v:=B[pB][2]+Y[pY][2];
      if v<>0 then Add(U,[B[pB][1],v]); fi;
      pB:=pB+1;
      pY:=pY+1;
    fi;
  od;

  while pB<=lB do
    Add(U,B[pB]);
    pB:=pB+1;
  od;
  while pY<=lY do
    Add(U,Y[pY]);
    AddSet(cobounds[n][Y[pY][1]],i);
    pY:=pY+1;
  od;

  return U;
end;

################################################################
# Replace an (n-1)-cell b by bnd in every n-cell containing b.
#
# Sparse boundaries in this application are generally short, so a
# tight linear search is faster than a general binary-search helper.
# The replacement boundary is negated only once when necessary.
################################################################
Replace:=function(n,b,bnd)
local cbnd,B,pos,c,Y,z,i,negBnd;

  cbnd:=cobounds[n][b];
  negBnd:=fail;

  for i in cbnd do
    B:=bounds[n][i];

    if B<>0 and B<>[] then
      pos:=PositionProperty(B,a->a[1]=b);

      if IsInt(pos) then
        c:=B[pos][2];

        if c=1 then
          Y:=bnd;
        elif c=-1 then
          if negBnd=fail then
            negBnd:=List(bnd,x->[x[1],-x[2]]);
          fi;
          Y:=negBnd;
        else
          Y:=List(bnd,x->[x[1],c*x[2]]);
        fi;

        Remove(B,pos);
        bounds[n][i]:=merge(B,Y,n,i);

        # This boundary has changed, so its cached pivot must be
        # recomputed.  This also handles a pivot being introduced
        # by the replacement.
        if Length(arg)=1 then
          unitpos[n][i]:=FindUnitPos(bounds[n][i]);
        elif Length(bounds[n][i])<=arg[2] then
          unitpos[n][i]:=FindUnitPos(bounds[n][i]);
        else
          unitpos[n][i]:=fail;
        fi;

      fi;
    fi;
  od;

  return true;
end;

################################################################
################################################################
# Workhorse.
################################################################
for n in [1..Length(bounds)] do
  for k in [1..Length(bounds[n])] do
    i:=unitpos[n][k];

    if IsInt(i) then
      # Remove the pivot entry without constructing two slices.
      bnd:=ShallowCopy(bounds[n][k]);
      b:=bnd[i];
      Remove(bnd,i);

      if n>1 then
        bounds[n-1][b[1]]:=0;
      fi;

      if n=1 then
        RemoveSet(ZeroCells,b[1]);
      fi;

      if b[2]=1 then
        Apply(bnd,x->[x[1],-x[2]]);
      fi;

      Replace(n,b[1],bnd);
      bounds[n][k]:=0;

      ################################################################
      # Remove k from every (n+1)-boundary containing it.
      #
      # Those boundaries are normalized, so binary search avoids a
      # full Filtered scan.
      ################################################################
      if n<Length(cobounds) then
        for j in cobounds[n+1][k] do
          B:=bounds[n+1][j];
          if B<>0 and B<>[] then
            i:=PositionProperty(B,x->x[1]=k);
            if IsInt(i) then
              Remove(B,i);
              # Removing an entry cannot create a new +/-1
              # coefficient, but it can destroy the cached one.
              if unitpos[n+1][j]=i then
                unitpos[n+1][j]:=FindUnitPos(B);
              elif IsInt(unitpos[n+1][j]) and unitpos[n+1][j]>i then
                unitpos[n+1][j]:=unitpos[n+1][j]-1;
              fi;
              if Length(arg)=2 and Length(B)>arg[2] then
                unitpos[n+1][j]:=fail;
              fi;
            fi;
          fi;
        od;
        cobounds[n+1][k]:=[];
      fi;

      cobounds[n][b[1]]:=[];
    fi;
  od;
od;

################################################################
# Build surviving-cell lists and compact boundaries.
################################################################
NewGens:=[];
NewGens[1]:=ZeroCells;

for n in [1..Length(bounds)] do
  NewGens[n+1]:=
    Filtered([1..Length(bounds[n])],k->not bounds[n][k]=0);
  bounds[n]:=Filtered(bounds[n],i->not i=0);
od;

Dimension:=function(n)
  if n<0 or n>=Length(NewGens) then return 0; fi;
  return Length(NewGens[n+1]);
end;

################################################################
# Build old-cell -> new-cell maps and final boundaries directly.
# This avoids the extra Boundary() closure and an intermediate pass.
################################################################
PNG:=[];
BoundaryRec:=[];
for n in [1..Length(bounds)] do
  PNG[n]:=List([1..C!.dimension(n-1)],i->0);
  for k in [1..Length(NewGens[n])] do
    PNG[n][NewGens[n][k]]:=k;
  od;

  BoundaryRec[n]:=[];
  for k in [1..Length(bounds[n])] do
    BoundaryRec[n][k]:=List(bounds[n][k],x->[PNG[n][x[1]],x[2]]);
  od;
od;

lb:=Length(bounds);

Boundary:=function(n,k)
  if n>lb then return []; fi;
  return BoundaryRec[n][k];
end;

Unbind(ZeroCells);

return Objectify(HapSparseChainComplex,
  rec(
    dimension:=Dimension,
    boundary:=Boundary,
    bounds:=bounds,
    cobounds:=cobounds,
    properties:=
      [["length",EvaluateProperty(C,"length")],
       ["type","chainComplex"],
       ["characteristic",EvaluateProperty(C,"characteristic")]]));
end);

################################################################
################################################################

InstallGlobalFunction(ContractedComplexViaChild,
function(arg)
local C,r,n,L,bool,file,tmpdir,t,cmd,D;

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
ChildCommand(cmd,t);
cmd:=Concatenation("D:=ContractedComplex(C,",String(r),");");
NextAvailableChild([t]);
ChildCommand(cmd,t);
NextAvailableChild([t]);
cmd:=Concatenation("HAPPrintTo(\"",file,"\",D);");
ChildCommand(cmd,t);
NextAvailableChild([t]);
D:=HAPRead(file);

ChildClose(t);
return D;
end);
################################################################
################################################################
