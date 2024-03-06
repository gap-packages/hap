#RT:=0;
########################################################
########################################################
InstallGlobalFunction(HAP_HomToIntModP_ChainComplex,
function(C,p)
local D, DimensionD, BoundaryD, BoundaryDrec, len, n, k,i, A, one;

if IsPrimeInt(p) then
one:=One(GF(p));
else
one:=One(Field(1));
fi;
################################
################################
if not IsHapChainComplex(C) then
Print("Error: input must be a chain complex.\n");
return fail;
fi;

if not EvaluateProperty(C,"characteristic")=0 then
Print("Error: chain complex must have characteristic 0.\n"); 
return fail;
fi;
################################
################################

len:=Length(C);
BoundaryDrec:=[];

################################
################################
DimensionD:=function(n);
if n<0 or n>len then return 0; fi;

return C!.dimension(n);
end;
################################
################################

#for n in [1..len] do
#for k in [1..C!.dimension(n)] do
#C!.boundary(n,k);
#od;
#od;

#######################
BoundaryD:=function(n,k);
if n<0 or n>=len then return [0]; fi;

if not IsBound(BoundaryDrec[n+1]) then
A:=BoundaryMatrix(C,n+1);
if Length(A)=0 then A:=List([1..C!.dimension(n)],x->[0]); fi;   
BoundaryDrec[n+1]:=A*one;
fi;

return BoundaryDrec[n+1][k];
end;
#######################

D:=Objectify(HapCochainComplex,
                rec(
                dimension:=DimensionD,
                boundary:=BoundaryD,
                properties:=
                [["length",len],
                ["type", "cochainComplex"],
                ["characteristic",p ] ]));

return D;

end);
###############################################################
###############################################################

########################################################
########################################################
InstallGlobalFunction(HAP_HomToIntModP_CochainComplex,
function(C,p)
local D, DimensionD, BoundaryD, BoundaryDrec, len, n, k,i, A, one;

if IsPrimeInt(p) then
one:=One(GF(p));
else
one:=One(Field(1));
fi;
################################
################################
if not IsHapCochainComplex(C) then
Print("Error: input must be a cochain complex.\n");
return fail;
fi;

if not EvaluateProperty(C,"characteristic")=0 then
Print("Error: cochain complex must have characteristic 0.\n");
return fail;
fi;
################################
################################

len:=Length(C);
BoundaryDrec:=[];

################################
################################
DimensionD:=function(n);
if n<0 or n>len then return 0; fi;

return C!.dimension(n);
end;
################################
################################

#for n in [1..len] do
#for k in [1..C!.dimension(n)] do
#C!.boundary(n,k);
#od;
#od;

#######################
BoundaryD:=function(n,k);
if n<0 or n>len then return [0]; fi;

if not IsBound(BoundaryDrec[n]) then
A:=BoundaryMatrix(C,n-1);
if Length(A)=0 then A:=List([1..C!.dimension(n)],x->[0]); fi;
BoundaryDrec[n]:=A*one;
fi;

return BoundaryDrec[n][k];
end;
#######################

D:=Objectify(HapChainComplex,
                rec(
                dimension:=DimensionD,
                boundary:=BoundaryD,
                properties:=
                [["length",len],
                ["type", "chainComplex"],
                ["characteristic",p ] ]));

return D;

end);
###############################################################
###############################################################


##########################################
##########################################
InstallGlobalFunction(HAP_HomToIntModP_ChainMap,
function(F,p)
local S, T, HS, HT, HF, HThomHS, zero,n,sparsemap, A, B, InitA, one;

if IsPrimeInt(p) then
one:=One(GF(p));
else one:=One(Field(1));
fi;
S:=Source(F);
HS:=HomToIntegersModP(S,p);
T:=Target(F);
HT:=HomToIntegersModP(T,p);

zero:=[];
for n in [1..Length(S)+1] do
zero[n]:=0*[1..S!.dimension(n-1)];
od;

#A:=List([0..Length(S)],i->IdentityMat(S!.dimension(i)));
A:=[];
#for n in [0..Length(A)-1] do
InitA:=function(n);
A[n+1]:=IdentityMat(S!.dimension(n));
B:=List(A[n+1], r->F!.mapping(r,n));
A[n+1]:=TransposedMat(B);
if Length(A[n+1])=0 then A[n+1]:=[[0*one]]; fi;
end;
#od;
#################
HThomHS:=function(v,n)
local rowA, ans, k;
if not IsBound(A[n+1]) then InitA(n); fi;
return v*A[n+1];
end;
#################

#################
sparsemap:=function(v,n)
local ans, rowA, k, f, x;

rowA:=StructuralCopy(zero[n+1]);
ans:=StructuralCopy(zero[n+1]);;

for k in [1..Length(ans)] do
rowA[k]:=1;
f:= F!.mapping(rowA,n);
rowA[k]:=0;
  for x in v do
  ans[k]:=ans[k]+x[2]*f[x[1]]*SignInt(x[1]);
  od;
od;
return ans;
end;
#################


return Objectify(HapCochainMap,
        rec(
           source:=HT,
           target:=HS,
           mapping:=HThomHS,
           sparseMap:=sparsemap,
           properties:=[ ["type","cochainMap"],
           ["characteristic", Maximum(
           EvaluateProperty(F!.source,"characteristic"),
           EvaluateProperty(F!.target,"characteristic"))]
           ]));

end);
##########################################
##########################################



##########################################
##########################################
InstallGlobalFunction(HAP_HomToIntModP_CochainMap,
function(F,p)
local char, S, T, HS, HT, HF, HThomHS, zero,n,sparsemap, A, B, InitA, one;

if IsPrimeInt(p) then
one:=One(GF(p));
else one:=One(Field(1));
fi;
S:=Source(F);
HS:=HomToIntegersModP(S,p);
T:=Target(F);
HT:=HomToIntegersModP(T,p);

zero:=[];
for n in [1..Length(S)+1] do
zero[n]:=0*[1..S!.dimension(n-1)];
od;

#A:=List([0..Length(S)],i->IdentityMat(S!.dimension(i)));
A:=[];
#for n in [0..Length(A)-1] do
####################
InitA:=function(n)
local B;
A[n+1]:=IdentityMat(S!.dimension(n));
#RT:=RT-Runtime();
B:=List(A[n+1], r->F!.mapping(r,n));
#RT:=RT+Runtime();
A[n+1]:=TransposedMat(B);
end;
####################
#od;
#################
HThomHS:=function(v,n);
if not IsBound(A[n+1]) then InitA(n); fi;
 return v*A[n+1];
end;
#################


return Objectify(HapChainMap,
        rec(
           source:=HT,
           target:=HS,
           mapping:=HThomHS,
           properties:=[ ["type","chainMap"],
                         ["characteristic", p] ]));

end);
##########################################
##########################################

