#RT:=0;
########################################################
########################################################
InstallGlobalFunction(HomToInt_ChainComplex,
function(C)
local D, DimensionD, BoundaryD, BoundaryDrec, len, n, k,i, A;

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
BoundaryDrec[n+1]:=A;
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
                ["characteristic",0 ] ]));

return D;

end);
###############################################################
###############################################################

########################################################
########################################################
InstallGlobalFunction(HomToInt_CochainComplex,
function(C)
local D, DimensionD, BoundaryD, BoundaryDrec, len, n, k,i, A;

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
BoundaryDrec[n]:=A;
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
                ["characteristic",0 ] ]));

return D;

end);
###############################################################
###############################################################


##########################################
##########################################
InstallGlobalFunction(HomToInt_ChainMap,
function(arg)
local F,S, T, HS, HT, HF, HThomHS, Frec, zero,n,sparsemap, A, B, InitA;

F:=arg[1];
S:=Source(F);
HS:=HomToIntegers(S);
T:=Target(F);
HT:=HomToIntegers(T);

zero:=[];
for n in [1..Length(S)+1] do
zero[n]:=0*[1..S!.dimension(n-1)];
od;

A:=[];
InitA:=function(n);
A[n+1]:=IdentityMat(S!.dimension(n));
B:=List(A[n+1], r->F!.mapping(r,n));
A[n+1]:=TransposedMat(B);
end;


Frec:=List([0..Length(S)],i->[]);
#################
HThomHS:=function(v,n)
local w, ww, i, j;

if S!.dimension(n)<20001 then               #THIS IS A VERY ARBITRARY CUT-OFF: IT's AN ATTEMPT TO AVOID LARGE MATRICES
if not IsBound(A[n+1]) then InitA(n); fi;
return v*A[n+1];
fi;

w:=[];
for i in [1..S!.dimension(n)] do
if not IsBound(Frec[n+1][i]) then
zero[n+1][i]:=1;
ww:=F!.mapping(zero[n+1],n);;
Add(w,v*ww);
zero[n+1][i]:=0;
Frec[n+1][i]:=[ Filtered([1..Length(ww)],a->not IsZero(ww[a])), Filtered(ww,a->not IsZero(a))  ];
else
ww:=0*v;
for j in [1..Length(Frec[n+1][i][1])] do
ww[Frec[n+1][i][1][j]]:=1*Frec[n+1][i][2][j];
od;
Add(w,v*ww);
fi;
od;
return w;

end;
#################

#################
sparsemap:=function(v,n)
local ans, rowA, k, f, x;

rowA:=1*zero[n+1];
ans:=1*zero[n+1];;

for k in [1..Length(ans)] do
rowA[k]:=1;
#RT:=RT-Runtime();
f:= F!.mapping(rowA,n);
#RT:=RT+Runtime();
rowA[k]:=0;

  for x in v do
  ans[k]:=ans[k]+x[2]*f[x[1]];  
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
