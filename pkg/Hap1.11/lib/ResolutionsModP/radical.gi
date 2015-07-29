

###############################################
###############################################
InstallGlobalFunction(RadicalSeriesOfResolution,
function(R)
local G, ordG, FG, FGdims, Tdims, T,L,prime, 
      Dimension, FDimension, Boundary, BoundRec, IntToPair, one, 
      PairToInt, cnt, Mult, invT, n, i,g, j, k, r,
      WordToVectorList, VectorListToWord ;

G:=R!.group;
ordG:=Order(G);
prime:=Factors(ordG)[1];

###################
###################
if not (IsPrime(EvaluateProperty(R,"characteristic"))
and
IsPrimePowerInt(ordG)) then
Print("The resolution must have prime characteristic for a prime-power group\n");
return fail;
fi;
###################
###################

one:=One(GF(prime));

###################
Mult:=function(g,h);
return Position(R!.elts,R!.elts[g]*R!.elts[h]);
end;
###################

###################
FG:=GroupAlgebraAsFpGModule(G);
L:=RadicalSeriesOfFpGModule(FG);;
L:=L{[1..Length(L)-1]};
L:=Reversed(L);
T:=List(L,x->GeneratorsOfFpGModule(x));
FGdims:=List(T,Length);
Tdims:=1*FGdims;
for i in [2..Length(FGdims)] do
FGdims[i]:=FGdims[i]+FGdims[i-1];
od;
T:=Concatenation(T); ##So T is the new basis

invT:=T^-1;
###################

###################
Dimension:=function(n);
if n>Length(R) then return 0; fi;
return ordG*R!.dimension(n);
end;
###################

###################
FDimension:=function(k,n);
if n>Length(R) then return 0; fi;
return FGdims[k]*R!.dimension(n);
end;
###################


##################################
IntToPair:=List([0..Length(R)],i->[]);
PairToInt:=List([0..Length(R)],i->List([1..R!.dimension(i)], j->[]));
# so i --> [r,g] where r is the FG summand and g is the (new basis) position 
# in summand

for n in [0..Length(R)] do
cnt:=0;
   for j in [1..Length(Tdims)] do
   for g in [1..Tdims[j]] do
   for r in [1..R!.dimension(n)] do
   cnt:=cnt+1;
   IntToPair[n+1][cnt]:=[r,g];
   if j>1 then IntToPair[n+1][cnt][2]:=IntToPair[n+1][cnt][2]+FGdims[j-1]; fi;
   PairToInt[n+1][IntToPair[n+1][cnt][1]][IntToPair[n+1][cnt][2]]:=cnt;
   od;
   od;
   od;
od;
################################

#####################################################################
WordToVectorList:=function(w,k) #w is a FG-word in R_k.
local v,x,r;                    #v is a list of vectors mod p.

v:=List([1..R!.dimension(k)],i->List([1..ordG],j->0*one) );

for x in w do
r:=AbsInt(x[1]);
v[r][x[2]]:=v[r][x[2]] + SignInt(x[1])*one;
od;

return v ;
end;
#####################################################################


#####################################################################
VectorListToWord:=function(v,n)  #returns a sparse word over F.
local w, r, g, vv;

w:=[];

for r in [1..Length(v)] do
for g in [1..Length(v[r])] do
if not IsZero(v[r][g]) then Add(w,[PairToInt[n+1][r][g],v[r][g]]); fi;
od;
od;

return w;
end;
#####################################################################


BoundRec:=List([1..Length(R)],i->[]);;

###################################
Boundary:=function(n,k)
local w, x,b, bnd, pr, gg,i;

if IsBound(BoundRec[n][AbsInt(k)]) then
if SignInt(k)>0 then return 1*BoundRec[n][k]; 
else
return NegateWord(BoundRec[n][k]);
fi;
fi;
 
pr:=IntToPair[n+1][AbsInt(k)];
r:=pr[1];
g:=pr[2];

w:=T[g];
gg:=[];
for i in [1..Length(w)] do
for j in [1..IntFFE(w[i])] do
Add(gg,i);
od;
od;

w:=R!.boundary(n,r);  #This is an FG-word
bnd:=[];
for i in gg do
b:=List(w, x->[x[1],Mult(i,x[2])]);
Append(bnd,b);
od;
bnd:=WordToVectorList(bnd,n-1);

for i in [1..Length(bnd)] do
if not IsZero(bnd[i]) then
bnd[i]:=bnd[i]*invT;
fi;
od;

bnd:=VectorListToWord(bnd,n-1);
BoundRec[n][AbsInt(k)]:=bnd;

if SignInt(k)=1 then return 1*bnd;
else return NegateWord(bnd); fi;

end;
###################


return 
Objectify(HapFilteredSparseChainComplex,
           rec(
           dimension:=Dimension,
           boundary:=Boundary,
           filteredDimension:=FDimension,
           properties:=[
           ["length",Length(R)],
           ["filtration_length",Length(FGdims)],
           ["type","FilteredChainComplex"],
           ["characteristic",prime]]
           ));


end);
###############################################
###############################################

