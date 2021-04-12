
##########################################################
##########################################################
InstallGlobalFunction(BarComplexOfMonoid,
function(M,N)
local Dimension,Boundary,Properties,ord,pos, invpos, MT;
#This function returns the first N terms of the *normalized* bar chain 
#complex of a finite monoid M with multiplication table MT. It is returned
#as a sparse chain complex.

if not (IsMonoid(M) and IsFinite(M)) then
Print("The first argument must be a finite monoid.\n");
return fail;
fi;
MT:=MultiplicationTable(M);

ord:=Length(MT)-1;
#ord=|M|-1.
#We'll "re-index" MT from 0 to |M|-1 so that elements of the monoid
#correspond to elements of the list [0..|M|-1].
#We assume that 0 is the identity element in the re-indexed MT.

###############################
pos:=function(x)
local len, p, i;
#Inputs a list [a,b,c,...] with 0<=a,b,c,... 
#and returns the position on MxMxMx... with position starting at 0.
len:=Length(x);
p:=0;
for i in [1..len] do
p:=p+x[i]*ord^(len-i);
od;
return p;
end;
###############################

###############################
invpos:=function(k,d)
local x, m, q, l,i;
#The inverse to pos.
m:=k;
x:=List([1..d],i->0);

for i in Reversed([0..d-1]) do
q:=EuclideanQuotient(m,ord^i);
x[i+1]:=q; 
m:=m-q*ord^i;
od;

return Reversed(x);
end;
###############################

###############################
###############################
Dimension:=function(n)
return ord^n;
end;
###############################

###############################
Boundary:=function(n,k)
local x, y, bnd, i, a;
bnd:=[];
if n<=1 then return bnd; fi;
x:=invpos(k-1,n);
y:=x{[2..n]};
Add(bnd,[1+pos(y),1]);
for i in [1..n-1] do
a:=MT[x[i]+2][x[i+1]+2]-2;
if not a<0 then
y:=Concatenation(x{[1..i-1]},[a],x{[i+2..n]});
##Remember: We have "re-indexed" so that elements of M correspond
##to the set [0..Ord-1], and an tuple containing a zero is degenerate
##and thus not included in the normalized bar complex.
Add(bnd,[1+pos(y),(-1)^i]);
fi;
od;
y:=x{[1..n-1]};
Add(bnd,[1+pos(y),(-1)^n]);
return bnd;
#return AlgebraicReduction(bnd);

end;
###############################



return 
Objectify(HapSparseChainComplex,
           rec(
           dimension:=Dimension,
           boundary:=Boundary,
           properties:=[
           ["length",N],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));


end);
##########################################################
##########################################################

