
###################################################################
###################################################################
InstallGlobalFunction(ChevalleyEilenbergComplexOfModule,
function(V,s)

local A,Sctab,B,BasV,Dim,n,i,d,j,Boundary,BD,BDI,bound,r,Comb,ITT,TTI,ONE,Charact,
      dimV, pair2int, int2pair, action, BoundaryRec;

dimV:=Dimension(V);                 #This is the dimension of the module V
A:=V!.LeftActingAlgebra;            #A is the Lie algebra
B:=Basis(A);
BasV:=Basis(V);
Sctab:=StructureConstantsTable(B);
d:=Length(B);                       #d = dimension of the Lie algebra A
Comb:=List([1..s],n->Combinations([1..d],n)); #The graded generators of C_n(A)
BoundaryRec:=List([0..s],i->[]);
#################################################################

# We'll try to use n for the degree of a chain complex.
# Let C_n(A) be the vector space of rank=Binomial(d,n) where the Lie algebra
# A has dimension d.
# Let -(x)- denote tensor product over ground field.
# Let c_p(x)v_q be the generator of C(n)=C_n(A)(x)V where c_p is 
# the p-th basis element of C_n(A) and v_q is the q-th basis element of V.
# Comb[n][j] corresponds to j-th generator of C_n(A)

ONE:=One(A!.LeftActingDomain);
################################################################
Charact:=function(M)
if not IsFinite(M!.LeftActingDomain) then 
	if Name(M!.LeftActingDomain)=Name(Integers) then 
		return 0;
	fi;
	if  Name(M!.LeftActingDomain)=Name(Rationals) then 
		return -(1/2);
	fi;
else
	return Characteristic(M!.LeftActingDomain); 
fi;
end;
################################################################


###############################################################
ITT:=function(n,j);
return StructuralCopy(Comb[n][j]);
end;
###############################################################

###############################################################
TTI:=function(c);
return PositionSorted(Comb[Length(c)],SSortedList(c));
end;
###############################################################

################################################################
Dim:=function(n);
if n>s or n<0 then return 0; else
return Binomial( d, n )*dimV; fi;
end;
###############################################################

#######################################
pair2int:=function(m,k);
return dimV*(m-1)+k;
end;
#######################################

#######################################
int2pair:=function(kk)
local m,k;
k:= kk mod dimV;
   if k=0 then k:=dimV; 
   m:=Int(kk/dimV);
   return [m,k];
   fi;
m:=1+Int(kk/dimV);
return [m,k];
end;
#######################################

#######################################
action:=function(a,v)   #NEEDS IMPROVING AS IT TAKES MOST OF THE TIME
local m, w, i, j, x, z, y;
m:=Length(v)/dimV;
w:=[];
for i in [0..m-1] do
x:=Zero(V);
for j in [1..dimV] do
z:=v[i*dimV+j];
if z<>0 then
x:=x + z*BasV[j];  
fi;
od;
y:=a*x;
if x=y then
Append(w,v{i*dimV+[1..dimV]});
else
x:=Coefficients(BasV,y);
Append(w,x);
fi;
od;
return w;
end;
#######################################

###############################################################
BDI:=function(n,i,j)
local bound, Q, pos;

if n>s then
        return [0]*ONE;
fi;
if j>Binomial(d,n) then
        return [0]*ONE;
fi;

if n=1 then return [1]; fi;
bound:=List([1..Binomial(d,n-1)], i->0);
Q:=ITT(n,j);
Remove(Q,i);
pos:=PositionSorted(Comb[n-1],Q);
if IsOddInt(i) then
bound[pos]:=-1;
else
bound[pos]:=1;
fi;
return bound;
end;
###############################################################

###############################################################
BD:=function(n,j)
#Boundary returns the image of the j-th generator under the map
#D:C_n(A)-->C_{n-1}(A) in the Chevalley-Eilenberg complex with trivial
#coefficients

local bound,a,b,x,p,q,R,m,Q,t;

if n>s then
        return [0]*ONE;
fi;
if j>Binomial(d,n) then
        return [0]*ONE;
fi;

bound:=List([1..Binomial(d,n-1)], i->0);
Q:=ITT(n,j);

for a in [1..n-1] do
for b in [a+1..n] do
        p:=Length(Sctab[Q[a]][Q[b]][1]);
        for m in [1..p] do
                R:=StructuralCopy(Q);
                Remove(R,a);
                Remove(R,b-1);
                if not Sctab[Q[a]][Q[b]][1][m] in R then
                AddSet(R,Sctab[Q[a]][Q[b]][1][m]);
                        q:=Position(SortedList(R),Sctab[Q[a]][Q[b]][1][m]);
                        t:=TTI(R);
                        bound[t]:=bound[t]+(-1)^(a+b+q-1)*Sctab[Q[a]][Q[b]][2][m];
                fi;
        od;
od;
od;

return bound;

end;
#############################################################

#############################################################
Boundary:=function(n,kk)
local bound, bnd, vec2, x, m, k, i, j, vec,Q;

# We'll implement the boundary map 
#   C_n(A) (x) V ---> C_{n-1}(A) (x) V 
#   c_p (x) v_q  |--> sum_i( x_i . DI(i)(c_p) (x) v_q ) + D(c_p (x) v_q) 
# where DI and D are defined above. This function returns the image of 
# the kk-th generator of C_n(A) (x) V under the boundary map.

if IsBound(BoundaryRec[n+1][AbsInt(kk)]) then
if kk>0 then return BoundaryRec[n+1][kk];
else
return -BoundaryRec[n+1][kk];
fi;
fi;

bound:=List([1..Dim(n-1)],i->0);

x:=int2pair(kk);
k:=x[2];
m:=x[1];
bnd:=BD(n,m);
for j in [1..Length(bnd)] do
x:=pair2int(j,k);
bound[x]:=bnd[j];
od;

Q:=ITT(n,m);

for i in [1..n] do
vec:=BDI(n,i,m);
vec2:=0*bound;
for j in [1..Length(vec)] do
if not vec[j]=0 then
x:=pair2int(j,k);
vec2[x]:=vec2[x]+vec[j];
vec2:=action(B[Q[i]],vec2);  
fi;
od;
bound:=bound+vec2;
od;

BoundaryRec[n+1][AbsInt(kk)]:= bound*ONE;
if kk>0 then return BoundaryRec[n+1][kk];
else
return -BoundaryRec[n+1][kk];
fi;
end;
#############################################################


return Objectify(HapChainComplex,rec(dimension:=Dim,
           boundary:=Boundary,
           properties:=
                [["length",s],
                 ["type","chainComplex"],
                 ["characteristic",Charact(A)]]));
end);
#############################################################
#############################################################


