
############################################
############################################
InstallGlobalFunction(HAP_4x4MatTo2x2Mat,
function(M,dd)
local d, fun, x, A, B, C, D;

if IsString(dd) then d:=EvalString(dd);
else d:=dd;
fi;

######Schoennenbeck transormation###
if  (d mod 4 =2) or (d mod 4 =3) then x:=Sqrt(d); fi;
if  (d mod 4 =1)   then x:=(1 + Sqrt(d))/2; fi;
####################################

###########################################
fun:=function(a,b,c,d)
local A, p, q;
#return p +q*x
A:=1*[a,b,c,d];
p:=a;
q:=SignInt(b)*Gcd(b,c,a-d);
return p + q*x;
end;
###########################################

A:=fun(M[1][1],M[1][2],M[2][1],M[2][2]) ;
B:=fun(M[1][3],M[1][4],M[2][3],M[2][4]) ;
C:=fun(M[3][1],M[3][2],M[4][1],M[4][2]) ;
D:=fun(M[3][3],M[3][4],M[4][3],M[4][4]) ;
return [[A,B],[C,D]];
end);
############################################
############################################

############################################
############################################
InstallGlobalFunction(HAP_nxnMatTo2nx2nMat,
function(M,d)
local x, fun, A, B, I, CM, rt, i, j;

######Check d#######################
if (d mod 4 =3) or (d mod 4 =2) then x:=Sqrt(d); fi;
if (d mod 4 =1)  then x:=(1 + Sqrt(d))/2; fi;
####################################

I:=IdentityMat(2);
CM:=CompanionMat(MinimalPolynomial(Rationals,x));
CM:=TransposedMat(CM);
rt:=Sqrt(-d)^-1;

if (d mod 4 =3) or (d mod 4 =2) then 
####################################
fun:=function(z)
local a,b;
b:=(ImaginaryPart(z)*rt);
a:=z-b*x;

return a*I+b*CM;
end;
####################################
fi;
if (d mod 4 =1) then 
####################################
fun:=function(z)
local a,b;
b:=(2*ImaginaryPart(z)*rt);
a:=z-b*x;

return a*I+b*CM;
end;
####################################
fi;

A:=0*IdentityMat(2*Length(M));
for i in [1..Length(M)] do
for j in [1..Length(M)] do
B:=fun(M[i][j]);
A[2*i-1][2*j-1]:=B[1][1];
A[2*i-1][2*j]:=B[1][2];
A[2*i][2*j-1]:=B[2][1];
A[2*i][2*j]:=B[2][2];
od;
od;

return A;
end);
############################################
############################################


########################################################
########################################################
InstallGlobalFunction(HomogeneousPolynomials_Bianchi,
function(arg)
local G,n,m,M,R,module,fn,x,y,u,v,z, gensM, gensG, ispsl, d;

G:=arg[1];
n:=arg[2];
m:=arg[3];

ispsl:=Name(G){[1,2,3]} = "PSL";
d:=G!.bianchiInteger;

if d<0 then
if d mod 4 =1 then R:=Field((1+Sqrt(d))/2); fi;
if d mod 4 =2 or d mod 4 = 3 then R:=Field(Sqrt(d)); fi;
fi;
if d>0 then R:=Field(Sqrt(d)); fi;

x:=Indeterminate(R,1);;
y:=Indeterminate(R,2);;
u:=Indeterminate(R,3);;
v:=Indeterminate(R,4);;
z:=Indeterminate(R,5);;

####################################
####################################
fn:=function(AA)
local A,A1,B,i,j,ii,jj,p,b,q,indet,k;
if ispsl then
A:=HAP_4x4MatTo2x2Mat(AA!.element[1],d);
else
A:=HAP_4x4MatTo2x2Mat(AA,d);
fi;
B:=[];
A1:=[ [ComplexConjugate(A[1][1]) , ComplexConjugate(A[1][2])], 
      [ComplexConjugate(A[2][1]) , ComplexConjugate(A[2][2])]];;
for i in [0..n] do
for j in [0..m] do
b:=List([0..n],i->[0..m]*0);
p:=x^(n-i)*y^i*u^(m-j)*v^j*z;
p:=Value(p,[x,y,u,v,z],
[ A[2][2]*x-A[1][2]*y , -A[2][1]*x+A[1][1]*y , 
  A1[2][2]*u-A1[1][2]*v , -A1[2][1]*u+A1[1][1]*v , z]);
p:=p!.ExtRepPolynomialRatFun;
for k in [1..Length(p)/2] do
q:=p[2*k-1];
##ii is the exponent of x.
ii:=Position(q{2*[1..Length(q)/2]-1}, 1 );
if ii=fail then ii:=0; else ii:=q[2*ii]; fi;
##jj is the exponent of u.
jj:=Position(q{2*[1..Length(q)/2]-1}, 3 );
if jj=fail then jj:=0; else jj:=q[2*jj]; fi;
b[ii+1][jj+1]:=p[2*k];
od;
Add(B,Flat(b));
od;
od;

return HAP_nxnMatTo2nx2nMat(Reversed(TransposedMat(B)),d);
end;
####################################
####################################

gensG:=GeneratorsOfGroup(G);

gensM:=List(gensG,fn);;

M:=Group( gensM );;
#module:=MappingByFunction(SL(4,Integers),M,x->fn(x));
module:=MappingByFunction(G,M,x->fn(x));

return module;
end);
########################################################
########################################################

