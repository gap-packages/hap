########################################################
########################################################
InstallGlobalFunction(HomogeneousPolynomials,
function(arg)
local G,n,m,M,module,fn,x,y,z;

G:=arg[1];
n:=arg[2];
if Length(arg)>2 then m:=arg[3]; else m:=0; fi;

if IsBound(G!.bianchiInteger) then
if IsInt(G!.bianchiInteger) then
return HomogeneousPolynomials_Bianchi(G,n,m); 
fi;
fi;
n:=n+m;

x:=Indeterminate(Rationals,1);;
y:=Indeterminate(Rationals,2);;
z:=Indeterminate(Rationals,3);;

####################################
####################################
fn:=function(A)
local B,i,ii,p,b,q,indet,k;
B:=[];

for i in [0..n] do
b:=[0..n]*0;
p:=x^(n-i)*y^i*z;
p:=Value(p,[x,y,z],[ A[2][2]*x-A[1][2]*y , -A[2][1]*x+A[1][1]*y , z]);
p:=p!.ExtRepPolynomialRatFun;
for k in [1..Length(p)/2] do
q:=p[2*k-1];
##ii is the exponent of x.
ii:=Position(q{2*[1..Length(q)/2]-1}, 1 );
if ii=fail then ii:=0; else ii:=q[2*ii]; fi;
b[ii+1]:=p[2*k];
od;
Add(B,b);
od;
return Reversed(TransposedMat(B));
end;
####################################
####################################

M:=Group( List(GeneratorsOfGroup(G),fn) );;
module:=GroupHomomorphismByFunction(G,M,x->fn(x));

return module;
end);
########################################################
########################################################
