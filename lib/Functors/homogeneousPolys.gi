
########################################################
########################################################
InstallGlobalFunction(HomogeneousPolynomials,
function(arg)
local G,n,m,M,module,fn,x,y,z,a,b,c,d,
      indsR,indsS,R,S,A,F,FF,BB,i,j,k,ln,fnn, bool,
      lnFF, lnFFk, In, pos, 2m, 2m1, EltsSS, Elts,recfnn;

bool:=ITER_POLY_WARN;
ITER_POLY_WARN:=false;
G:=arg[1];
n:=arg[2];
if Length(arg)>2 then m:=arg[3]; else m:=0; fi;
2m:=List([1..n+10],i->2*i);
2m1:=List([1..n+10],i->2*i-1);

if IsBound(G!.bianchiInteger) then
if IsInt(G!.bianchiInteger) then
return HomogeneousPolynomials_Bianchi(G,n,m);
fi;
fi;
n:=n+m;

#indsR:= ["a","b","c","d"];;
#indsS:=["x","y","z"];;
#R:=PolynomialRing(Rationals,indsR);;
#S:=PolynomialRing(R,indsS);;
R:=PolynomialRing(Rationals,4);;
S:=PolynomialRing(R,5);;
indsR:=R!.IndeterminatesOfPolynomialRing;
indsS:=S!.IndeterminatesOfPolynomialRing;
a:=indsR[1];
b:=indsR[2];
c:=indsR[3];
d:=indsR[4];
x:=indsS[1];
y:=indsS[2];
z:=indsS[3];
A:=[[a,b],[c,d]];;

ITER_POLY_WARN:=bool;

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

F:=fn(A);
FF:=NullMat(Length(F),Length(F[1]));
for i in [1..Length(F)] do
for j in [1..Length(F[1])] do
if F[i][j]=0 then FF[i][j]:=[[],0];
else
pos:=Position([a^n,b^n,c^n,d^n],F[i][j]);
if IsBound(F[i][j]!.IndeterminateNumberOfUnivariateRationalFunction) then
FF[i][j] := [ [pos,n],1 ];

else 
FF[i][j]:=F[i][j]!.ExtRepPolynomialRatFun;

fi;
fi;
od;
od;

lnFF:=NullMat(Length(F),Length(F[1]));
for i in [1..Length(F)] do
for j in [1..Length(F[1])] do
lnFF[i][j]:=[1..Length(FF[i][j])/2];
od;
od;

lnFFk:=NullMat(Length(F),Length(F[1]));
for i in [1..Length(F)] do
for j in [1..Length(F)] do
lnFFk[i][j]:=[];
for k in lnFF[i][j] do
lnFFk[i][j][k]:=[1..Length(FF[i][j][2*k-1])/2];
od;
od;
od;


BB:=NullMat(Length(F),Length(F[1]));
ln:=Length(BB);
In:=[1..ln];
Elts:=[];
EltsSS:=[];
recfnn:=[];

####################################
fnn:=function(B)
local L, C, i, j, k, m, cf,aa,2k1,2k,pos;

if not B in EltsSS then AddSet(EltsSS,B); Add(Elts,B);pos:=Length(Elts); 

C:=0*BB;
L:=[B[1][1],B[1][2],B[2][1],B[2][2]];

for i in In do
for j in In do
cf:=0;

#for k in [1..Length(FF[i][j])/2] do
for k in lnFF[i][j] do
#2k:=2*k;
#2k1:=2*k-1;
2k:=2m[k];
2k1:=2m1[k];
aa:=1;
#for m in [1..Length(FF[i][j][2*k-1])/2] do
for m in lnFFk[i][j][k] do
aa:=aa*L[  FF[i][j][2k1][2m1[m]]   ]^FF[i][j][2k1][2m[m]]  ;
od;
aa:=aa*FF[i][j][2k];
cf:=cf+aa;
od;

C[i][j]:=cf;
od;
od;

recfnn[pos]:=C;
else
pos:=Position(Elts,B);
fi;

return recfnn[pos];
end;
####################################


M:=Group( List(GeneratorsOfGroup(G),fnn) );;
module:=GroupHomomorphismByFunction(G,M,x->fnn(x));

return module;
end);
########################################################
########################################################
