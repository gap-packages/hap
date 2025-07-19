#Throughout we try to denote Sqrt(-d) by S

############################################################
############################################################
InstallGlobalFunction(HAP_BianchiFundamentalRectangle, #####
function(OQ)                                           #####
local s1, s2, t1, t2;                       
#The fundamental rectangle is [s1,t1*S]x[s2,t2*S].
#The numbers s1,t1,s2,t2 are returned. We also use the 
#subregion Q=[0,t1*S]x[0,t2*S].

if not -OQ!.bianchiInteger mod 4 = 3  then
s1:=-1/2;t1:=1/2;s2:=-1/2;t2:=1/2;
else
s1:=-1/2;t1:=1/2;s2:=-1/4;t2:=1/4;
fi;
return [s1,t1,s2,t2];                                  #####
end);                                                  #####
############################################################ 
############################################################


############################################################
############################################################
InstallGlobalFunction(UnimodularPairStandardForm,      #####
function(x)                                            #####
local y;                                    
#The unimodular pairs [a,b] and [-a,-b] are equivalent.
#One of these pairs is returned.
y:=1*x;
y[1]:=-y[1];
y[2]:=-y[2];
return Maximum(x,y);;                                  #####
end);                                                  #####
############################################################
############################################################


####################################################
####################################################
InstallGlobalFunction(QuadraticIntegersByNorm, #####
function(OQ,N)      
local BI, NW, L, x,y,z;                      
#Returns all quadratic integers a with 
#-N <= Norm(OQ,a) <= N where OQ is a ring of integers

if N<0 then return []; fi;

BI:=OQ!.bianchiInteger;

if IsBound(OQ!.normomega) then NW:=OQ!.normomega; 
else
    if -BI mod 4 =3 then NW:=(1+AbsInt(BI))/4; OQ!.normomega:=NW;
    else NW:=AbsInt(BI); OQ!.normomega:=NW;
    fi;
fi;

L:=[]; #To be a list of integers of appropriate norm.

if -BI mod 4 =3 then
for x in [-(1+Int(Sqrt(1.0*N)))..1+Int(Sqrt(1.0*N))] do
for y in [-(1+Int(Sqrt(1.0*N/NW)))..1+Int(Sqrt(1.0*N/NW))] do
z:=QuadraticNumber((2*x+y)/2,y/2,BI);  
if HAPNorm(OQ,z)<=N then Add(L,z); fi;
od;
od;
else
for x in [-(1+Int(Sqrt(1.0*N)))..1+Int(Sqrt(1.0*N))] do
for y in [-(1+Int(Sqrt(1.0*N/NW)))..1+Int(Sqrt(1.0*N/NW))] do
z:=QuadraticNumber(x,y,BI);
if HAPNorm(OQ,z)<=N then Add(L,z); fi;
od;
od;
fi;

SortBy(L,x->HAPNorm(OQ,x));

return L;                                    #####
end);                                        #####
##################################################
##################################################


########################################################
########################################################
InstallGlobalFunction(HAP_UnimodularComplements,   #####
function(OQ,a)                                     #####
local Q, nrma,BI, rat, irrat, N, L, U, b, x,s1, t1, s2, t2;
#Returns the list of unimodular pairs (a,b) in OQxOQ that 
#intersect with the fundamental rectangle Q. The pair (a,b) 
#is returned as [a,b,Norm(a)]. Here we think of (a,b) as 
#being represented by a circular disk with centre b/a and radius |a| 
#in the complex plane.


BI:=OQ!.bianchiInteger;

Q:=HAPQuadratic(a);
nrma:=HAPNorm(OQ,a);

L:=QuadraticIntegersByNorm(OQ,(1-BI)*nrma);   

L:=Filtered(L,x->not x=0);
U:=[];
for b in L do
if IsQQUnimodularPair(OQ,[a,b,nrma]) then

if Norm(QuadraticIdeal(OQ,[QuadraticToCyclotomic(a),QuadraticToCyclotomic(b)])) = 1 
then Add(U,UnimodularPairStandardForm([a,b,nrma])); fi;

fi;

#fi;
od;
#Add(U,[1,0,1]);      #Added July 2025
return SSortedList(U);                       #####
end);                                        #####
##################################################
##################################################


#######################################################
#######################################################
InstallGlobalFunction(HAP_SqrtInequality,         #####
function(a,b,c);                                  #####
#Inputs three positive rational numbers and       
#returns Sqrt(a) + Sqrt(b) <= Sqrt(c).
if a+b>c then return false; fi;
return  4*a*b <= (c-a-b)^2;                       #####
end);                                             #####
#######################################################
#######################################################


#######################################################
#######################################################
InstallGlobalFunction(HAP_SqrtStrictInequality,   #####
function(a,b,c);                                  #####
#Inputs three positive rational numbers and      
#returns Sqrt(a) + Sqrt(b) <= Sqrt(c).
if a+b>=c then return false; fi;
return  4*a*b < (c-a-b)^2;                        #####
end);                                             #####
#######################################################
#######################################################


##########################################################
##########################################################
InstallGlobalFunction(HAP_IsRedundantUnimodularPair, #####
function(OQ,L,x)                                     #####
local y,l,m,lm,nrm,a,b,ba,rm,ra,nrmm,nrma;        
#Inputs a ring of integers OQ, a list L of unimodular 
#pairs in OQ, and a unimodular pair x. Returns true if 
#x is 'redundant' and false otherwise. We think of a
#unimodular pair [a,b] as a sphere of radius 1/|a| and centre
#b/a in the complex plane. So we assume not a=0. A sphere is
#redundant in this function if it lies underneath some other 
#sphere.

l:=x[2];
m:=x[1];
nrmm:=x[3];

for y in L do
a:=y[1];
b:=y[2];
nrma:=y[3];

nrm:=HAPNorm(OQ,l*a - b*m)/(nrmm*nrma);
if not nrm > 1/nrma then #so circle y might contain circle x
if HAP_SqrtInequality(nrm,1/nrmm,1/nrma) then return true; fi;
fi;
od;

return false;                                #####
end);                                        #####
##################################################
##################################################


######################################################
######################################################
InstallGlobalFunction(UnimodularPairCoordinates, #####
function(OQ,x)                                   #####
local a,b,Q,ax,ay,bx,by, D,X,Y,BI;          
#inputs a unimodular pair [a,b,1/r^2] and returns [x1,x2,r^2] 
#where r is the radius and [x1,x2*S] the centre of the 
#corresponding hemisphere. Here S:=Sqrt(BianchiInteger)

if Length(x)>3 then return [x[4],x[5],x[6]]; fi;

a:=x[2];
b:=x[1];
BI:=OQ!.bianchiInteger;

Q:=HAPQuadratic(2*a);
ax:=Q.a/(2*Q.d);
ay:=Q.b/(2*Q.d);
Q:=HAPQuadratic(2*b);
bx:=Q.a/(2*Q.d);
by:=Q.b/(2*Q.d);

D:=bx^2 - by^2*BI;
X:=ax*bx-ay*by*BI;
Y:=ay*bx-ax*by;

D:=D^-1;
x[4]:=X*D;
x[5]:=Y*D;
x[6]:=1/x[3];

return  [x[4],x[5],x[6]];
                                             #####
end);                                        #####
##################################################
##################################################


##################################################
##################################################
InstallGlobalFunction(UnimodularPairs,       #####
function(arg)                                #####
local OQ,N,bool,BOOL,L,LL,a,b,fn,U,UU,u,w,   
K,i,j,nrm,NRM;
#Returns a list of pairs [a,b] with Norm(a)<=N and (a,b)
#a unimodular pair intersecting the closed rectangle Q=
#[0,t1*S]x[0,t2*S]. With arg[3]=false no redundancies are 
#eliminated. With arg[4]=LL a list of unimodular pairs up 
#to a given norm -- some work is avoided.

OQ:=arg[1];
N:=arg[2];
if Length(arg)>2 then bool:=arg[3]; else bool:=true; fi;
#bool=true means that some obvious redundancies will be removed.
if Length(arg)>3 then LL:=List(arg[4],x->1*x); else LL:=[]; fi;
#LL is a list of all unimodular pairs of a neighbourhood of the rectangle 
#up to a certain norm.
if Length(arg)>4 then BOOL:=arg[5]; else BOOL:=false; fi;
#BOOL=true means a slightly larger neighbourhood of the fundamental recangle is considered.


L:=QuadraticIntegersByNorm(OQ,N);
L:=Filtered(L,x->not x=0);
if Length(LL)>0 then
NRM:=Maximum(List(LL,x->x[3]));
L:=Filtered(L,x->HAPNorm(OQ,x)>NRM);
else NRM:=0;
fi;

for a in L do
U:=HAP_UnimodularComplements(OQ,a);

if bool then
UU:=[];
for u in U do
  if not HAP_IsRedundantUnimodularPair(OQ,LL,u) then
  Add(UU,u);
  fi;
od;
Append(LL,UU);
else
Append(LL,U);
fi;
od;

#LL:=SSortedList(LL);
LL:=Classify(LL,x->x);;
LL:=List(LL,x->x[1]);
UU:= Filtered(LL,x->x[3]>NRM);
Append(LL,UU);

if bool then
for i in Reversed([2..Length(LL)]) do
nrm:=LL[i][3];
for j in Reversed([1..i-1]) do
if LL[j][3] < nrm then
if HAP_IsRedundantUnimodularPair(OQ,LL{[1..j]},LL[i]) then
LL[i]:=0;
fi;
break;
fi;
od;
od;
fi;

LL:=Filtered(LL,x->not x=0);
LL:=List(LL,x->[x[1],x[2],HAPNorm(OQ,x[1])]);
#LL:=SSortedList(LL);
LL:=Classify(LL,x->x);;
LL:=List(LL,x->x[1]);

SortBy(LL,x->x[3]);
#List(LL,x->UnimodularPairCoordinates(OQ,x));
return LL;                                   #####
end);                                        #####
##################################################
##################################################


#######################################################
#######################################################
InstallGlobalFunction(UnimodularIntersectingLine, #####
function(OQ,a,b)                                  #####
local ca,cb,u,v,w,S,ra2,rb2,u2,lambda;
#Inputs two unimodular pairs a, b representing circles
#in the complex plane. It is assumed that the circles 
#intersect. A pair [v,w] is returned where w is the 
#direction vector of the radical axis and v is the point 
#where the radical acis intersects the line joining 
#the centres of a and b. 

S:=QuadraticNumber(0,1,-OQ!.bianchiInteger);

ca:=UnimodularPairCoordinates(OQ,a);
ca[2]:=ca[2]*S;
cb:=UnimodularPairCoordinates(OQ,b);
cb[2]:=cb[2]*S;
u:=cb-ca;
u:=u{[1,2]};
w:=[u[2],-u[1]];

if u=[0,0] then return [ca{[1,2]},w]; fi;

ra2:=ca[3];
rb2:=cb[3];
#u2:=Norm(OQ,u[1])+Norm(OQ,u[2]);
u2:=u[1]^2+u[2]^2;

ca:=ca{[1,2]};
lambda:=(ra2-rb2+u2)/(2*u2);

v:=ca+lambda*u;
return [v,w];                                     #####
end);                                             #####
#######################################################
#######################################################


#########################################################
#########################################################
InstallGlobalFunction(UnimodularIntersectingPoint,  #####
function(OQ,a,b,c,Lab,Lac)                          #####
local v,w,vv,ww,D,lambda;
#Inputs three unimodular pairs a,b,c thought of as
#circles in the complex plane. It is assumed that the
#circles intersect pairwise. The radical centre v
#is returned.

#Lab:=UnimodularIntersectingLine(OQ,a,b);
#Lac:=UnimodularIntersectingLine(OQ,a,c);

v:=Lab[1];
w:=Lab[2];
vv:=Lac[1];
ww:=Lac[2];
D:=w[1]*(-ww[2])+w[2]*ww[1];
if D=0 then
return fail; fi;
D:=D^-1;
lambda:=D*( -ww[2]*(vv[1]-v[1]) + ww[1]*(vv[2]-v[2]) );

return v+lambda*w;                                #####
end);                                             #####
#######################################################
#######################################################


############################################################
############################################################
InstallGlobalFunction(HAP_HeightOfPointOnSphere,       #####
function(OQ,z,p)                                       #####
local c,S;                                          
#Inputs the real and imaginary coordinates [x,y] of a point
#z=x+i*y in the complex plain and unimodular pair p. 
#The pair p determines a hemisphere centred at some point 
#on the complex plane. The square h^2 of the height h of 
#the lift of z to the hemisphere is returned. If z does 
#not lie under the hemisphere then no lift exists and a 
#negative square height will be returned. 

S:=QuadraticNumber(0,1,-OQ!.bianchiInteger);

c:=UnimodularPairCoordinates(OQ,p);
S:=c[3]-(z[1]-c[1])^2-(z[2]-S*c[2])^2;
return S;                                             #####
end);                                                 #####
###########################################################
###########################################################


################################################################
################################################################
InstallGlobalFunction(HAP_AreIntersectingUnimodularPairs,  #####
function(OQ,a,b)                                           #####
local A, B, BI;                                        
#Inputs a ring of integers OQ of an imaginary quadratic field
#together with two unimodular pairs.
#It returns true if the two circles in the commplex plane
#corresponding to the unimodular pairs intersect.  

BI:=AbsInt(OQ!.bianchiInteger);

A:=UnimodularPairCoordinates(OQ,a);
B:=UnimodularPairCoordinates(OQ,b);
return not HAP_SqrtStrictInequality(A[3],B[3], (A[1]-B[1])^2 + BI*(A[2]-B[2])^2 );
                                                       #####
end);                                                  #####
############################################################
############################################################


################################################################
################################################################
InstallGlobalFunction(HAP_AreStrictlyIntersectingUnimodularPairs,  #####
function(OQ,a,b)                                           #####
local A, B, BI;
#Inputs a ring of integers OQ of an imaginary quadratic field
#together with two unimodular pairs.
#It returns true if the two circles in the commplex plane
#corresponding to the unimodular pairs intersect.

BI:=AbsInt(OQ!.bianchiInteger);

A:=UnimodularPairCoordinates(OQ,a);
B:=UnimodularPairCoordinates(OQ,b);
return not HAP_SqrtInequality(A[3],B[3], (A[1]-B[1])^2 + BI*(A[2]-B[2])^2 );
                                                       #####
end);                                                  #####
############################################################
############################################################



#################################################################
#################################################################
InstallGlobalFunction(HAP_Are3IntersectingUnimodularPairs,  #####
function(arg)                                               #####
local OQ, a,b,c, Lab, Lac, Lbc, A, B, C, v, w;                
#Inputs a ring of integers OQ of an imaginary quadratic field
#together with three unimodular pairs that intersect pairwise.
#It returns true if the three circles in the commplex plane
#corresponding to the unimodular pairs have common region 
#of overlap. The common intersections may be a point. 
#Otherwise false is returned.

OQ:=arg[1];
a:=arg[2];
b:=arg[3];
c:=arg[4];
Lab:=arg[5];
Lac:=arg[6];
Lbc:=arg[7];

###########################
A:=UnimodularPairCoordinates(OQ,a);
B:=UnimodularPairCoordinates(OQ,b);
C:=UnimodularPairCoordinates(OQ,c);

if HAP_HeightOfPointOnSphere(OQ,Lab[1],c)>=0 then return [true]; fi;

if HAP_HeightOfPointOnSphere(OQ,Lac[1],b)>=0 then return [true]; fi;

v:=Lab[2]{[2,1]};
v[1]:=-v[1];
if v*Lac[2]=0 then return [true]; fi;    #circles are colinear

#v:=UnimodularIntersectingLine(OQ,c,b);
if HAP_HeightOfPointOnSphere(OQ,Lbc[1],a)>=0 then return [true]; fi;

v:=UnimodularIntersectingPoint(OQ,a,b,c,Lab,Lac);

if HAP_HeightOfPointOnSphere(OQ,v,a)<0 then return [false]; fi;
if HAP_HeightOfPointOnSphere(OQ,v,b)<0 then return [false]; fi;
if HAP_HeightOfPointOnSphere(OQ,v,c)<0 then return [false]; fi;

return [true];                               #####
end);                                        #####
##################################################
##################################################


##################################################
##################################################
InstallGlobalFunction(HAP_PrintFloat,        #####
function(xx)                                 #####
local x,s,m, p, l, sgn;                      
x:=xx;
sgn:=SignFloat(x);
x:=AbsoluteValue(x);

s:=SplitString(String(x),'e');;
if Length(s)=1 then
s:=SplitString(String(xx),'e');;
return s[1]{[1..Minimum(5,Length(s[1]))]}; fi;

p:=EvalString(s[2]);;
l:=Length(s[1]);
m:=s[1]{[2..l]};
if p<0 then
s:=Concatenation(ListWithIdenticalEntries(-p,"0"));
s:=Concatenation(".",s,m{[1..l-1]});
s:=s{[1..5]};
if sgn<0 then s:=Concatenation("-",s); fi;
return s;
fi;
m:=EvalString(m);
s:=10^p*m;
s:=Concatenation(String(s),"00000");
s:=Concatenation(s{[1..p]},".",s{[p+1..Minimum(l,p+5)]});
if sgn<0 then s:=Concatenation("-",s); fi;
return s;                                    #####
end);                                        #####
##################################################
##################################################


###################################################
###################################################
InstallGlobalFunction(DisplayUnimodularPairs, #####
function(arg)                                #####
local OQ,L, C, file, i, x, w, BI,ABI, S1,S,SS,s1,t1,s2,t2;

OQ:=arg[1];
L:=arg[2];
if Length(arg)=3 then C:=arg[3]; 
else
C:=[];
fi;


BI:=OQ!.bianchiInteger;
ABI:=AbsInt(BI);
S:=Sqrt(1.0*ABI);
S1:=Sqrt(1.0*ABI);
if -BI mod 4 =3 then S:=0.25*S; else S:=0.5*S; fi;
SS:=HAP_PrintFloat(S);
SS:=SS{[1..Minimum(Length(SS),5)]};

C:=List(C,x->[x[1], x[2]*(1/QuadraticNumber(0,1,ABI))]);

C:=List(C,x->[x[1]!.rational,x[2]!.rational]);

C:=List(C,x->[1.0*x[1],S1*x[2]]);

file:="tmpp.asy";
PrintTo(file,"size(0,200);\n");
AppendTo(file,"defaultpen(0.2);\n");
AppendTo(file,"pen colour=blue+opacity(0.50);\n");
AppendTo(file,"real S=");
AppendTo(file,SS);
AppendTo(file,";\n");

s1:=HAP_PrintFloat(-0.5);t1:=HAP_PrintFloat(0.5);s2:=HAP_PrintFloat(-S);t2:=HAP_PrintFloat(S);

AppendTo(file,"path g=(-0.5,-S)--(0.5,-S)--(0.5,S)--(-0.5,S)--cycle;\n");

AppendTo(file,"picture fd;\n");

S:=Sqrt(1.0*ABI);
for x in L do
w:=UnimodularPairCoordinates(OQ,x);
w:=1.0*w;
AppendTo(file,"real r=");
AppendTo(file,"sqrt(");
AppendTo(file,HAP_PrintFloat(w[3]));
AppendTo(file,");\n");
AppendTo(file,"pair z=(");
AppendTo(file,HAP_PrintFloat(w[1]));
AppendTo(file,",");
AppendTo(file,HAP_PrintFloat(w[2]*S));
AppendTo(file,");\n");
AppendTo(file,"path c=circle(z,r);\n");
AppendTo(file,"fill(fd,c,colour);");
AppendTo(file,"draw(fd,c);\n");
od;

#AppendTo(file,"clip(fd,g);\n");
AppendTo(file,"draw(fd,g,white+linewidth(.3));\n");

for x in C do
AppendTo(file,"real r=0.01");
AppendTo(file,";\n");
AppendTo(file,"pair z=(");
AppendTo(file,HAP_PrintFloat(x[1]));
AppendTo(file,",");
AppendTo(file,HAP_PrintFloat(x[2]));
AppendTo(file,");\n");
AppendTo(file,"path c=circle(z,r);\n");
AppendTo(file,"fill(fd,c,red);");
AppendTo(file,"draw(fd,c);\n");
od;

AppendTo(file,"add(fd);\n");


Exec("asy -f pdf tmpp.asy; evince tmpp.pdf");     #####
end);                                              #####
#######################################################
#######################################################



#####################################################
#####################################################
InstallGlobalFunction(Display3DUnimodularPairs, #####
function(OQ,LL)                                 #####
local  file, x, w, BI,ABI, S,SS, L,R,T,B,sz,V;

V:=HAP_VertexHeights(OQ,LL)[2];
V:=V!.Heights;

V:=Filtered(V,x->x>0);

V:=List(V,x->x!.rational);

V:=1.0*Minimum(V);
V:=Sqrt(V);
V:=Minimum(0.05,V/2);


BI:=OQ!.bianchiInteger;
ABI:=AbsInt(BI);
S:=Sqrt(1.0*ABI);
if -BI mod 4 = 3 then S:=0.25*S; else S:=0.5*S; fi;
SS:=HAP_PrintFloat(S);
SS:=SS{[1..Minimum(Length(SS),5)]};
sz:=SizeScreen();
SizeScreen([200,24]);
file:="tmpp.asy";
PrintTo(file,"import solids;\n import graph3;\n import three;\n size(1000);\n");
AppendTo(file,"currentprojection=orthographic(5,4,10);\n");
AppendTo(file,"currentlight=Headlamp;\n");
AppendTo(file,"nslice=4*nslice;\n");
AppendTo(file,"revolution b=sphere(O,1);\n");

L:=[];R:=[];T:=[];B:=[];
S:=Sqrt(1.0*ABI);
for x in LL do
w:=UnimodularPairCoordinates(OQ,x);
w:=1.0*w;
w[2]:=w[2]*S;
Add(L,w[1]-Sqrt(w[3]));
Add(R,w[1]+Sqrt(w[3]));
Add(T,w[2]+Sqrt(w[3]));
Add(B,w[2]-Sqrt(w[3]));
AppendTo(file,"real r=");
AppendTo(file,"sqrt(");
AppendTo(file,HAP_PrintFloat(w[3]));
AppendTo(file,");\n");
AppendTo(file,"draw(shift(",HAP_PrintFloat(w[1]),",",HAP_PrintFloat(w[2]),",0)*scale3(r)*surface(b),blue+opacity(1));\n");
od;


AppendTo(file,"real S=");
AppendTo(file,SS);
AppendTo(file,";\n");
AppendTo(file,"real V=");
AppendTo(file,V);
AppendTo(file,";\n");

#if -BI mod 4 = 3 then
AppendTo(file,"path3 g=(-0.5,-S,1)--(0.5,-S,1)--(0.5,S,1)--(-0.5,S,1)--cycle;\n");
#else
#AppendTo(file,"path3 g=(0,0,1)--(1,0,1)--(1,S,1)--(0,S,1)--cycle;\n");
#AppendTo(file,"path3 gg=(0,0,V)--(1,0,V)--(1,S,V)--(0,S,V)--cycle;\n");
#fi;

AppendTo(file,"draw(g,white+linewidth(0.3));");

AppendTo(file,"real L=");
AppendTo(file,Minimum(L));
AppendTo(file,";\n");

AppendTo(file,"real R=");
AppendTo(file,Maximum(R));
AppendTo(file,";\n");

AppendTo(file,"real T=");
AppendTo(file,Maximum(T));
AppendTo(file,";\n");

AppendTo(file,"real B=");
AppendTo(file,Minimum(B));
AppendTo(file,";\n");


AppendTo(file,"triple p0 =(L,T,V), p1=(R,T,V), p2 = (R,B,V), p3 = (L,B,V);\n");
AppendTo(file,"triple q0 =(L,T,-1), q1=(R,T,-1), q2 = (R,B,-1), q3 = (L,B,-1);\n");

AppendTo(file,"skeleton s;\n");
AppendTo(file,"b.transverse(s,reltime(b.g,0.5),P=currentprojection);\n");
AppendTo(file,"draw(surface(p0--p1--p2--p3--cycle),white);\n");
AppendTo(file,"draw(surface(q0--q1--q2--q3--cycle),white);\n");
AppendTo(file,"draw(surface(p0--p1--q1--q0--cycle),white);\n");
AppendTo(file,"draw(surface(q1--p1--p2--q2--cycle),white);\n");
AppendTo(file,"draw(surface(q3--q2--p2--p3--cycle),white);\n");
AppendTo(file,"draw(surface(p0--q0--q3--p3--cycle),white);\n");

SizeScreen(sz);

Exec("asy -V tmpp.asy");                           #####
end);                                              #####
########################################################
########################################################


###################################################
###################################################
InstallGlobalFunction(IsQUnimodularPair,     #####
function(arg)                                 #####
local OQ,x,BI,ABI,L,s1,t1,s2,t2,T,S;                
#Inputs a unimodular pair x=[a,b] and returns true
#if the associated hemisphere touches the 
#fundamental rectangle [s1,t1*S]x[s2,t2*S] in the 
#complex plane.


OQ:=arg[1];
BI:=OQ!.bianchiInteger;
ABI:=AbsInt(BI);
x:=arg[2];
L:=HAP_BianchiFundamentalRectangle(OQ);
s1:=L[1];t1:=L[2];; s2:=L[3]; t2:=L[4];

T:=UnimodularPairCoordinates(OQ,x);

if T[1]>=0 and
(not  HAP_SqrtStrictInequality(t1^2,T[3],T[1]^2) ) and
 T[2]>=s2 and T[2]<=t2
then return true; fi;

if T[2]>=0 and 
(not HAP_SqrtStrictInequality(ABI*t2^2,T[3],ABI*T[2]^2) ) and
 T[1]>=s1 and T[1]<=t1
then return true; fi;


if T[1]<=0 and #=
(not  HAP_SqrtStrictInequality(s1^2,T[3],T[1]^2)  ) and
 T[2]>=s2 and T[2]<=t2
then return true; fi;


if T[2]<=0 and   #=
(not  HAP_SqrtStrictInequality(ABI*t2^2,T[3],ABI*T[2]^2)  ) and
 T[1]>=s1 and T[1]<=t1
then return true; fi;


if (T[1]-s1)^2 + ABI*(T[2]-s2)^2<=T[3]
then return true; fi;

if (T[1]-s1)^2 + ABI*(T[2]-t2)^2<=T[3]
then  return true; fi;

if (T[1]-t1)^2 + ABI*(T[2]-s2)^2<=T[3]
then  return true; fi;

if (T[1]-t2)^2 + ABI*(T[2]-t2)^2<=T[3]
then  return true; fi;

return false;
                                             #####
end);                                        #####
##################################################
##################################################


###################################################
###################################################
InstallGlobalFunction(IsQQUnimodularPair,     #####
function(arg)                                 #####
local OQ,x,BI,ABI,L,s1,t1,s2,t2,T,S;
#Inputs a unimodular pair x=[a,b] and returns true
#if the associated hemisphere touches the 
#fundamental Q-rectangle [0,t1]x[0,t2] in the complex 
#plane.


OQ:=arg[1];
BI:=OQ!.bianchiInteger;
ABI:=AbsInt(BI);
x:=arg[2];
L:=HAP_BianchiFundamentalRectangle(OQ);
s1:=0;t1:=L[2];; s2:=0; t2:=L[4];
#s1:=L[1];t1:=L[2];; s2:=L[3]; t2:=L[4];

T:=UnimodularPairCoordinates(OQ,x);

if T[1]>=0 and   #=
(not  HAP_SqrtStrictInequality(t1^2,T[3],T[1]^2) )
and T[2]>=0 and T[2]<=t2
then return true; fi;

if T[2]>=0 then    #=
if (not HAP_SqrtStrictInequality(ABI*t2^2,T[3],ABI*T[2]^2) )
and T[1]>=0 and T[1]<=t1
then return true; fi;
fi;

if T[1]<=0 then #=
if T[1]^2<=T[3]
and T[2]>=0 and T[2]<=t2
then return true; fi;
fi;

if T[2]<=0 then  #=
if ABI*T[2]^2<=T[3]
and T[1]>=0 and T[1]<=t1
then return true; fi;
fi;

if (T[1]-s1)^2 + ABI*(T[2]-t1)^2<=T[3]
then return true; fi;

if (T[1]-s1)^2 + ABI*(T[2]-t2)^2<=T[3]
then return true; fi;

if (T[1]-s2)^2 + ABI*(T[2]-t1)^2<=T[3]
then return true; fi;

if (T[1]-s2)^2 + ABI*(T[2]-t2)^2<=T[3]
then return true; fi;

return false;
                                             #####
end);                                        #####
##################################################
##################################################


################################################################
################################################################
InstallGlobalFunction(SimplicialComplexOfUnimodularPairs,  #####
function(OQ,L)                                             #####
local E,S,K,Lines,Points,bool,p,i,j,k,x;                  
#Inputs a list of unimodular pairs corresponding to disks 
#in the complex plane. The nerve of these disks is
#returned as an abstract simplicial complex.

Lines:=[];
K:=List(L,x->[Position(L,x)]);;
E:=[];
for i in [1..Length(L)] do
Lines[i]:=[];
for j in [i+1..Length(L)] do
if HAP_AreIntersectingUnimodularPairs(OQ,L[i],L[j])
then
Add(E,[i,j]);
fi;
od;od;

E:=SSortedList(E);

Append(K,E);

Points:=[];
for x in E do
i:=x[1]; j:=x[2];
for k in [j+1..Length(L)] do
if [i,k] in E and [j,k] in E then
  if not IsBound(Lines[i][j]) then
     Lines[i][j]:=UnimodularIntersectingLine(OQ,L[i],L[j]);
  fi;
  if not IsBound(Lines[i][k]) then
     Lines[i][k]:=UnimodularIntersectingLine(OQ,L[i],L[k]);
  fi;
  if not IsBound(Lines[j][k]) then
     Lines[j][k]:=UnimodularIntersectingLine(OQ,L[j],L[k]);
  fi;
  p:= HAP_Are3IntersectingUnimodularPairs(OQ,L[i],L[j],L[k],
                                  Lines[i][j],Lines[i][k],Lines[j][k]);
  if p[1] then
     Add(K,[i,j,k]);
     if not IsBound(p[2]) then
         p[2]:=UnimodularIntersectingPoint(OQ,
                  L[i],L[j],L[k],Lines[i][j],Lines[i][k]);
     fi;
     Add(Points,p[2]);
  fi;
fi;
od;od;


K:=MaximalSimplicesToSimplicialComplex(K);
K!.Lines:=Lines;
K!.Points:=Points;
K!.unimodularPairs:=L;
return K;                                      #####
end);                                          #####
####################################################
####################################################


#####################################################
#####################################################
InstallGlobalFunction(QQNeighbourhoodOfUnimodularPairs,
function(arg)                                   #####
local OQ,L,U,F,FF, Q, BI, x, p, u, S,St;      #####
#Inputs a list of unimodular pairs and returns a list 
#of all pairs/hemispheres in L that intersect
#the rectangle Q=[0,t1]x[0,t2]. 
OQ:=arg[1];
L:=arg[2];
F:=[];
Q:= Filtered(L,x->IsQQUnimodularPair(OQ,x));   
Q:=List(Q,x->UnimodularPairStandardForm(x));
BI:=OQ!.bianchiInteger;
if -BI mod 4 = 3 then
St:=[QuadraticNumber(1/2,1/2,BI),QuadraticNumber(-1/2,1/2,BI)];
else
St:=[QuadraticNumber(0,1,BI)];
fi;

for x in Q do
for S in St do
Add(F,x);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2]),x[3]]);
Add(F,[x[1],-x[2],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2]),x[3]]);
##############
Add(F,[x[1],x[2]-x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])+QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])-QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+x[1],x[3]]);
##############
Add(F,[x[1],x[2]-S*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])-S*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],x[2]-(1+S)*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])+(1-S)*QuadraticNumberConjugate(x[1]),x[3]]);
##############
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])+S*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+S*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])+(-1+S)*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+(1+S)*x[1],x[3]]);
od;
od;

F:=Classify(F,x->UnimodularPairStandardForm(x){[1,2]});;
F:=List(F,x->x[1]);
F:=List(F,x->x{[1,2,3]});

U:=Filtered(F,x->IsQQUnimodularPair(OQ,x));     
FF:=[];
for x in F do
for u in U do
if HAP_AreIntersectingUnimodularPairs(OQ,x,u) then
Add(FF,x); break; fi;
od;od;

return FF;                                       #####
end);                                           #####
#####################################################
#####################################################


#####################################################
#####################################################
InstallGlobalFunction(QNeighbourhoodOfUnimodularPairs,
function(arg)                                   #####
local OQ,L,U,F,FF, Q, BI, x, p, u, S,St;      #####
#Inputs a list of unimodular pairs and returns a list 
#of the 'images' of pairs in L that intersect the  
#fundamental rectangle [s1,t1]x[s2,t2]. 
OQ:=arg[1];
L:=arg[2];
F:=[];
Q:= Filtered(L,x->IsQQUnimodularPair(OQ,x));   
Q:=List(Q,x->UnimodularPairStandardForm(x));
BI:=OQ!.bianchiInteger;
if -BI mod 4 = 3 then
St:=[QuadraticNumber(1/2,1/2,BI),QuadraticNumber(-1/2,1/2,BI)];
else
St:=[QuadraticNumber(0,1,BI)];
fi;

for x in Q do
for S in St do
Add(F,x);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2]),x[3]]);
Add(F,[x[1],-x[2],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2]),x[3]]);
##############
Add(F,[x[1],x[2]-x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])+QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])-QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+x[1],x[3]]);
##############
Add(F,[x[1],x[2]-S*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])-S*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],x[2]-(1+S)*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),-QuadraticNumberConjugate(x[2])+(1-S)*QuadraticNumberConjugate(x[1]),x[3]]);
##############
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])+S*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+S*x[1],x[3]]);
Add(F,[QuadraticNumberConjugate(x[1]),QuadraticNumberConjugate(x[2])+(-1+S)*QuadraticNumberConjugate(x[1]),x[3]]);
Add(F,[x[1],-x[2]+(1+S)*x[1],x[3]]);
od;
od;

F:=Classify(F,x->UnimodularPairStandardForm(x){[1,2]});;
F:=List(F,x->x[1]);
F:=List(F,x->x{[1,2,3]});
#return F;    

U:=Filtered(F,x->IsQUnimodularPair(OQ,x));

return U;                                       #####
end);                                           #####
#####################################################
#####################################################



#####################################################
#####################################################
InstallGlobalFunction(HAP_VertexHeights,        #####
function(OQ,LL)                                 #####
local L,Lines, Points, Heights, Cusps, K, 
 s, c, D,i,bool,S,ABI,a,b;

ABI:=AbsInt(OQ!.bianchiInteger);
L:=QQNeighbourhoodOfUnimodularPairs(OQ,LL);    
K:=SimplicialComplexOfUnimodularPairs(OQ,L);  
Lines:=K!.Lines;
Points:=K!.Points;
S:=K!.simplicesLst[3];

Heights:=[];
for i in [1..Length(S)]  do
if Points[i] = fail then Heights[i]:=-infinity; else
Heights[i]:= HAP_HeightOfPointOnSphere(OQ,Points[i],L[S[i][1]]);
fi;
od;
K!.Heights:=Heights;

Cusps:=[];
for i in [1..Length(Points)] do
if Heights[i]=0 then
bool:=true;
for s in L do 
c:=UnimodularPairCoordinates(OQ,s);
b:=Points[i][2]*(1/QuadraticNumber(0,1,ABI));
a:=(c[1]-Points[i][1])^2 + ABI*(c[2]-b)^2;
if a < c[3] then
bool:=false; break; fi;
od;
if bool then Add(Cusps, Points[i]); fi;
fi;
od;

Cusps:=SSortedList(Cusps);
K!.Cusps:=Cusps;

return [Cusps,K];                               #####
end);                                           #####
#####################################################
#####################################################


###################################################
###################################################
InstallGlobalFunction(UnimodularPairsReduced_NN, #####
function(OQ,LL)                                  #####
local L, Lines, Points, K, K3, E,EE, V, VV,
fn, fn2, v, w, h, dist, bool, BOOL, C,D,DD,i,F,S,M;

#L:=QQNeighbourhoodOfUnimodularPairs(OQ,LL);
if LL=[] then return [LL,[]]; fi;
L:=LL;

K:=SimplicialComplexOfUnimodularPairs(OQ,L);
Lines:=K!.Lines;
Points:=K!.Points;

K3:=K!.simplicesLst[3];

DD:=Filtered([1..K!.nrSimplices(2)], v-> IsQQUnimodularPair(OQ,L[K3[v][1]])
or IsQQUnimodularPair(OQ,L[K3[v][2]])
or IsQQUnimodularPair(OQ,L[K3[v][3]]) );
#Need to check that this is OK here

S:=K3{DD};
Points:=Points{DD};

V:=List([1..Length(S)],x->[S[x][1],S[x][2],S[x][3], Points[x] ] );

V:=Filtered(V,v->not v[4]=fail);

V:=List(V,v-> [v[1],v[2],v[3],v[4],HAP_HeightOfPointOnSphere(OQ,v[4],L[v[1]])]);
V:=Filtered(V,v->v[5]>=0);

VV:=[];
for v in V do
#h:=HAP_HeightOfPointOnSphere(OQ,v[4],L[v[1]]) ;
h:=v[5];
bool:=true; ;
#NEXT LINE IS TIME CONSUMING
#DD:=Filtered(K!.simplicesLst[3],x-> (v[1] in x and v[2] in x) or (v[1] in x and v[3] in x) or (v[3] in x and v[2] in x));  #?????????
DD:=Filtered(K!.simplicesLst[3],x-> (v[1] in x and v[2] in x) );

DD:=Concatenation(DD);
DD:=SSortedList(DD);
for i in DD do
w:=L[i];
#for w in L do
if  HAP_HeightOfPointOnSphere(OQ,v[4],w) - h > 0  then
bool:=false; break; fi;
od;
if bool then Add(VV,v); fi;
od;



V:=List(VV,x->[SortedList([x[1],x[2],x[3]]),x[4],x[5]]);
D:=K!.vertices; ##??
DD:=[];
for v in D do
F:=Filtered([1..Length(V)],i-> v in V[i][1]);
F:=List(F,i->[V[i][2],V[i][3]]);
F:=SSortedList(F);
if Length(F)>2 then Add(DD,v);  fi;
od;

L:=L{DD};
L:=QQNeighbourhoodOfUnimodularPairs(OQ,L); 
L:=SSortedList(L);


VV:=Filtered(VV,v->v[1] in DD and v[2] in DD and v[3] in DD);
M:=List(VV,v->Position(K!.simplicesLst[3],[v[1],v[2],v[3]]));
M:=Difference([1..Length(K!.Points)], M);
for i in M do
K!.Points[i]:=fail;
#K!.Heights:=-infinity;
od;


return [L,K];                                    #####
end);                                            #####
##################################################
##################################################

####################################################
####################################################
InstallGlobalFunction(UnimodularPairsReduced,
function(OQ,LL)
local K,L;

L:=QNeighbourhoodOfUnimodularPairs(OQ,LL);
if OQ!.bianchiInteger=-5 then
L:=Filtered(L,x->IsQUnimodularPair(OQ,x));
else
L:=Filtered(L,x->IsQQUnimodularPair(OQ,x));
fi;
K:=UnimodularPairsReduced_NN(OQ,L);
K:=QQNeighbourhoodOfUnimodularPairs(OQ,K[1]);
K:=UnimodularPairsReduced_NN(OQ,K);
return K;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(CoverOfUnimodularPairs,  #####
function(arg)                                  #####
local OQ,N, NRMS, L, LL, K, Y, bool, BI, A,     
ABI, LST, UU, u, w, a, b, I, i, pos, HAPRECORD;

OQ:=arg[1];
N:=1;
BI:=OQ!.bianchiInteger;
ABI:=AbsInt(BI);
###############################

HAPRECORD:=[
[ 1, 1 ],[ 2, 1 ],[ 3, 1 ],[ 5, 20 ],[ 6, 24 ],[ 7, 1 ],[ 10, 40 ],[ 11, 1 ],[ 13, 52 ],[ 14, 56 ],[ 15, 15 ],[ 17, 68 ],[ 19, 5 ],[ 21, 84 ],[ 22, 88 ],[ 23, 18 ],[ 26, 121 ],[ 29, 169 ],[ 30, 169 ],[ 31, 20 ],[ 33, 528 ],[ 34, 225 ],[ 35, 35 ],[ 37, 289 ],[ 38, 4 ],[ 39, 39 ],[ 41, 361 ],[ 42, 672 ],[ 43, 9 ],
#[ 46, 25 ],
[ 47, 36 ],[ 51, 51 ],
#[ 53, 25 ],
[ 55, 56 ],
#[ 57, 49 ],
[ 59, 36 ],[ 67, 23 ],[ 71, 64 ],[ 79, 64 ],[ 83, 36 ],[ 87, 87 ],[ 91, 91 ],[ 95, 95 ],[ 103, 64 ],[ 107, 64 ],[ 111, 111 ],[ 115, 115 ],
#[ 119, 25 ],
[ 123, 123 ],[ 127, 121 ],[ 131, 100 ],[ 139, 100 ],[ 143, 576 ],[ 151, 110 ],[ 155, 155 ],[159,240],[ 163, 53 ],[167,384],[ 179, 144], [199,196]
];;
                 #HAPRECORD contains some precomuted [d,r] for which Swan's
                 #criterion is satisfied for the given norm r.
#HAPRECORD:=[];  #Uncomment this line to calculate HAPRECORD entries from 
                 #scratch.  

if Length(arg)=2 then pos:=arg[2];
else

pos:=Position(List(HAPRECORD,x->x[1]),ABI);
if pos =fail or pos=0 then pos:=infinity;
if Length(arg)=1 then
Print("Try \n  P:=BianchiPolyhedron(OQ,N);\nfor some guessed positive integer value of N and then try\n  SwanBianchiCriterion(P);\nto test if the value of N was large enough. If the test returns false then you'll need to try a larger value of N.\n\nA successful value of N can be stored as a pair [d,N] in the list HAPRECORD which can be edited manually in the file hap/lib/Congruence/bianchi.gi .\n\n");
fi;
else pos:=HAPRECORD[pos][2]; fi;

fi;
###############################
NRMS:=QuadraticIntegersByNorm(OQ,1000);  #We'll never need anything higher!!!
NRMS:=List(NRMS,x->HAPNorm(OQ,x));
NRMS:=SSortedList(NRMS);
bool:=true;

L:=[];

while bool do
N:=N+1; 
if Length(arg)>1 then
#Print("Adding hemispheres of squared radius ",1/NRMS[N],"\n");
fi;

L:=UnimodularPairs(OQ,NRMS[N],true,L);
if pos=infinity then
K:=SimplicialComplexOfUnimodularPairs(OQ,
                 QNeighbourhoodOfUnimodularPairs(OQ,L));
Y:=RegularCWComplex(K);
CocriticalCellsOfRegularCWComplex(Y,1);
if Homology(Y,0)=[0] and Homology(Y,1)=[] then break; fi;
fi;
if NRMS[N]>=pos then break; fi;
od;

K:=UnimodularPairsReduced(OQ,L);

L:=K[1]; K:=K[2];

###########CHECK SWAN'S CONDITION REALLY IS MET
if pos=infinity then 
A:=BianchiPolyhedron(OQ!.bianchiInteger,L);
A:=SwanBianchiCriterion_alt(A);
if not Length(A)=0 then 
return CoverOfUnimodularPairs(OQ,N+1,true);
fi;
fi;
##################################################

L:=Filtered(L,x->IsQQUnimodularPair(OQ,x));   
L:=QQNeighbourhoodOfUnimodularPairs(OQ,L);
L:=Filtered(L,x->IsQQUnimodularPair(OQ,x));
return [L,K];                                    #####
end);                                        #####
##################################################
##################################################

##################################################
##################################################
InstallGlobalFunction(IsUnimodularCollection,
function(OQ,L)
local K,Y;
K:=SimplicialComplexOfUnimodularPairs(OQ,L);
Y:=RegularCWComplex(K);
CocriticalCellsOfRegularCWComplex(Y,1);
if not Homology(Y,0)=[0] then return false; fi;
if not Homology(Y,1)=[] then return false; fi;
return true; 
end);
##################################################
###############################################

##################################################
##################################################
InstallGlobalFunction(SwanBianchiCriterion,
function(P)
local mr, mv, R, L,LL, i, u, Points;

mr:=P!.minRadius;
mv:=P!.minVertexHeight;
if IsInt(1/mv) then R:=1/mv;
else R:= 1+Int(1/mv);
fi;

L:=UnimodularPairs(P!.ring,R);;
LL:=UnimodularPairs(P!.ring,mr);;
L:=Difference(L,LL);;
L:=Filtered(L,x->IsQQUnimodularPair(P!.ring,x));;
L:=Filtered(L,x->not HAP_IsRedundantUnimodularPair(P!.ring,P!.unimodularPairs,x));

Points:=Filtered([1..Length(P!.points)],i->P!.visibleHeights[i]>0);

for u in L do
for i in Points do
if P!.visibleHeights[i] < HAP_HeightOfPointOnSphere(P!.ring,P!.points[i],u) then return false; fi;
od;
od;

return true;

end);
##################################################
##################################################


##################################################
##################################################
InstallGlobalFunction(SwanBianchiCriterion_alt,
function(arg)
local  PP,P,F,A,K,C,CC,D,H, nrmL, L, S, S1, v,b,x,i,h,s1,s2,t1,t2;

P:=arg[1];

if not IsUnimodularCollection(P!.ring,P!.unimodularPairs) then 
return [false]; fi;

H:=P!.heights;
H:=Filtered(H,x->x>0);
H:=1/Minimum(H);
if IsHapQuadraticNumber(H) then H:=H!.rational; fi;

L:=P!.unimodularPairs;
L:=List(L,x->HAPNorm(P!.ring,x[1]));;
nrmL:=Maximum(L);

if nrmL>=H then return []; fi;

L:=UnimodularPairs(P!.ring,nrmL);

if Length(arg)=2 then 
   H:=UnimodularPairs(P!.ring,arg[2]);
   H:=QNeighbourhoodOfUnimodularPairs(P!.ring,H);
else
   H:=UnimodularPairs(P!.ring,H);
   H:=Filtered(H,h->not h in L);
   H:=QNeighbourhoodOfUnimodularPairs(P!.ring,H);
fi;


S:=QuadraticNumber(0,1,-P!.bianchiInteger);
S1:=1/S;

L:=HAP_BianchiFundamentalRectangle(P!.ring);
s1:=0; t1:=L[2]; s2:=0; t2:=L[4];
#s1:=L[1];t1:=L[2];; s2:=L[3]; t2:=L[4];


D:=[];
for i in [1..Length(P!.points)] do
if not P!.points[i]=fail then
if P!.heights[i]>0 
and P!.points[i][1]>=s1 and P!.points[i][1]<=t1 and P!.points[i][2]*S1>=s2
and P!.points[i][2]*S1<=t2
then Add(D,i);
fi;
fi;
od;

C:=[]; CC:=[];
for i in D do
for h  in H do
if HAP_HeightOfPointOnSphere(P!.ring,P!.points[i],h) > P!.heights[i]
then 
Add(C,h);
Add(CC,i);
fi;
od;
od;
CC:=SSortedList(CC);

K:=P!.simplicialComplex;
A:=[];
for v in K!.vertices do
F:=Filtered(K!.simplicesLst[3],x-> v in x);
F:=Concatenation(F);
F:=SSortedList(F);
if IsSubset(CC,F) then Add(A,K!.unimodularPairs[v]); fi;
od;
if Length(A)>0 then return A; fi;

C:=Filtered(C,x->x[3]>nrmL);
PP:=BianchiPolyhedron(P!.bianchiInteger,Concatenation(P!.unimodularPairs,C));
return Filtered(PP!.unimodularPairs,x->not x in P!.unimodularPairs);

end);
##################################################
##################################################

#####################################################
#####################################################
InstallGlobalFunction(BianchiPolyhedron,
function(arg)
local d,K, P,C, OQ, F, H, L,M,i,x,N,s;

d:=arg[1];
if d>=0 then
Print("Input must be a negative square free integer.\n");
return fail;
fi;

F:=Factors(AbsInt(d));

if not Length(F)=Length(SSortedList(F)) then
Print("Input is not a square free integer.\n");
return fail;
fi;

P:=Objectify(HapBianchiPolyhedron,rec());
P!.bianchiInteger:=d;
OQ:=RingOfIntegers(QuadraticNumberField(d));
P!.ring:=OQ;

if AbsInt(d) in [1,3] then
Print("Warning: the polyhedron is a union of fundamental domains for PSL( ",OQ, ").\n");
fi; 

if Length(arg)=1 then
    K:=CoverOfUnimodularPairs(OQ); P!.unimodularPairs:=K[1]; K:=K[2];
else

    if IsInt(arg[2]) then
        K:=CoverOfUnimodularPairs(OQ,arg[2]); P!.unimodularPairs:=K[1]; K:=K[2];
    else
        K:=UnimodularPairsReduced(OQ,arg[2]); P!.unimodularPairs:=K[1]; K:=K[2];
    fi;
fi;

N:=4;
while Length(P!.unimodularPairs)=0 do 
#Print("Try again with a larger choice of limit on the norms or a larger list of unimodular pairs\n");
N:=N+1;     #Bit of a botch job!
K:=CoverOfUnimodularPairs(OQ,N); P!.unimodularPairs:=K[1]; K:=K[2];
od;


K:=HAP_VertexHeights(OQ,P!.unimodularPairs);

P!.unimodularPairs:=K[2]!.unimodularPairs;
P!.simplicialComplex:=K[2];
P!.vertices:=K[2]!.Points;
P!.heights:=K[2]!.Heights;
C:=[];
for x in K[1] do
Add(C,x);
Add(C,[-x[1],x[2]]);
Add(C,[x[1],-x[2]]);
Add(C,[-x[1],-x[2]]);
od;
C:=SSortedList(C);
P!.cusps:=C;




L:=K[2]!.unimodularPairs;
H:=P!.heights;
############
for i in [1..Length(H)] do
if H[i]>0 then
for s in L do 
if HAP_HeightOfPointOnSphere(OQ,K[2]!.Points[i],s)!.rational>H[i] then
H[i]:=-infinity;
fi;
od;
fi;
od;
############
P!.visibleHeights:=1*H;
H:=Filtered(H,x->x>0);
if Length(H)>0 then
H:=Minimum(H);
else
H:=0;
fi;

L:=P!.unimodularPairs;
L:=List(L,x->HAPNorm(P!.ring,x[1]));;
L:=1/Maximum(L);

P!.minVertexHeight:=H;
P!.minRadius:=L;
P!.points:=K[2]!.Points;

P!.unimodularPairs:=QNeighbourhoodOfUnimodularPairs(P!.ring,P!.unimodularPairs);
return P;

end);
####################################################
####################################################

