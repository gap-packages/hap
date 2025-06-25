
##############################################################
##############################################################
InstallGlobalFunction(IsHapQuadraticInteger,
function(OQ,x)
local a,b,d;

d:=OQ!.bianchiInteger;

a:=x!.rational;
b:=x!.irrational;
if d mod 4 = 1 then a:=a-b; b:=2*b;fi;

if IsInt(a) and IsInt(b) then return true;fi;
return false;

end);
##############################################################
##############################################################

##############################################################
##############################################################
InstallGlobalFunction(HAP_BianchiTransformations,
function(OQ,V,W)
local PairToQuadInt,
      HAPComplexConjugate,
      Potential_c,
      IntegerCircle,
      Potential_d,
      Potential_b,
      Standardize,
      IsIntMat,
      one,zero, 
      c,x,y,a,b,d,cx,D,G,z,zz,t,tt,cc,rt,rtt,pqv1,pqw1,
      pot_c,pot_dc,pot_abcd;

#returns the matrices in SL(2,OQ) that send V to W. 
#If V=W then these are returned as a group. Here V,W are lists of length 3.

if IsHapQuadraticNumber(V[1]) then V[1]:=V[1]!.rational; fi;
if IsHapQuadraticNumber(W[1]) then W[1]:=W[1]!.rational; fi;
z:=[V[1],V[2]]; t:=V[3];
zz:=[W[1],W[2]]; tt:=W[3];
if t*tt=0 and not t=tt then return []; fi;

#Let z:=[x,y]. Then [z,t] is viewed as a point in upper half-space, t>0. 
#So too is [zz,tt]. We sometimes write v:=[z,t] and w:=[zz,tt].
#We want to find all matrices M:=[[a,b],[c,d]] in SL(2,OQ) such 
#that M.[z,t] = [zz,tt].   
#Here x and y are rational numbers (IsRat has value true) and t is a cyclotomic real number.
##So too for zz and tt.

one:=QuadraticNumber(1,0,OQ!.bianchiInteger);
zero:=QuadraticNumber(0,0,OQ!.bianchiInteger);
##################################################
IsIntMat:=function(M);
if not IsHapQuadraticInteger(OQ,M[1]) then return false; fi;
if not IsHapQuadraticInteger(OQ,M[2]) then return false; fi;
if not IsHapQuadraticInteger(OQ,M[3]) then return false; fi;
if not IsHapQuadraticInteger(OQ,M[4]) then return false; fi;
return true;
end;
##################################################

##################################################
PairToQuadInt:=function(z);
if IsRat(z[2]) then
return z[1]+z[2]*HAPSqrt(OQ!.bianchiInteger);   #for case d=-1
else
return z[1]+z[2]!.irrational*HAPSqrt(OQ!.bianchiInteger);
fi;
end;
##################################################



##################################################
HAPComplexConjugate:=function(z);
if not IsHapQuadraticNumber(z) then return z; fi;   #case z is rational
return QuadraticNumber(z!.rational, -z!.irrational, z!.bianchiInteger);
end;
##################################################

##################################################
Potential_c:=function(v,w)
local ttt;

#ttt:=(t*tt)^2;
ttt:=t*tt;
if ttt=0 then ttt:=1/1000;
else
ttt:=Rat(Sqrt(Float(ttt)));
fi;

return QuadraticIntegersByNorm(OQ,1+1/(ttt));

end;
##################################################


##################################################
IntegerCircle:=function(OQ,c,r)
local C,cc,L,o,x,xx,a,b,BI,k;
#return all integers on the circle of radius=Sqrt(r) and centre c where c
#is a quadratic integer. 
if not IsRat(r) then return []; fi;
if r<0 then return []; fi;

BI:=-OQ!.bianchiInteger;
L:=Sqrt(1.0*r);
L:=2+Int(L);
o:=Int(c!.rational);
C:=[];


if not BI mod 4 = 3 then     
######################################
for x in [-L..L] do
xx:=x+o;
b:=r-(xx-c!.rational)^2;

k:=b/BI;
if not  (IsSquareInt(DenominatorRat(k)) and IsSquareInt(NumeratorRat(k)))   then continue; fi;
b:=c!.irrational+Sqrt(k);

if IsInt(b) then
Add(C,QuadraticNumber(xx,b,-BI));
fi;
od;
return C;
######################################
else
######################################
for x in [-2*L..2*L] do
xx:=(x/2)+o;
b:=r-(xx-c!.rational)^2;
     if IsSquareInt(NumeratorRat(b/BI)) and IsSquareInt(DenominatorRat(b/BI)) 
     then   ##############  I guess this is mathematically redundant!
b:=c!.irrational+Sqrt(b/BI);
if IsInt(xx-b) and IsInt(2*b) then
Add(C,QuadraticNumber(xx,b,-BI)); 
fi;
     fi;    ######################

od;
return C;
######################################
fi;






#STUPID ALTERNATIVE IMPLEMENTATION BELOW!  

cc:=QuadraticNumber(Int(c!.rational),Int(c!.irrational),c!.bianchiInteger);

#return all integers on the circle of radius Sqrt(r) and centre c.
C:=QuadraticIntegersByNorm(OQ,r-OQ!.bianchiInteger+1);;
C:=Difference(C,QuadraticIntegersByNorm(OQ,r+OQ!.bianchiInteger-1));
C:=C+cc;
C:=Filtered(C,x->HAPNorm(OQ,-c+x)=r);
return C;

end;
##################################################


#if tt=0 then rtt:=infinity; else rtt:=t*Sqrt(1/tt^2); fi;
#rt:=t^2;

if tt=0 then rtt:=infinity; fi;

if tt>0 then rtt:=t/tt;
  if IsSquareInt(NumeratorRat(rtt)) and IsSquareInt(DenominatorRat(rtt))
  then rtt:=Sqrt(rtt); else return []; fi;   #Is this return option valid?
fi;

rt:=t;


pqv1:=PairToQuadInt(z);
pqw1:=PairToQuadInt(zz);
##################################################
Potential_d:=function(v,w,c)
local r, czd,z, L;     
# tt =  t/(|cz+d|^2 + |c|^2t^2)
# t/tt = |cz+d|^2 + |c|^2t^2
#|cz+d|^2 = t/tt - |c|^2t^2

if w[2]=0 and not v[2]=w[2] then return []; fi;

if w[2]=0 then #v[2]=w[2] then
   #r:=1-HAPNorm(OQ,c)*v[2]^2;
   r:=1-HAPNorm(OQ,c)*rt;
else
   #r:=v[2]/w[2]-HAPNorm(OQ,c)*v[2]^2;
   #r:=Sqrt(v[2]^2/w[2]^2)-HAPNorm(OQ,c)*v[2]^2;
   #r:=v[2]*Sqrt(1/w[2]^2)-HAPNorm(OQ,c)*v[2]^2;
    r:=rtt-HAPNorm(OQ,c)*rt;
fi;


# want to return the set L  of all quadratic integers d such 
#that Norm(c*z+d) =r


#z:=PairToQuadInt(v[1]);
z:=pqv1;


L:=IntegerCircle(OQ,-c*z,r);   #Takes over 25% of the time


return List(L,xx->[c,xx]);
end;
##################################################


##################################################
Potential_b:=function(v,w,c,d)
local D,z,cczd,cct2,a, b, zz; #t,tt;    
if d=0 then return []; fi;

#z:=PairToQuadInt(v[1]);
z:=pqv1;
#t:=v[2];
#zz:=PairToQuadInt(w[1]);
zz:=pqw1;
#tt:=w[2];
#D:= HAPNorm(OQ,c*z+d) +HAPNorm(OQ,c)*t^2;
D:= HAPNorm(OQ,c*z+d) +HAPNorm(OQ,c)*rt;
cczd:=HAPComplexConjugate(c*z+d);
#cct2:=HAPComplexConjugate(c)*t^2;
cct2:=HAPComplexConjugate(c)*rt;

#So now we have
#(az+b)*cczd + a*cct2 = D*zz  ,
#ad-cb=1    a=(1+bc)/d .
#So
#( (1+bc)*z/d  +b ) * cczd  + (1+bc)/d * cct2  =D*zz
#  b*(c*z + d )* cczd / d + b(c*cct2)/d= D*zz - z*cczd/d - cct2/d
#  b*(c*z + d )* cczd  + b(c*cct2)= D*zz*d - z*cczd - cct2
#  b*( ((c*z + d )* cczd)  + c*cct2 ) = D*zz*d - z*cczd - cct2
#  b:= (D*zz*d - z*cczd - cct2) * ( (c*z + d )* cczd + c*cct2 )^-1

b:= (D*zz*d -z*cczd - cct2)*  ( (c*z + d )* cczd + c*cct2 )^-1  ;
if not IsHapQuadraticInteger(OQ,b) then return []; fi;
a:=(1+b*c)*(d^-1);
if not IsHapQuadraticInteger(OQ,a) then return []; fi;


return [a,b];
end;
##################################################

pot_c:=Potential_c([z,t],[zz,tt]);
pot_dc:=[];


for c in pot_c do
Append(pot_dc,Potential_d([z,t],[zz,tt],c) );
od;


##################################################
Standardize:=function(x);

return Maximum(x,-x);
end;
##################################################



pot_abcd:=[];
for x in pot_dc do
if not x=0 then
y:=Potential_b([z,t],[zz,tt],x[1],x[2]);
if Length(y)>0 then
Add(pot_abcd,Standardize([y[1],y[2],x[1],x[2]]));
fi;
fi;
od;




####################Case d=0#################
d:=0;c:=1;b:=-1; 
d:=zero; c:=one;b:=-one;
#(a*z - 1)*cz + a*t^2 = zz *D
# a*z*cz - cz +a*t^2 = zz *D  
# a*(z*cz + t^2) = zz * D +cz
# a = (D*zz +cz)/(Norm(z) + t^2);   #Let's assume t^2>0.

#x:=PairToQuadInt(z);
x:=pqv1;
cx:=HAPComplexConjugate(x);
#y:=PairToQuadInt(zz);
y:=pqw1;
#D:=HAPNorm(OQ,x)+t^2;
D:=HAPNorm(OQ,x)+t;
#a:= (D*y +cx)*(HAPNorm(OQ,x) + t^2)^-1;
a:= (D*y +cx)*(HAPNorm(OQ,x) + rt)^-1;
if IsHapQuadraticInteger(OQ,a) then
Add(pot_abcd,[a,b,c,d]);
fi;

if OQ!.bianchiInteger=-1 then
#(a*z + b)*cc*cz + a*cc*t^2 = zz *D
#a*z*cc*cz + b*cc*cz +a*cc*t^2 = D*zz
#a*(z*cc*cz +cc*t^2) = D*zz - b*cc*cz

d:=0;c:=HAPSqrt(-1);b:=c; cc:=-c;
d:=zero;
#a:=(D*y -b*cc*cx)*(HAPNorm(OQ,x)*cc + cc*t^2)^-1;
a:=(D*y -b*cc*cx)*(HAPNorm(OQ,x)*cc + cc*t)^-1;

if IsHapQuadraticInteger(OQ,a) then
Add(pot_abcd,[a,b,c,d]);
fi;

fi;;

if OQ!.bianchiInteger=-3 then
Print("Need to finish the case d=-3. units!!\n");
fi;;

#############################################


#pot_abcd:=Concatenation(pot_abcd,-pot_abcd);


pot_abcd:=List(pot_abcd,Standardize);

pot_abcd:=DuplicateFreeList(pot_abcd);  

#Print("pot_abcd",pot_abcd,"\n");

pot_abcd:=Filtered(pot_abcd,IsIntMat); 


Apply(pot_abcd, x->[[x[1],x[2]],[x[3],x[4]]]);

#Print("pot_c",pot_c,"\n");
#Print("pot_dc",pot_dc,"\n");

return pot_abcd;
end);
##############################################################
##############################################################

