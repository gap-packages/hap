##############################################################
##############################################################
InstallGlobalFunction(HAP_BianchiTransformations,
function(OQ,V,W)
local PairToQuadInt,
      HAPquadtoquad,
      HAPComplexConjugate,
      Potential_c,
      IntegerCircle,
      Potential_d,
      Potential_b,
      Standardize,
      IsIntMat,
      one,zero,
      c,x,y,a,b,d,cx,D,G,z,zz,t,tt,cc,
      pot_c,pot_dc,pot_abcd;

#returns the matrices in SL(2,OQ) that send V to W. 
#If V=W then these are returned as a group. Here V,W are lists of length 3.

z:=[V[1],V[2]]; t:=V[3];
zz:=[W[1],W[2]]; tt:=W[3];

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
if not HAPquadtoquad(M[1]) in OQ then return false; fi;
if not HAPquadtoquad(M[2]) in OQ then return false; fi;
if not HAPquadtoquad(M[3]) in OQ then return false; fi;
if not HAPquadtoquad(M[4]) in OQ then return false; fi;
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
HAPquadtoquad:=function(x);
if IsRat(x) then return x; fi;
return x!.rational+x!.irrational*Sqrt(OQ!.bianchiInteger);
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
ttt:=(t*tt)^2;
ttt:=1+Int(ttt);
return QuadraticIntegersByNorm(OQ,1+1/(ttt));
return QuadraticIntegersByNorm(OQ,1/(t*tt));
end;
##################################################

##################################################
IntegerCircle:=function(OQ,c,r)
local C;
                   #STUPID IMPLEMENTATION!
#return all integers on the circle of radius Sqrt(r) and centre c.
C:=QuadraticIntegersByNorm(OQ,30);;
C:=Filtered(C,x->HAPNorm(OQ,-c+x)=r);
return C;
end;
##################################################

##################################################
Potential_d:=function(v,w,c)
local r, czd,z, L;

# tt =  t/(|cz+d|^2 + |c|^2t^2)
# t/tt = |cz+d|^2 + |c|^2t^2
#|cz+d|^2 = t/tt - |c|^2t^2

r:=v[2]/w[2]-HAPNorm(OQ,c)*v[2]^2;

# want to return the set L  of all quadratic integers d such 
#that Norm(c*z+d) =r

z:=PairToQuadInt(v[1]);
L:=IntegerCircle(OQ,-c*z,r);

return List(L,xx->[c,xx]);
end;
##################################################


##################################################
Potential_b:=function(v,w,c,d)
local D,z,t,cczd,cct2,a, b,zz,tt;
if d=0 then return[]; fi;

z:=PairToQuadInt(v[1]);
t:=v[2];
zz:=PairToQuadInt(w[1]);
tt:=w[2];
D:= HAPNorm(OQ,c*z+d) +HAPNorm(OQ,c)*t^2;
cczd:=HAPComplexConjugate(c*z+d);
cct2:=HAPComplexConjugate(c)*t^2;

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
a:=(1+b*c)*(d^-1);



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

x:=PairToQuadInt(z);
cx:=HAPComplexConjugate(x);
y:=PairToQuadInt(zz);
D:=HAPNorm(OQ,x)+t^2;
a:= (D*y +cx)*(HAPNorm(OQ,x) + t^2)^-1;
Add(pot_abcd,[a,b,c,d]);

if OQ!.bianchiInteger=-1 then
#(a*z + b)*cc*cz + a*cc*t^2 = zz *D
#a*z*cc*cz + b*cc*cz +a*cc*t^2 = D*zz
#a*(z*cc*cz +cc*t^2) = D*zz - b*cc*cz

d:=0;c:=HAPSqrt(-1);b:=c; cc:=-c;
d:=zero;
a:=(D*y -b*cc*cx)*(HAPNorm(OQ,x)*cc + cc*t^2)^-1;
Add(pot_abcd,[a,b,c,d]);

fi;;

if OQ!.bianchiInteger=-3 then
Print("Need to finish the case d=-3. units!!\n");
fi;;

#############################################

pot_abcd:=Concatenation(pot_abcd,-pot_abcd);
pot_abcd:=DuplicateFreeList(pot_abcd);
pot_abcd:=Filtered(pot_abcd,IsIntMat);

Apply(pot_abcd, x->[[x[1],x[2]],[x[3],x[4]]]);

return pot_abcd;
G:=Group(pot_abcd);


return G;
end);
##############################################################
##############################################################

