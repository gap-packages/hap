#######################################################
#######################################################
HAP_QuadToCyclotomic:=function(z);
return z!.rational+z!.irrational*Sqrt(z!.bianchiInteger);
end;
#######################################################
#######################################################

#######################################################
#######################################################
HAP_QuadToCyclotomicMat:=function(A)
local B,R,r,x;
B:=[];

for r in A do
R:=[];
for x in r do
Add(R,HAP_QuadToCyclotomic(x));
od;
Add(B,R);
od;

return B;
end;
#######################################################
#######################################################

#######################################################
#######################################################
HAP_QuadToCyclotomicGroup:=function(G)
local gens;

gens:=GeneratorsOfGroup(G);
gens:=List(gens,HAP_QuadToCyclotomicMat);

return Group(gens);
end;
#######################################################
#######################################################

#######################################################
#######################################################
QuadraticToCyclotomicCoefficients:=function(RR)
local R,BI,n,k,d,Q,OQ,I,G,L,pos;

for n in [1..Length(RR)] do
for k in [1..RR!.dimension(n)] do
RR!.boundary(n,k);
od;
od;

R:=Objectify(HapResolution,
                rec(
                dimension:=RR!.dimension,
                boundary:=RR!.boundary,
                homotopy:=RR!.homotopy,
                elts:=List(RR!.elts,HAP_QuadToCyclotomicMat),
                group:=fail,
                properties:= RR!.properties));

d:=One(RR!.elts[1][1][1]);;
d:=d!.bianchiInteger;
Q:=QuadraticNumberField(d);;
OQ:=RingOfIntegers(Q);;
I:=QuadraticIdeal(OQ,1);;
G:=HAP_CongruenceSubgroupGamma0(I);;
G!.tree:=true;
R!.group:=G;

return R;

end;
#######################################################
#######################################################


#######################################################
#######################################################
HAP_BianchiAction:=function(A,PP)
local P,ii, a,b,c,d,z,t,zz,tt,BI,B1,B2,ans,D,cc;
#We need P[1] and P[2] to be rationals and P[3] to be a cyclotomic.

P:=[];B1:=false;B2:=false;
if IsRat(PP[1]) then P[1]:=PP[1]; else B1:=true; P[1]:=PP[1]!.rational; fi;
if IsRat(PP[2]) then P[2]:=PP[2]; else B2:=true; P[2]:=PP[2]!.irrational; fi;
P[3]:=Sqrt(PP[3]);
BI:=A[1][1]!.bianchiInteger;

ii:=Sqrt(BI);
z:=P[1]+P[2]*ii; t:=P[3];


a:=A[1][1]; a:=HAP_QuadToCyclotomic(a);
b:=A[1][2]; b:=HAP_QuadToCyclotomic(b);
c:=A[2][1]; c:=HAP_QuadToCyclotomic(c); 
d:=A[2][2]; d:=HAP_QuadToCyclotomic(d);

cc:=ComplexConjugate(c*z+d);

#zz:=(a*z+b)*ComplexConjugate(c*z+d) +a*ComplexConjugate(c)*t^2;
zz:=(a*z+b)*cc +a*ComplexConjugate(c)*t^2;

D:=(c*z+d)*cc +(c*t)*ComplexConjugate(c*t);
if D=0 then return infinity; fi;

#zz:=zz/((c*z+d)*ComplexConjugate(c*z+d) +(c*t)*ComplexConjugate(c*t));
zz:=zz/D;


#tt:=t/((c*z+d)*ComplexConjugate(c*z+d) +(c*t)*ComplexConjugate(c*t));
tt:=t/D;



ans:=[];
if B1 then ans[1]:=QuadraticNumber(RealPart(zz),0,-BI);
else ans[1]:=RealPart(zz); fi;
if B2 then ans[2]:=QuadraticNumber(0,ImaginaryPart(zz)/Sqrt(-BI),-BI);
else ans[2]:=ImaginaryPart(zz)/Sqrt(-BI); fi;
ans[3]:=tt;

return ans;;

end;
#######################################################
#######################################################


#######################################################
#######################################################
HAP_BianchiRepresentativesAction:=function(Y,n,j,k)
local BianchiTrans, rep,Rec, Rec2, Rec3, Rec4, pos, OQ, T, TT, ans, e, f, A, ee, ff, fff, i,ii,iii,u,v;

OQ:=Y!.ring;

######################################
rep:=function(i);
return PositionProperty(Y!.ORBS[1],x->i in x);
end;
######################################

Rec:=Y!.BianchiRepresentativesActionRecords[1];
Rec2:=Y!.BianchiRepresentativesActionRecords[2];
Rec3:=Y!.BianchiRepresentativesActionRecords[3];
Rec4:=Y!.BianchiRepresentativesActionRecords[4];
pos:=Position(Rec3[n+1],[j,k]);
if not pos=fail then return Rec4[n+1][pos]; fi;
Add(Rec3[n+1],[j,k]);
pos:=Length(Rec3[n+1]);

######################################
BianchiTrans:=function(OQ,i,j)
local pos,p,q, A, T; 
pos:=Position(Rec,[i,j]);
if pos=fail then

####################
Add(Rec,[i,j]);
pos:=Length(Rec);

#if Length(Flat(Y!.ORBS[1]))<Y!.nrCells(0)  then 
if Length(Flat(Y!.ORBS[1]))<Length(Y!.TMP[1])  then
Rec2[pos]:=HAP_BianchiTransformations(OQ,Y!.points[i],Y!.points[j]);
else 

#p:=PositionProperty(Y!.ORBS[1],x->i in x);
#q:=PositionProperty(Y!.ORBS[1],x -> j in x);
p:=rep(i);
q:=rep(j);
if not p=q then
Rec2[pos]:=[];
else
#Rec2[pos]:=HAP_BianchiTransformations(OQ,Y!.points[i],Y!.points[j]);

if true then
##########################################
p:=Y!.ORBS[1][p][1];

if p=i then
A:=HAP_BianchiTransformations(OQ,Y!.points[p],Y!.points[i]);
else
A:=Rec2[Position(Rec,[p,i])];
fi;


if p=j then
T:=HAP_BianchiTransformations(OQ,Y!.points[p],Y!.points[j]);
else
T:=Rec2[Position(Rec,[p,j])];
fi;

Rec2[pos]:=T*A[1]^-1;;
##########################################
fi;


fi;
fi;
####################

Rec2[pos]:=Concatenation(Rec2[pos],-Rec2[pos]);             #??????
Rec2[pos]:=DuplicateFreeList(Rec2[pos]);                          #??????
fi;

return Rec2[pos];
end;
######################################

######################################
if n=0 then
if Y!.points[j][3]=0 and not Y!.points[k][3]=0 then 
Rec4[n+1][pos]:=[];
return []; fi;
if Y!.points[k][3]=0 and not Y!.points[j][3]=0 then 
Rec4[n+1][pos]:=[];
return []; fi;
T:=BianchiTrans(OQ,j,k); 
Rec4[n+1][pos]:=T;
    return T;
fi;
######################################

######################################
if n=1 then
ans:=[];
e:=Y!.boundaries[2][j];
e:=e{[2,3]};
f:=Y!.boundaries[2][k];
f:=f{[2,3]};

if not SortedList(List(e,x->rep(x))) = SortedList(List(f,x->rep(x))) then
ans:=[]; 

else 

#if Order(Group(BianchiTrans(OQ,e[1],e[1])))*Order(Group(BianchiTrans(OQ,e[2],e[2])))=infinity then Print(-BianchiTrans(OQ,e[2],e[2])[50] in (BianchiTrans(OQ,e[2],e[2])) ,"\n\n\n  "); fi;

fff:=List(f,i->Y!.points[i]);
T:=BianchiTrans(OQ,e[1],f[1]);
TT:=BianchiTrans(OQ,e[2],f[2]);
ans:=Intersection(T, TT );
#Print([Length(T),Length(TT),Length(Intersection(T,TT))], "  ");

#for A in T do
#if HAP_BianchiAction(A,Y!.points[e[2]]) in fff
#then Add(ans, A); 
#fi;
#od;

T:=BianchiTrans(OQ,e[1],f[2]);
TT:=BianchiTrans(OQ,e[2],f[1]);
Append(ans,Intersection(T, TT));
#Print([Length(T),Length(TT),Length(Intersection(T,TT)),Length(ans)], "  ");

#for A in T do
#if HAP_BianchiAction(A,Y!.points[e[2]]) in fff
#then Add(ans, A);
#fi;
#od;

fi;

Rec4[n+1][pos]:=ans;
return ans;

fi;
######################################

######################################
if n=2 then
ans:=[];
e:=Y!.boundaries[3][j];
e:=e{[2..Length(e)]};
f:=Y!.boundaries[3][k];
f:=f{[2..Length(f)]};
if (not Length(e)=Length(f)) then 
Rec4[n+1][pos]:=[];

return []; fi;

ee:=[];
for i in e do
Append(ee,Y!.boundaries[2][i]{[2,3]});
od;
ee:=SSortedList(ee);

ff:=[];
for i in f do
Append(ff,Y!.boundaries[2][i]{[2,3]});
od;
ff:=SSortedList(ff);

if not SortedList(List(ee,x->rep(x))) = SortedList(List(ff,x->rep(x)))
 then

Rec4[n+1][pos]:=[];

return []; fi;


fff:=List(ff,i->Y!.points[i]);

   for i in ff do
   for ii in ff do
   for iii in ff do
   if Length(SSortedList([i,ii,iii]))=3 then
      T:=BianchiTrans(OQ,ee[1],i);
      if Length(T)>0 then
      T:=Intersection(T, BianchiTrans(OQ,ee[2],ii));
      if Length(T)>0 then
      T:=Intersection(T,BianchiTrans(OQ,ee[3],iii));
      Append(ans, T);
      fi;fi;
      #for A in T do
      #   u:=HAP_BianchiAction(A,Y!.points[ee[2]]); 
      #   v:=HAP_BianchiAction(A,Y!.points[ee[3]]);

      #   if u in fff and v in fff then Add(ans, A); fi;
      #od;
   fi;
   od;
   od;
   od;
ans:=DuplicateFreeList(ans);

Rec4[n+1][pos]:=ans;
return ans;
fi;
######################################

return fail;
end;
#######################################################
#######################################################

#######################################################
#######################################################
BianchiGcomplex:=function(d)
local R,P,OQ,Y,Dimension,Stabilizer,Action,STABS,stb,rot,rotREC,Boundary,rep,ELTS,
G,K,S,T,i,V,n,B,BB,k,BoundaryRec,TMP, EquivSpheres;

if not IsInt(d) then
Print("input must be a negative square free integer.\n");
return fail; fi;
if not d<0 then
Print("input must be a negative square free integer.\n");
return fail; fi;
if not Length(Factors(d))=Length(SSortedList(Factors(d))) then
Print("input must be a negative square free integer.\n");
return fail; fi;

P:=BianchiPolyhedron(d);
Y:=HAP_BianchiRegularCWComplex(P!.ring,P!.unimodularPairs);
OQ:=Y!.ring;
STABS:=[[],[],[]];
ELTS:=[];

Y!.ORBS:=[[],[],[]];
Y!.BianchiRepresentativesActionRecords:=[[],[],[[],[],[]],[[],[],[]]];


###############################################
rep:=function(n,k)
local pos;
return PositionProperty(Y!.ORBS[n+1],x->k in x);
end;
###############################################

###############################################
stb:=function(n,i)
local G,gens,g,c,rnk,GENS;
if n=0 and Y!.points[i][3]=0 then 
G:=HAP_BianchiRepresentativesAction(Y,n,i,i);
   ###################################
   gens:=[[0,0]]; GENS:=[];
   rnk:=0;
   for g in G do
      if g^2<>One(g) then
      c:=g[2][1];
      if not IsZero(c) then
         c:=[c!.rational,c!.irrational];
         if Rank(Concatenation(gens,[c]))>rnk then
         Add(gens,c); Add(GENS,g); rnk:=rnk+1;
         fi;
      fi;
      fi;
   if rnk=2 then break; fi;
   od;
   Add(GENS,-One(GENS[1]));
   ###################################
STABS[1][i]:= Group(GENS); STABS[1][i]!.Order:=infinity; STABS[1][i]!.Size:=infinity;
return infinity; fi;
G:=HAP_BianchiRepresentativesAction(Y,n,i,i);

STABS[n+1][i]:= Group(G); 
return IdGroup(STABS[n+1][i]);
end;
###############################################

###############Orbits pre-computation##########
TMP:=[[],[],[]];

if (d mod 4) <>1 then
#########################
EquivSpheres:=function(u,v)
local w;
w:=u-v;
if IsInt(w[1]) and IsInt(w[2]) then return true; fi;
return false;
end;
#########################
else
#########################
EquivSpheres:=function(u,v)
local w;
w:=u-v;
w:=[w[1],w[2]];
if IsInt(w[1]-w[2]) and IsInt(2*w[2]) then return true; fi;
return false;
end;
#########################
fi;

V:=[];
S:=[1..Y!.nrCells(2)];
while Length(S)>0 do
T:=[S[1]];
     for i in S{[2..Length(S)]} do
     if EquivSpheres(Y!.sphereCentres[S[1]],Y!.sphereCentres[i]) then Add(T,i); fi;
     od;
Add(V,T);
S:=Difference(S,T);
od;
V:=List(V,x->Minimum(x));
V:=SSortedList(V);
TMP[3]:=V;

V:=List(V,f->Y!.boundaries[3][f]);
V:=List(V,v->v{[2..v[1]+1]});
V:=Flat(V);
V:=SSortedList(V);
TMP[2]:=V;

V:=List(V,f->Y!.boundaries[2][f]);
V:=List(V,v->v{[2,3]});
V:=Flat(V);
V:=SSortedList(V);
TMP[1]:=V;

Y!.TMP:=TMP;


###############Orbits pre-computation##########


for n in [0,1,2] do
V:=[];
#S:=[1..Y!.nrCells(n)];
S:=TMP[n+1];
while Length(S)>0 do
T:=[S[1]];
     for i in S{[2..Length(S)]} do
        if n=0 then
        G:=HAP_BianchiRepresentativesAction(Y,n,S[1],i);
        if Length(G)>0
           then Add(T,i);
        fi;
        fi;

        if n=1 then
        if Length(HAP_BianchiRepresentativesAction(Y,n,S[1],i))>0
           then Add(T,i);
        fi;
        fi;

        if n=2 then
        B:=Y!.boundaries[3][S[1]];
        B:=B{[2..B[1]]};
        B:=List(B,x->Y!.boundaries[2][x]{[2,3]});
        B:=Set(Flat(B));
        
        BB:=Y!.boundaries[3][i];
        BB:=BB{[2..BB[1]]};
        BB:=List(BB,x->Y!.boundaries[2][x]{[2,3]});
        BB:=Set(Flat(BB));

        if Length(HAP_BianchiRepresentativesAction(Y,n,S[1],i))>0
           then Add(T,i);
        fi;
        fi;
     od;
Add(V,T);
S:=Difference(S,T);
od;
Y!.ORBS[n+1]:=V;
od;

#######################################
Dimension:=function(n);
if not n in [0,1,2] then return 0; fi;
return Length(Y!.ORBS[n+1]);
end;
#######################################

for n in [0,1,2] do
for i in Y!.ORBS[n+1] do
stb(n,i[1]);od;od;



#######################################
Stabilizer:=function(n,i);
return STABS[n+1][Y!.ORBS[n+1][i][1]];
end;
#######################################

rotREC:=[[],[]];
#######################################
rot:=function(n,i)
local S, R, g, prm,bnd,bbnd,x,y,H,j,P,PP;

if IsBound(rotREC[n][i]) then return rotREC[n][i]; fi;

S:= STABS[n+1][Y!.ORBS[n+1][i][1]];  


########################
if n=1 then 
bnd:=Y!.boundaries[n+1][Y!.ORBS[n+1][i][1]];
bnd:=bnd{[2..bnd[1]+1]};
R:=[]; 
stb(n-1,bnd[1]); 
H:=STABS[n][bnd[1]];
Order(H);  #Next step needs to have the order already computed!!
for x in Elements(S) do
if x in H then Add(R,x); fi;
od;
R:=Group(R);
rotREC[n][i]:=R;
fi;
########################

########################
if n=2 then
j:=Y!.ORBS[n+1][i][1];
P:=1*Y!.sphereCentres[j];
P[3]:=0;
R:=[];
for x in Elements(S) do
PP:= HAP_BianchiAction(x,P);
if PP=P then Add(R,x); fi;
od;
R:=Group(R);
rotREC[n][i]:=R;
fi;
########################

return rotREC[n][i];;
end;
#######################################


#######################################
Action:=function(n,k,g)
local id, r, u, ans, abk, H;

if n=0 then  return 1; fi;

abk:=AbsInt(k);
H:=Stabilizer(n,abk);

if Order(H)=infinity then return 1; fi; #Assuming infinite groups act trivially??

id:=CanonicalRightCosetElement(H,Identity(H));
r:=CanonicalRightCosetElement(H,ELTS[g]^-1);
r:=id^-1*r;
u:=r*ELTS[g];

if u in rot(n,abk) then  ans:= 1;
else ans:= -1; fi;

return ans;   

end;
#######################################

#######################################
Boundary:=function(n,kk)
local k,e, i,x, bnd, ans,p,g, gg,pos, z;

k:=AbsInt(kk);
e:=Y!.ORBS[n+1][k][1];

ans:=[];
bnd:=1*Y!.boundaries[n+1][e];
bnd:=bnd{[2..bnd[1]+1]};

for i in [1..Length(bnd)] do
x:=bnd[i];
p:=rep(n-1,x);
gg:=HAP_BianchiRepresentativesAction(Y,n-1,Y!.ORBS[n][p][1],x);
g := gg[1]; 

#g:=CanonicalRightCosetElement(Stabilizer(n-1,p),g^-1)^-1;

pos:=Position(ELTS,g);
if pos=fail then Add(ELTS,g); pos:=Length(ELTS); fi;

Add(ans,[p,pos]);
od;


if SignInt(kk)>0 then return ans;
else return NegateWord(ans); fi;
end;
#######################################

BoundaryRec:=[];
for n in [1..2] do
BoundaryRec[n]:=[];
for k in [1..Dimension(n)] do
BoundaryRec[n][k]:=Boundary(n,k);
od;
od;

#######################################
Boundary:=function(n,k);
if k>0 then return BoundaryRec[n][k];
else return NegateWord(BoundaryRec[n][-k]); fi;
end;
#######################################



R:= Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=ELTS,
            group:=Group(ELTS),
            stabilizer:=Stabilizer,
            action:=Action,
            cwSpace:=Y,
            properties:=
            [["length",100000],
             ["characteristic",0],
             ["type","resolution"],
             ["reduced",false]]  ));

RecalculateIncidenceNumbers_NonFreeRes(R);
return R;
end;
#######################################################
#######################################################



