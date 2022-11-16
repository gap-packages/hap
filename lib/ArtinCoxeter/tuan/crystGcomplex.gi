InstallGlobalFunction(CrystGcomplex,
function(gens,basis,check)
local i,x,k,combin,n,j,r,m,vect,c,
      B,G,T,S,Bt,Action,Sign,FinalBoundary,BoundaryList,
      L,kcells,cells,w,StabGrp,ActionRecord,lnth,PseudoRotSubGroup,RotSubGroupList,
      Dimension,SearchOrbit,pos,StabilizerOfPoint,PseudoBoundary,RotSubGroup,
      Elts,Boundary,Stabilizer,DVF,DVFRec,Homotopy,rmult,FinalHomotopy;
B:=basis[1];
c:=basis[2];
vect:=c-Sum(B)/2;

vect:=0*vect;

G:=AffineCrystGroup(gens);
T:=TranslationSubGroup(G);
Bt:=T!.TranslationBasis;
S:=RightTransversal(G,T);
n:=DimensionOfMatrixGroup(G)-1;
Elts:=[One(G)];
Append(Elts,gens);
lnth:=1000;


if check=1 then    # B is the G-full basis
#######################
L:=[];
for k in [0..n] do
L[k+1]:=[];

###   list all centers of k-cells

kcells:=[];
combin:=Combinations([1..n],k);
for x in combin do
  w:=[];
  for i in [1..n] do
    if i in x then
      Add(w,[1/2]);
    else Add(w,[0,1]);
    fi;
  od;
  cells:=Cartesian(w);
  Append(kcells,cells*B+vect);
od;
###  search for k-orbits 
Add(L[k+1],kcells[1]);
for i in [2..Length(kcells)] do
r:=0;
for j in [1..Length(L[k+1])] do
  if IsList(IsCrystSameOrbit(G,Bt,S,kcells[i],L[k+1][j])) then
    break;
  fi;
  r:=r+1;
od;

if r=Length(L[k+1]) then Add(L[k+1],kcells[i]);fi;

od;


od;
##########



##############
elif check=0 then   #slice the fundamental cell into 2^n parts to get a proper action of G on R^n
B:=List(B,x->x/2); 
L:=[];
for k in [0..n] do
L[k+1]:=[];

###   list all centers of k-cells

kcells:=[];
combin:=Combinations([1..n],k);
for x in combin do
  w:=[];
  for i in [1..n] do
    if i in x then
      Add(w,[1/2,3/2]);
    else Add(w,[0,1,2]);
    fi;
  od;
  cells:=Cartesian(w);
  Append(kcells,cells*B+vect);
od;

###  search for k-orbits 
Add(L[k+1],kcells[1]);
for i in [2..Length(kcells)] do
r:=0;
for j in [1..Length(L[k+1])] do
  if IsList(IsCrystSameOrbit(G,Bt,S,kcells[i],L[k+1][j])) then
    break;
  fi;
  r:=r+1;
od;

if r=Length(L[k+1]) then Add(L[k+1],kcells[i]);fi;

od;


od;
#########



else 
Print("check is either 1 for B is G-full basis and 0 for proper action", "\n");
return fail;
fi;
#######################
Dimension:=function(k)
if k>n then return 0;fi;
return Length(L[k+1]);
end;
#######################
pos:=function(g)
local p;
p:=Position(Elts,g);
if p=fail then 
Add(Elts,g);
return Length(Elts);
else return p;
fi;
end;
#######################
SearchOrbit:=function(g,k)
local i,p,h;
for i in [1..Length(L[k+1])] do
  p:=IsCrystSameOrbit(G,Bt,S,L[k+1][i],g);
  if IsList(p) then 
    h:=pos(p);
    return [i,h];fi;
od;
end;



ActionRecord:=[];
for m in [1..lnth+1] do
ActionRecord[m]:=[];
for k in [1..Dimension(m-1)] do
ActionRecord[m][k]:=[];
od;
od;


#############
rmult:=function(L,k,g)
local x,w,t,h,y,vv;
vv:=[];
for x in [1..Length(L)] do
w:=Elts[L[x][2]]*Elts[g];
L[x][1]:=Sign(k,L[x][1],pos(w))*L[x][1];
w:=CanonicalRightCosetElement(StabGrp[k+1][AbsInt(L[x][1])],w);
t:=pos(w);
Add(vv,[Sign(k,L[x][1],t)*L[x][1],t]);
od;
return vv;
end;
#############
#######################
Action:=function(m,k,g)
local id,r,u,H,abk,ans,x,h,l,i;

abk:=AbsInt(k);

if not IsBound(ActionRecord[m+1][abk][g]) then 
H:=StabGrp[m+1][abk];

if Order(H)=infinity then ActionRecord[m+1][abk][g]:=1;
#So we are assuming that any infinite stabilizer group acts trivially!!
else
######
id:=CanonicalRightCosetElement(H,Identity(H));
 r:=CanonicalRightCosetElement(H,Elts[g]^-1);
 r:=id^-1*r;
 u:=r*Elts[g];
# r:=CanonicalRightCosetElement(H,Elts[g]);
 #r:=id^-1*r;
# u:=r*Elts[g]^-1*id;
########

if u in RotSubGroupList[m+1][abk] then  ans:= 1;
else ans:= -1; fi;

ActionRecord[m+1][abk][g]:=ans;
fi;
######
fi;

return ActionRecord[m+1][abk][g];

end;
#######################
PseudoBoundary:=function(k,s)
local f,x,bdry,i,Fnt,Bck,j,ss;
ss:=AbsInt(s);
f:=L[k+1][ss];
if k=0 then return [];fi;
#x:=f*B^-1;
x:=(f-vect)*B^-1;
bdry:=[];
j:=0;
for i in [1..n] do
Fnt:=StructuralCopy(x);
Bck:=StructuralCopy(x);
if not IsInt(x[i]) then
j:=j+1;
Fnt[i]:=Fnt[i]-1/2;
Bck[i]:=Bck[i]+1/2;
#Fnt:=Fnt*B;
#Bck:=Bck*B;

Fnt:=Fnt*B+vect;
Bck:=Bck*B+vect;
Append(bdry,[SearchOrbit(Fnt,k-1),SearchOrbit(Bck,k-1)]);
#Append(bdry,[SearchOrbit(Fnt,k-1),SearchOrbit(Bck,k-1)]);

fi;
od;
return bdry;
end;
#######################
Sign:=function(m,k,g)
local x,h,p,r,c,i,y,f,s,kk,e,B1,B2,w;
kk:=AbsInt(k);
if m=0 then return 1;fi;
h:=Elts[g];
p:=CrystFinitePartOfMatrix(h);
e:=L[m+1][kk];
#x:=e*B^-1;
x:=e*B^-1;
r:=[];
for i in [1..Length(x)] do
if not IsInt(x[i]) then
Add(r,i);
fi;
od;
B1:=B{r};
B1:=B1*p;

e:=Flat(e);
Add(e,1);
f:=e*h;
Remove(f);
y:=f*B^-1;
c:=[];
for i in [1..Length(y)] do
if not IsInt(y[i]) then
Add(c,i);
fi;
od;

B2:=B{c};
s:=[];
for i in [1..Length(B2)] do
Add(s,SolutionMat(B1,B2[i]));
od;
#Print(s);

return SignInt(Determinant(s));
end;
#######################
Boundary:=function(k,s)
local psbdry,j,w,bdry;
psbdry:=PseudoBoundary(k,s);
bdry:=[];
for j in [1..Length(psbdry)] do
w:=psbdry[j];
if (j mod 4 = 3) or (j mod 4 = 2) then
#if IsEvenInt(j) then
Add(bdry,Negate([Sign(k-1,w[1],w[2])*w[1],w[2]]));
else Add(bdry,[Sign(k-1,w[1],w[2])*w[1],w[2]]);
fi;

od;


if s<0 then return NegateWord(bdry);
else
return bdry;
fi;

end;
########################
BoundaryList:=[];
for i in [1..n] do
BoundaryList[i]:=[];
for j in [1..Dimension(i)] do
BoundaryList[i][j]:=Boundary(i,j);
od;
od;
#######################
FinalBoundary:=function(n,k)
if k>0 then return BoundaryList[n][k];
else return NegateWord(BoundaryList[n][AbsInt(k)]);
fi;
end;

##################################################
StabilizerOfPoint:=function(g)
local H,stbgens,i,h,p;
g:=Flat(g);
Add(g,1);
stbgens:=[];
for i in [1..Length(S)] do
h:=g*S[i]-g;
Remove(h);
p:=h*Bt^-1;
if IsIntList(p) then Add(stbgens,S[i]*VectorToCrystMatrix(h)^-1);fi;
od;
H:=Group(stbgens);
return H;
end;
###
StabGrp:=[];
for i in [1..(n+1)] do
StabGrp[i]:=[];
for j in [1..Length(L[i])] do
StabGrp[i][j]:=StabilizerOfPoint(L[i][j]);
od;
od;
###

Stabilizer:=function(m,k)
local kk;
kk:=AbsInt(k);
return StabGrp[m+1][k];
end;
##########################
PseudoRotSubGroup:=function(m,k)
local x,kk,l,h,i,w,r,y,H,id,eltsH,g,RotSbGrp;
kk:=AbsInt(k);
RotSbGrp:=[];
H:=StabGrp[m+1][k];
eltsH:=Elements(H);

for g in eltsH do
if Sign(m,k,pos(g))=1 then Add(RotSbGrp,g);fi;
od;
RotSubGroupList[m+1][kk]:=Group(RotSbGrp);
return Group(RotSbGrp);
end;
#######################
RotSubGroupList:=[];
for i in [1..(n+1)] do
RotSubGroupList[i]:=[];
for j in [1..Length(L[i])] do
RotSubGroupList[i][j]:=PseudoRotSubGroup(i-1,j);
od;
od;
#######################
RotSubGroup:=function(m,k)
local kk;
kk:=AbsInt(k);
return RotSubGroupList[m+1][kk];
end;
######################

DVFRec:=[];
for k in [1..n+1] do
DVFRec[k]:=[];
for i in [1..Length(L[k])] do
DVFRec[k][i]:=[];
od;od;
######################
if check=1 then


DVF:=function(k,w)    #input an n-cell acts like the starting point of an arrow
		    # the function returns n+1-cell acts like the end point of the above arrow
		    #those cells presented by its center
local f,x,g,i,y,ww,s,b,j;
ww:=[AbsInt(w[1]),w[2]];
if not IsBound(DVFRec[k+1][ww[1]][ww[2]]) then

x:=StructuralCopy(L[k+1][ww[1]]);
Add(x,1);
x:=x*Elts[ww[2]];
Remove(x);

f:=(x-vect)*B^-1;
#Print("test  ",f);
for i in [1..n] do
  if not f[i]=0 then
    if not IsInt(f[i]) then 
    DVFRec[k+1][ww[1]][ww[2]]:=[];
    return DVFRec[k+1][ww[1]][ww[2]];
    else 
	s:=SignInt(f[i]);
	f[i]:=f[i]-s*1/2;
	x:=f*B;
	y:=SearchOrbit(x,k+1);
	y[2]:=pos(CanonicalRightCosetElement(StabGrp[k+2][y[1]],Elts[y[2]]));

	DVFRec[k+1][ww[1]][ww[2]]:=y;
    return DVFRec[k+1][ww[1]][ww[2]];
    fi;
  fi;
od;
DVFRec[k+1][ww[1]][ww[2]]:=[];
return DVFRec[k+1][ww[1]][ww[2]];
else
return DVFRec[k+1][ww[1]][ww[2]];
fi;
end;
########
Homotopy:=function(k,w)
local h,d,x,y,i,ww,b,p1,p2,s1,s2,v,s,p,t,a,u;

if w=[] then return [];fi;
a:=Sign(AbsInt(k),w[1],w[2]);
d:=[];
w[2]:=pos(CanonicalRightCosetElement(StabGrp[k+1][AbsInt(w[1])],Elts[w[2]]));
w[1]:=a*Sign(k,w[1],w[2])*w[1];
ww:=[AbsInt(w[1]),w[2]];
h:=StructuralCopy(DVF(k,ww));

if h=[] then return [];fi;

x:=PseudoBoundary(k+1,h[1]);

u:=List(x,v->[v[1],Elts[v[2]]*Elts[h[2]]]);
u:=List(u,v->[v[1],pos(CanonicalRightCosetElement(StabGrp[k+1][AbsInt(v[1])],v[2]))]);
p:=Position(u,ww);

s:=1;;

b:=StructuralCopy(FinalBoundary(k+1,h[1]));

b:=rmult(b,k,h[2]);

c:=StructuralCopy(b);
t:=SignInt(b[p][1]);
Remove(c,p);

Add(d,h);
for i in [1..Length(c)] do
	Append(d,NegateWord(Homotopy(k,c[i])));
od;

if w[1]*t<0 then return NegateWord(d);
else
return d;
fi;


end;

##########
else 
DVF:=fail;
FinalHomotopy:=fail;

fi;
###########################################

###########################################
return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=FinalBoundary,
	    PseudoBoundary:=PseudoBoundary,
	    dvf:=DVF,
	    CellList:=L,
	    Sign:=Sign,
            homotopy:=Homotopy,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
	    RotSubGroup:=RotSubGroup,
            properties:=
            [["length",100],
             ["characteristic",0],
             ["type","resolution"]]  ));

end);


###############################################################




InstallGlobalFunction(ResolutionCubicalCrystGroup,
function(G,n)
local gens,B,C,R,Gram, pos, Homotopy,Cnew;
Gram:=GramianOfAverageScalarProductFromFiniteMatrixGroup(PointGroup(G));
if Gram=IdentityMat(DimensionOfMatrixGroup(PointGroup(G))) then
gens:=GeneratorsOfGroup(G);
G:=AffineCrystGroup(gens);
B:=CrystGFullBasis(G);
if IsList(B) then


C:=CrystGcomplex(gens,B,1);
Cnew:=CrystGcomplex(gens,B,1);
Apply(Cnew!.elts,x->x^-1);

##########
pos:=function(L,g)
local p;
p:=Position(L,g);
if p=fail then 
Add(L,g);
return Length(L);
else return p;
fi;
end;
###########
Homotopy:=function(n,w)
local p,h;
p:=pos(C!.elts,Cnew!.elts[w[2]]^-1);
h:=StructuralCopy(C!.homotopy(n,[w[1],p]));
Apply(h,x->[x[1],pos(Cnew!.elts,C!.elts[x[2]]^-1)]);
return h;
end;
###########
Cnew!.homotopy:=Homotopy;

R:=FreeZGResolution(Cnew,n);
return R;
else return fail;
fi;
else 
Print("Gramian matrix is not identity \n");
return fail;
fi;
end);

#############################################################
InstallGlobalFunction(BredonChainComplex,
function(C)
local StabIrrTable,i,j,N,
      Dimension,PairToTriple,BoundaryMatrix,Boundary,
      TripleToPair,StabGrp,BoundaryRec;

####
############
StabGrp:=[];
i:=0;
while C!.dimension(i)>0 do
    StabGrp[i+1]:=[];
    for j in [1..C!.dimension(i)] do
	Add(StabGrp[i+1],C!.stabilizer(i,j));
    od;
    i:=i+1;
od;
##############
StabIrrTable:=[];
i:=0;
while C!.dimension(i)>0 do
    StabIrrTable[i+1]:=[];
    for j in [1..C!.dimension(i)] do
	Add(StabIrrTable[i+1],OrdinaryCharacterTable(StabGrp[i+1][j]));
    od;
    i:=i+1;
od;


N:=i-1;

############
Dimension:=function(k)
local d,i;
d:=0;
for i in [1..C!.dimension(k)] do
    d:=d+Size(StabIrrTable[d+1][i]);
od;
return d;
end;
############
PairToTriple:=function(i,j)
local k,x;
k:=j;
x:=1;
while k>Size(StabIrrTable[i+1][1]) do
    k:=k-Size(StabIrrTable[i+1][x]);
    x:=x+1;
od;
return [i,x,k];
end;
############
TripleToPair:=function(i,j,k)
local d,x;
d:=0;
for x in [1..(j-1)] do
    d:=d+Size(StabIrrTable[i+1][x]);
od;
d:=d+k;
return [i,k];
end;
############
BoundaryMatrix:=function(n,k)
local bdry,x,Coeffs,Mat,W,A,B,i,xx;
bdry:=C!.boundary(n,k);
Mat:=[];
for i in [1..Length(bdry)] do
    x:=bdry[i][1];
    xx:=AbsInt(x);
    B:=StabGrp[n][xx];
    A:=OrdinaryCharacterTable(ConjugateGroup(B,C!.elts[bdry[i][2]]^-1));
    W:=Induced(StabIrrTable[n+1][k],A,Irr(StabIrrTable[n+1][k]));
    Coeffs:=MatScalarProducts(A,Irr(A),W);     
    Add(Mat,[SignInt(x),xx,Coeffs]);
od;
return Mat;
end;
#############
BoundaryRec:=[];
for i in [1..N] do
    BoundaryRec[i]:=[];
    for j in [1..C!.dimension(i)] do
	Add(BoundaryRec[i],BoundaryMatrix(i,j));
    od;
od;
################################
Boundary:=function(n,k)
local w,x,y;

w:=PairToTriple(n,k);
x:=BoundaryRec[n][k];


end;
#####################################
return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundarymatrix:=BoundaryMatrix,
            boundary:=Boundary,
	    #PseudoBoundary:=PseudoBoundary,

            homotopy:=fail,
            #elts:=Elts,
            group:=Integers,
            #stabilizer:=Stabilizer,
            #action:=Action,
	    #RotSubGroup:=RotSubGroup,
            properties:=
            [["length",1000],
             ["characteristic",0],
             ["type","resolution"]]  ));

end);
