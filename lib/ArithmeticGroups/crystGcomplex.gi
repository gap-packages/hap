InstallGlobalFunction(CrystGcomplex,
function(gens,basis)
local i,x,k,combin,n,j,r,m,
      B,G,T,S,Bt,Action,Sign,FinalBoundary,BoundaryList,
      L,kcells,cells,w,StabGrp,ActionRecord,lnth,PseudoRotSubGroup,RotSubGroupList,
      Dimension,SearchOrbit,pos,StabilizerOfPoint,PseudoBoundary,RotSubGroup,
      Elts,Boundary,Stabilizer;
B:=basis;
G:=AffineCrystGroup(gens);
T:=TranslationSubGroup(G);
Bt:=T!.TranslationBasis;
S:=RightTransversal(G,T);
n:=DimensionOfMatrixGroup(G)-1;
Elts:=[One(G)];
Append(Elts,gens);
lnth:=1000;
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
  Append(kcells,cells*B);
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


#######################
# Action:=function(n,k,g)
# local x,kk,l,h,i,w,r,y,H,id;
# kk:=AbsInt(k);
# h:=Elts[g];
# x:=(L[n+1][kk])*B^-1;
# l:=[];
# for i in [1..Length(x)] do
# if not IsInt(x[i]) then
# Add(l,i);
# fi;
# od;
# w:=h{l}{l};
# if IsMatrix(w) and Determinant(w)=-1 then return -1;
# else return 1;
# fi;
# end;
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
x:=f*B^-1;
bdry:=[];
j:=0;
for i in [1..n] do
Fnt:=StructuralCopy(x);
Bck:=StructuralCopy(x);
if not IsInt(x[i]) then
j:=j+1;
Fnt[i]:=Fnt[i]-1/2;
Bck[i]:=Bck[i]+1/2;
Fnt:=Fnt*B;
Bck:=Bck*B;
Append(bdry,[SearchOrbit(Fnt,k-1),SearchOrbit(Bck,k-1)]);

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

return SignRat(Determinant(s));
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
###########################################

return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=FinalBoundary,
	    PseudoBoundary:=PseudoBoundary,
#	    RotSubGroupList:=RotSubGroupList,
	    CellList:=L,
	    Sign:=Sign,
            homotopy:=fail,
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

