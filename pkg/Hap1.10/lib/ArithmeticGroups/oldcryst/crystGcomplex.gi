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
#  r:=CanonicalRightCosetElement(H,Elts[g]^-1);
#  r:=id^-1*r;
#  u:=r*Elts[g];
r:=CanonicalRightCosetElement(H,Elts[g]);
r:=id^-1*r;
u:=Elts[g]*r^-1;
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
return StabGrp[m+1][k];
end;
# for i in [1..n+1] do
# for j in [1..Length(L[i])] do
# Append(Elts,Elements(StabGrp[i][j]));
# od;
# od;
#######################
PseudoRotSubGroup:=function(m,k)
local bdry,RotSbGrp,i,K,H,A;
if m=0 then 
#Print([m,k],StabGrp[m+1][k]);
return StabGrp[m+1][k];
else
bdry:=PseudoBoundary(m,k);
#Print([m,k],bdry,"\n");
bdry:=StructuralCopy(List(bdry,w->ConjugateGroup(RotSubGroupList[m][w[1]],Elts[w[2]])));
#Print([m,k],bdry,"\n");
RotSbGrp:=bdry[1];
#Print([m,k],RotSbGrp,"\n");
for i in [2..Length(bdry)] do
#H:=bdry[i];
#Print([m,k],H,"\n");
#K:=StructuralCopy(RotSbGrp);
#Print([m,k],K,"\n");
RotSbGrp:=Group(Intersection(Elements(RotSbGrp),Elements(bdry[i])));
#A:=Group(Intersection(Elements(H),Elements(K)));
#Print([m,k],A,"\n");
od;
return RotSbGrp;
fi;
end;
#####
RotSubGroupList:=[];
for i in [1..(n+1)] do
RotSubGroupList[i]:=[];
for j in [1..Length(L[i])] do
RotSubGroupList[i][j]:=PseudoRotSubGroup(i-1,j);
od;
od;
#####
RotSubGroup:=function(m,k)
return RotSubGroupList[m+1][k];
end;
#######################



return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=FinalBoundary,
	    PseudoBoundary:=PseudoBoundary,
	    RotSubGroupList:=RotSubGroupList,
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

# 
# InstallGlobalFunction(RecalculateIncidenceNumbers,
# function(R)
# local
#         NewBoundaryList,
#         NewBoundary,
#         S,N,pos,
#         cnt,n,i,b,e,V,v,L,LL,x,xx;
# R!.elts:=List(R!.elts,x->x^-1);
# 
# NewBoundaryList:=[];
# for n in [1..Length(R)] do
# NewBoundaryList[n]:=[];
# for i in [1..R!.dimension(n)] do
# b:=StructuralCopy(R!.boundary(n,i));
# #b:=List(b,x->[AbsInt(x[1]),x[2]]);
# Add(NewBoundaryList[n],b);
# od;
# od;
# 
# 
# ##############################
# NewBoundary:=function(n,i);
# if i>0 then
# return NewBoundaryList[n][i];
# else
# return NegateWord(NewBoundaryList[n][-i]);
# fi;
# end;
# ##############################
# 
# ########First Incidence Numbers###########
# for b in NewBoundaryList[1] do
# if b[1][1]*b[2][1]>0 then
# b[1][1]:=-1*b[1][1];
# fi;
# od;
# ##########################################
# 
# 
# R!.oldBoundary:=R!.boundary;
# R!.boundary:=NewBoundary;
# 
# for N in [2..Length(NewBoundaryList)] do
# ######## N-th Incidence Numbers##########
# cnt:=0;
# for e in NewBoundaryList[N] do
#   cnt:=cnt+1;
#   #########################
#   if Length(ResolutionBoundaryOfWord(R,N-1,e))>0 then
#   #########################
#     V:=List(e,x->[AbsInt(x[1]),x[2]]);
#     b:=[V[1]]; V:=V{[2..Length(V)]}; 
#     while Length(V)>0 do
#       L:=ResolutionBoundaryOfWord(R,N-1,b);
#       LL:=List(L,x->[AbsInt(x[1]),x[2]]);
#       for v in V do
#         x:=ResolutionBoundaryOfWord(R,N-1,[v]);
#         xx:=List(x,a->[AbsInt(a[1]),a[2]]);  
#         if Length(Intersection(xx,LL))>0 then
#           if Length(Intersection(x,L))>0 then Add(b,[-v[1],v[2]]); 
#              else Add(b,[v[1],v[2]]); 
#           fi;
#           pos:=Position(V,v);
#           V[pos]:=0;
#           V:=Filtered(V,v->not v=0); 
#           break; 
#         fi;
#       od;
#     od;
#     NewBoundaryList[N][cnt]:=b;
#   #########################
#   fi;
#   #########################
# od;
# ##########################################
# od;
# 
# end);
# ################################################