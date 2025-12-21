###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalGamma0SubgroupsIdeal,
function(G,H)
local  N, L, R, T, F, poscan, t, one, c;

N:=H!.level;
T:=RightTransversal(N);
R:=[];
for t in T do
Add(R, [[1 mod N,0 mod N],[t mod N,1 mod N]]);
od;

Add(R, [[0 mod N,1 mod N],[-1 mod N,0 mod N]]);

R:=SSortedList(R);

L:=Length(R);
F:=AssociatedRing(N) mod N;;
one:=One(F);

##########################################
poscan:=function(x)
local i,xx, a, b, c, d, A, B, C, D, pos;
xx:=x^-1;
#for i in [1..L]  do
#if R[i]*xx in H then return i; fi;
#od;
if not IsZero(x[2][2]*one) then
   c:=(x[2][1]*one)*(x[2][2]*one)^-1; c:=c![1];
   #return Position(R,[[1 mod N,0],[c,1 mod N]]);  
   return Position(T,c);
fi;
   #return Position(R,[[0 mod N,1 mod N],[-1 mod N,0 mod N]]); 
   return L;
end;
##########################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=R,
               poscan:=poscan ));

end);
###################################################################
###################################################################

############################################################
############################################################
InstallOtherMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsHapSL2ZSubgroup,IsHapSL2ZSubgroup],
0,  #There must be a better way to ensure this method is not used!
function(H,HH)
Print("Using   HAP_TransversalCongruenceSubgroups\n");
return HAP_TransversalCongruenceSubgroups(H,HH);

end);
############################################################
############################################################

if false then
###################################################################
###################################################################
InstallGlobalFunction(HAP_TransversalCongruenceSubgroups,
function(G,H)
local gensG, gensH, N, GN, iso, HN, R, R2, x, sln, one, epi, epi2, poscan;

#AT PRESENT THIS APPROACH IS SLOW AND SO RANKED VERY LOW AS A METHOD
Print("HAP should not be using this method.\n");
gensH:=GeneratorsOfGroup(H);
for x in gensH do
if not x in G then
Print("The second argument is not a subgroup of the first.\n");
return fail;fi;
od;
gensG:=GeneratorsOfGroup(G);
N:=H!.level;

sln:=SL(2,Integers mod N);
one:=One(sln);
GN:=Group(gensG*one);
iso:=IsomorphismPermGroup(GN);
epi:=GroupHomomorphismByImagesNC(G,Image(iso),gensG,List(gensG*one,x->Image(iso,x)));
epi2:=GroupHomomorphismByFunction(G,GN,x->Image(iso,x*one));
HN:=Group(List(gensH*one,x->Image(iso,x)));

R:=RightTransversal(Image(iso),HN);
R2:=List(R,x->PreImagesRepresentative(epi,x));

##########################################
poscan:=function(x);
return PositionCanonical(R,ImagesRepresentative(epi2,x));
end;
##########################################
return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=R2,
               poscan:=poscan ));
end);
###################################################################
###################################################################
fi;


############################################################
############################################################
InstallGlobalFunction(HAP_RightTransversalSL2ZSubgroups,
function(H,HH)
local F, rels, S, T, G, FhomG, Q, gensQ, epi, ElementToWord,
gensH, gensHH, QH, QHH, R, R1, R2, s, t, poscan, iso;

G:=SL(2,Integers);
F:=FreeGroup(2);s:=F.1;t:=F.2;
rels:=[s^4,(s*t)^3*s^-2];
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
Q:=F/rels;
FhomG:=GroupHomomorphismByImages(F,G,[s,t],[S,T]);
gensQ:=GeneratorsOfGroup(Q);
epi:=GroupHomomorphismByImages(F,Q,GeneratorsOfGroup(F),GeneratorsOfGroup(Q));

###########################################################
###########################################################
ElementToWord:=function(g)
local A, d,i, b,L, bool;

#bool is true if S^n was the last added element
bool:=false;

L:=[One(F)];
A:=1*g;

while not A[2][1]=0 do
if AbsInt(A[2][1])>AbsInt(A[1][1]) then
   A:=S*A; Add(L,s);  bool:=true;
fi;
if  not A[2][1]=0 then
  if A[1][1]/A[2][1]>=0 then
    A:=T^-1*A;
    Add(L,t^-1);  bool:=false;
  else
    A:=T*A;  bool:=false;
    Add(L,t);
  fi;
fi;
od;




d:=-A[1][1]*A[1][2];


if not bool then
   A:=T^SignInt(d)*A;

   if SignInt(d)=1 then
      Add(L,t);
   fi;
   if SignInt(d)=-1 then
      Add(L,t^-1);
   fi;
else
   A:=T^SignInt(d)*A;
   bool:=false;
   if SignInt(d)=1 then
       Add(L,t);
   fi;
   if SignInt(d)=-1 then
       Add(L,t^-1);
   fi;
fi;

for i in [2..AbsInt(d)] do
A:=T^SignInt(d)*A;

if SignInt(d)=1 then
Add(L,t);
fi;
if SignInt(d)=-1 then
 Add(L,t^-1);
fi;
od;

if IsEvenInt(Length(L)) then

if A[1][1]=-1 then
A:=S^2*A;
Add(L,s^2);
fi;

else

if A[1][1]=-1 then
A:=S^2*A;
Add(L,s^2);
fi;

fi;

#Print(A,"\n");
return Product(List(L,x->x^-1));
end;
#########################################################
#########################################################

gensH:=List(GeneratorsOfGroup(H),x->ElementToWord(x));
gensH:=List(gensH,x->Image(epi,x));
QH:=Subgroup(Q,gensH);

####################
if HH=false then

QH!.epimorphism:=epi;
QH!.ElementToWord:=ElementToWord;
return QH;

iso:=IsomorphismFpGroup(QH);
QH!.isoFpGroup:=GroupHomomorphismByFunction(H, Image(iso), x->
                Image(iso,Image(epi,ElementToWord(x)) ) );

fi;
####################

gensHH:=List(GeneratorsOfGroup(HH),x->ElementToWord(x));
gensHH:=List(gensHH,x->Image(epi,x));
QHH:=Subgroup(QH,gensHH);
R1:=RightTransversal(QH,QHH);;
R2:=List(R1,x->PreImagesRepresentative(epi,x));
Apply(R2,x->Image(FhomG,x));

#####################################################
#####################################################
poscan:=function(g)
local w, a;

w:=ElementToWord(g);
w:=Image(epi,w);
return PositionCanonical(R1,w);

end;
#####################################################
#####################################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := H,
               subgroup := HH,
               cosets:=R2,
               poscan:=poscan ));

end);;
###########################################################
###########################################################

