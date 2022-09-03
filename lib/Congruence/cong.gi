
############################################################
############################################################
InstallMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsMatrixGroup,IsMatrixGroup],
function(H,HH);
if not (IsHapSL2ZSubgroup(HH) or IsHapSL2ZSubgroup(H)) then TryNextMethod(); fi;
#Print("HAP_RightTransversalSL2ZSubgroups\n");
return HAP_RightTransversalSL2ZSubgroups(H,HH);
end);
############################################################
############################################################

############################################################
############################################################
InstallGlobalFunction(HAP_RightTransversalSL2ZSubgroups,
function(H,HH)
local F, rels, S, T, G, FhomG, Q, gensQ, epi, ElementToWord,
gensH, gensHH, QH, QHH, R, R1, R2, s, t, poscan, iso;

#Print("start\n");
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

#Print("stop\n");
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

###########################################################
###########################################################
InstallOtherMethod(PositionCanonical,
"for HapRightTransversals of subrougs in SL(2,Z) or SL(2,O)",
[IsHapRightTransversalSL2ZSubgroup,IsObject],
function(R,x)
return R!.poscan(x);
end);
###########################################################
###########################################################



################################################
################################################
InstallOtherMethod(IndexNC,
"index for HapSLOSubgroups",
[IsMatrixGroup,IsHapSL2ZSubgroup],
function(G,H);
return Length(RightTransversal(G,H));
end);
################################################
################################################

################################################
################################################
InstallOtherMethod(IndexNC,
"index for HapSLOSubgroups",
[IsMatrixGroup,IsHapSL2OSubgroup],
function(G,H);
return Length(RightTransversal(G,H));
end);
################################################
################################################


