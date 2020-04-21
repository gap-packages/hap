
############################################################
############################################################
InstallMethod(RightTransversal,
"right transversal for finite index subgroups of SL(2,Integers)",
[IsMatrixGroup,IsMatrixGroup],
function(H,HH)
local F, rels, S, T, G, FhomG, Q, gensQ, epi, ElementToWord,
gensH, gensHH, QH, QHH, R, R1, R2, s, t, poscan;

G:=SL(2,Integers);
if not IsSubgroup(G,H) and IsSubgroup(H,HH) then
TryNextMethod(); fi;

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
gensHH:=List(GeneratorsOfGroup(HH),x->ElementToWord(x));
gensHH:=List(gensHH,x->Image(epi,x));
#gensH:=Concatenation(gensH,gensHH);
QH:=Subgroup(Q,gensH);
QHH:=Subgroup(QH,gensHH);
R1:=RightTransversal(QH,QHH);;
R2:=List(R1,x->PreImagesRepresentative(epi,x));
Apply(R2,x->Image(FhomG,x));

#####################################################
#####################################################
poscan:=function(g)
local w;
w:=ElementToWord(g);
w:=Image(epi,w);
return PositionCanonical(R1,w);
end;
#####################################################
#####################################################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2Subgroup and IsList and
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
"for HapRightTransversals of subrougs in SL(2,Z)",
[IsHapRightTransversalSL2Subgroup,IsObject],
function(R,x)
return R!.poscan(x);
end);
###########################################################
###########################################################
