
############################################################
############################################################
InstallGlobalFunction(AsWordInSL2Z,
function(g)
local  S, T, G,  ElementToWord, s, t;

G:=SL(2,Integers);

S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
s:=S;
t:=T;

###########################################################
###########################################################
ElementToWord:=function(g)
local A, d,i, b,L, bool;

#bool is true if S^n was the last added element
bool:=false;


L:=[IdentityMat(2)];
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
return List(L,x->x^-1);
end;
#########################################################
#########################################################

return ElementToWord(g);
end);
###########################################################
###########################################################
