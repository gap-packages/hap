###########################################################
###########################################################
InstallGlobalFunction(ElementsLazy,
function(G)
local fn, L, Elts, posfn,LT,mult,multrec,inv,invrec;

L:=[];
LT:=NewDictionary(GeneratorsOfGroup(G)[1],true,G);
#############################################
fn:=function(i);
if i>0 and i<=Length(L) then return L[i]; fi;
return fail;
end;
#############################################
posfn:=function(x)
local pos;
if not x in G then return fail; fi;
pos:=LookupDictionary(LT,x);

if pos=fail then
Add(L,x);
pos:=Length(L);
AddDictionary(LT,x,pos);
fi;
return pos;

end;

Elts:=LazyList(fn);
Elts!.posfun:=posfn;

multrec:=[];
#############################################
mult:=function(i,j);
if not IsBound(multrec[i]) then multrec[i]:=[]; fi;
if not IsBound(multrec[i][j]) then
multrec[i][j]:=posfn(Elts[i]*Elts[j]);
fi;
return multrec[i][j];
end;
#############################################

invrec:=[];
#############################################
inv:=function(g);
if not IsBound(invrec[g]) then
invrec[g]:= posfn(Elts[g]^-1);
fi;
return invrec[g];
end;
#############################################

Elts!.mult:=mult;
Elts!.inv:=inv;

return Elts;

end);
###########################################################
###########################################################

###########################################################
###########################################################
InstallGlobalFunction(HAP_ElementsSL2Zfn,
function()
local DecimalToBinary,
      BinaryToDecimal,
      ElementToDecimal,
      DecimalToElement,
      InSL2Z,
      Elts;
      

###########################################################
###########################################################
DecimalToBinary:=function(DD)
local D,L,log, n,i,m;
if DD<0 then return fail; fi;
if DD=0 then return [0]; fi;
D:=1*DD;
log:=Log(D,2);
n:=log; 
L:=[];
for i in [0..n] do 
m:=2^(n-i);
if m<=D then Add(L,1);
D:=D-m;
else Add(L,0); 
fi;
od;

return L;
end;
###########################################################
###########################################################


###########################################################
###########################################################
BinaryToDecimal:=function(L)
local D,ln,i;
ln:=Length(L);
D:=0;
for i in [1..ln] do
if L[i]>0 then
D:=D+ 2^(ln-i);
fi;
od;
return D;

end;
###########################################################
###########################################################

###########################################################
###########################################################
InSL2Z:=function(g);
if not IsList(g) then return false; fi;
if not Length(g)=2 then return false; fi;
if not IsList(g[1]) and IsList(g[2]) then return false; fi;
if not Length(g[1])=2 and Length(g[2])=2 then return false; fi;
if IsInt(g[1][1])
and IsInt(g[1][2])
and IsInt(g[2][1])
and IsInt(g[2][2])
and g[1][1]*g[2][2]-g[2][1]*g[1][2] = 1
then return true; 
else return false; fi;
end;
###########################################################
###########################################################


###########################################################
###########################################################
ElementToDecimal:=function(g)
local S,T,A, d,i, b,L, bool;

#bool is true if S^n was the last added element
bool:=false;

if not InSL2Z(g) then
return fail;
fi;

L:=[];
S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
A:=1*g; 

while not A[2][1]=0 do
if AbsInt(A[2][1])>AbsInt(A[1][1]) then
   A:=S*A; Add(L,1);  bool:=true;
   else Add(L,0);  bool:=true;
fi;
if  not A[2][1]=0 then 
  if A[1][1]/A[2][1]>=0 then
    A:=T^-1*A; 
    Add(L,1);  bool:=false;
  else
    A:=T*A;  bool:=false;
    Add(L,0);
  fi;
fi;
od;




d:=-A[1][1]*A[1][2];


if not bool then
   A:=T^SignInt(d)*A;

   if SignInt(d)=1 then
      Add(L,0); Add(L,0);  
   fi;
   if SignInt(d)=-1 then
      Add(L,0); Add(L,1); 
   fi;
else
   A:=T^SignInt(d)*A;
   bool:=false;
   if SignInt(d)=1 then
       Add(L,0);  
   fi;
   if SignInt(d)=-1 then
       Add(L,1);  
   fi;
fi;

for i in [2..AbsInt(d)] do
A:=T^SignInt(d)*A;

if SignInt(d)=1 then
Add(L,0); Add(L,0); 
fi;
if SignInt(d)=-1 then
 Add(L,0); Add(L,1); 
fi;
od;

if IsEvenInt(Length(L)) then

if A[1][1]=-1 then 
A:=S^2*A;  
 L[1+Length(L)]:= 2; 
else
L[1+Length(L)]:= 0; 
fi;

else

if A[1][1]=-1 then
A:=S^2*A;  
 L[Length(L)]:= L[Length(L)]+2; 
fi;

fi;

b:=DecimalToBinary(L[Length(L)]);
if Length(b)=1 then b[2]:=b[1]; b[1]:=0; fi;
b:=Reversed(b);
L[Length(L)]:=b[1];
L[1+Length(L)]:=b[2];
Add(L,1);  #just to record the length of L

L:=Reversed(L);

return BinaryToDecimal(L);
end;
###########################################################
###########################################################

###########################################################
###########################################################
DecimalToElement :=function(n)
local g, b, S, T, A, i;

if n<=0 then return 0; fi;
if n=1 then return [[1,0],[0,1]]; fi;

S:=[[0,-1],[1,0]];
T:=[[1,1],[0,1]];
A:=[[1,0],[0,1]];;

b:=DecimalToBinary(n);
if IsEvenInt(Length(b)) then 
b:=Reversed(b); Add(b,1); b:=Reversed(b);
fi;

A:=S^BinaryToDecimal([b[2],b[3]]);

for i in [2..(Length(b)-1)/2] do
if b[2*i]=0 then
     A:=A*T; else A:=A*T^-1;
fi;
if b[2*i+1] =1 then
     A:=A*S;
fi;
od;

return A^-1;
end;
###########################################################
###########################################################

Elts:=LazyList(DecimalToElement);
Elts!.posfun:=ElementToDecimal;

return Elts;
end);

###########################################################
###########################################################
InstallGlobalFunction(ResolutionSL2Z_alt,
function(n)
local C,R;

C:=ContractibleSL2ZComplex_alt();
R:=ResolutionGTree(C,n);

return R;
end);
###########################################################
###########################################################
