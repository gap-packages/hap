#(C) Graham Ellis 2009

#########################################################
#########################################################
InstallGlobalFunction(RipsChainComplex,
function(arg)
local 
	M,epsilon,BOOL,
	BinaryMatrix,
        FindFreeEdge,
	MaximalTree,
	MaximalForest,
	ExtendedForest,
	SelectTree,
	ChnCmplx,
	colour,
	MM,	
	bool,
	Edges,
	Faces,
	Dimension,
	Boundary,
	answer,
	zro, zroone,
        BettiZero, BettiZeroA, BettiZeroB,
	i,j,k,c;

M:=arg[1];
epsilon:=arg[2];
if Length(arg)=3 then BOOL:=true; else BOOL:=false; fi;

################################################
BinaryMatrix:=function(M,t)
local i,j, len;
#M is a square symmetric matrix and t a number.
#The matrix M will be destroyed: it will become a binary matrix.

len:=Length(M);
for i in [1..len] do
for j in [1..len] do
if M[i][j]>t then M[i][j]:=0;
else M[i][j]:=1;
fi;
od;
od;

for i in [1..len] do
M[i][i]:=0;
od;

end; 
################################################

################################################
FindFreeEdge:=function(M)
local i,j,len;

len:=Length(M);
for i in [1..len] do
for j in [i..len] do
if M[i][j]=1 then 
if Maximum(M[i])<=1 and Maximum(M[j])<=1 then
return [i,j]; 
fi;
fi;
od;
od;

return false;
end;
################################################

################################################
MaximalTree:=function(M,v)
local i,j,colour,len,leaves,newleaves,vertices;
#Find a maximal tree containing the free edge v.
#Colour the edges of the tree.

len:=Length(M);
colour:=Maximum(List(M,x->Maximum(x)))+1;
M[v[1]][v[2]]:=colour;
M[v[2]][v[1]]:=colour;
leaves:=v;
vertices:=SSortedList(v);


while Size(leaves)>0 do
newleaves:=[];
for i in leaves do
for j in [1..len] do
#if M[i][j]=1 and not colour in M[j] then 
if M[i][j]=1 and not j in vertices then
M[i][j]:=colour; M[j][i]:=colour; AddSet(vertices,j);
Add(newleaves,j);
fi;
od;
od;
leaves:=newleaves;
od;

return colour;
end;
################################################

################################################
MaximalForest:=function(M)
local v,T,colour;

v:=FindFreeEdge(M);
if v=false then return 0; fi;

while not v=false do
colour:=MaximalTree(M,v);
v:=FindFreeEdge(M);
od;

return colour-1;
end;
################################################

################################################
ExtendedForest:=function(M)
local colour,i,j,k,Triangle,bool;

colour:=MaximalForest(M)+2;

####################
Triangle:=function(i,j,k)
local lst;

lst:=[ M[i][j], M[i][k], M[j][k] ];
if 0 in lst then return false; fi;
if not Length(Filtered(lst,a->a=1))=1 then return false; fi;
if lst[1]=1 then M[i][j]:=colour; M[j][i]:=colour; fi;
if lst[2]=1 then M[i][k]:=colour; M[k][i]:=colour; fi;
if lst[3]=1 then M[j][k]:=colour; M[k][j]:=colour; fi;
return true;
end;
####################

bool:=false;

while bool do
bool:=false;
for i in [1..Length(M)] do
for j in [i..Length(M)] do
for k in [j..Length(M)] do
if Triangle(i,j,k) then bool:=true; fi;
od;od;od;
od;


return colour-2;

end;
################################################

MM:=StructuralCopy(M);
BinaryMatrix(MM,epsilon);
BettiZeroA:=ExtendedForest(MM);
BettiZeroB:=0;
for i in [1..Length(MM)] do
if Sum(MM[i])=0 then BettiZeroB:=BettiZeroB+1; fi;
od;
BettiZero:=BettiZeroA + BettiZeroB;

################################################
SelectTree:=function(MM,c)
local M,i,j;
M:=StructuralCopy(MM);

for i in [1..Length(M)] do
for j in [i..Length(M)] do
if not M[i][j]=0 then
if not c=M[i][j] and (not c in M[i]) and not (c in M[j]) then
M[i][j]:=0; M[j][i]:=0;
fi;
fi;
od;
od;

return M;
end;
################################################

#####################################
ChnCmplx:=function(MM,c)
local  
        bool,
        Edges,
        Faces,
        Dimension,
        Boundary,
        zro, zroone,
        i,j,k;
 

Edges:=[];
for i in [1..Length(MM)] do
for j in [i+1..Length(MM)] do
if MM[i][j]=1 then Add(Edges,[i,j]); fi;
od;
od;
Edges:=SSortedList(Edges);
Faces:=[];

for i in [1..Length(MM)] do
for j in [i..Length(MM)] do
for k in [j..Length(MM)] do
if (not MM[i][j]=0) and (not MM[i][k]=0) and (not MM[j][k]=0) then
if MM[i][j]=1 or MM[i][k]=1 or MM[j][k]=1  then Add(Faces,[i,j,k]); fi;
fi;
od;od;od;


#################################
Dimension:=function(n);
if n=0 and c=1 then return BettiZero; fi;
if n=0 then return 0; fi;
if n=1 then return Length(Edges); fi;
if n=2 then return Length(Faces); fi;
if n>2 then return 0; fi;
end;
#################################

zro:=List([1..BettiZero],i->0);
zroone:=List([1..Length(Edges)],i->0);
#################################
Boundary:=function(n,i)
local a,x,y,z;
if n=2 then 
  a:=zroone*0;
  x:=Position(Edges,[Faces[i][1],Faces[i][2]]); 
  if not x=fail then a[x]:=1; fi;
  y:=Position(Edges,[Faces[i][2],Faces[i][3]]);
  if not y=fail then a[y]:=1; fi;
  z:=Position(Edges,[Faces[i][1],Faces[i][3]]);
  if not z=fail then a[z]:=-1; fi;
  return a;
fi;
if n=1 then return zro; fi;
if n=0 then return [0]; fi;
end;
#################################

return
 Objectify(HapChainComplex,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=
                [["length",3],
                ["connected",true],
                ["type", "chainComplex"],
                ["characteristic",0] 
   		]));
end;
#####################################

if BOOL=false then 
return ChnCmplx(MM,1);
fi;

answer:=[];
for c in [1..BettiZeroA] do
Add(answer,ChnCmplx(SelectTree(MM,c+1),c));
od;

return answer;
end);
####################################################
####################################################







###########################################
###########################################
ReadBioData:=function(file)
local L,M,fun,T,t,it;

it:=InputTextFile(file);
Sleep(1);
L:=ReadLine(it);;
M:=[];

#################
fun:=function(x);
if x='\t' then return ',';
else return x; fi;
end;
#################

while true do
L:=ReadLine(it);;
if L=fail then break; fi;
t:=Position(L,'\t');
L:=L{[t..Length(L)-2]};
T:=List(L,fun);
T[1]:='[';
Add(T,']');
T:=EvalString(T);
Add(M,T);
od;


return M;

end;
#######################################
#######################################


#######################################
#######################################
InstallGlobalFunction(VectorsToSymmetricMatrix,
function(arg)
local M, Distance,S,i,j;

#############################
M:=arg[1];
if Length(arg)=1 then
	Distance:=function(u,v);
	return Sum(List(u-v,x->AbsInt(x)));
	end;
else Distance:=arg[2];
fi;
#############################

S:=[];
for i in [1..Length(M)] do
S[i]:=[];
for j in [1..Length(M)] do
S[i][j]:=Distance(M[i],M[j]);
od;
od;

return S;
end);
#######################################
#######################################
