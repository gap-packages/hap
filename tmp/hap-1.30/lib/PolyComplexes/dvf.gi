######################################################
######################################################
InstallGlobalFunction(ContractCubicalComplex,
function(K);

if not IsHapCubicalComplex(K) then 
Print("This function must be applied to a cubical complex.\n");
fi;

if not IsBound(K!.vectors) then
K!.vectors:=0*K!.binaryArray;
fi;

if Dimension(K)=2 then
ContractCubicalComplex_dim2(K);
fi;

if Dimension(K)=3 then
ContractCubicalComplex_dim3(K);
fi;

end);
######################################################
######################################################

######################################################
######################################################
InstallGlobalFunction(ContractCubicalComplex_dim2,
function(K)
local A,IsFree,i,j,Rows,RRows,Cols,RCols,Toggle,Vectors,Rewrite,FinalRewrite,pos,L,v;

if not Dimension(K)=2 then 
Print("This function works only for 2-dimensional cubical complexes.\n");
return fail;
fi;


A:=K!.binaryArray;
A:=FrameArray(A);
Vectors:=K!.vectors;

#####################
#####################
IsFree:=function(i,j)
local S;

###
if A[i][j]=0 then return false; fi;
###

###
if IsOddInt(i) and IsOddInt(j) then
return false;
fi;
###

###
if IsEvenInt(i) and IsEvenInt(j) then
S:=A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1];
if S=1 then
   pos:=Position([A[i-1][j],A[i+1][j],A[i][j-1],A[i][j+1]],1);
   L:=[[i-1,j],[i+1,j],[i,j-1],[i,j+1]];
   v:=L[pos];
   Vectors[i-1][j-1]:=[v[1]-1,v[2]-1];
   A[i][j]:=0;
   A[v[1]][v[2]]:=0;
   return true;
else return false;
fi;
fi;
###

###
if IsEvenInt(i) and IsOddInt(j) then
S:=A[i-1][j]+A[i+1][j];
if S=1 then
   pos:=Position([A[i-1][j],A[i+1][j]],1);
   L:=[[i-1,j],[i+1,j]];
   v:=L[pos];
   Vectors[i-1][j-1]:=[v[1]-1,v[2]-1];
   A[i][j]:=0;
   A[v[1]][v[2]]:=0;
   return true;
else return false;
fi;
fi;
###

###
if IsOddInt(i) and IsEvenInt(j) then
S:=A[i][j-1]+A[i][j+1];
if S=1 then
   pos:=Position([A[i][j-1],A[i][j+1]],1);
   L:=[[i,j-1],[i,j+1]];
   v:=L[pos];
   Vectors[i-1][j-1]:=[v[1]-1,v[2]-1];
   A[i][j]:=0;
   A[v[1]][v[2]]:=0;
   return true;
else return false;
fi;
fi;
###


end;
#####################
#####################


#####################
Toggle:=true;
while Toggle do
Toggle:=false;

Rows:=Filtered([1..Length(A)],i->1 in A[i]);
RRows:=Reversed(Rows);
Cols:=[2..Length(A[1])-1];
RCols:=Reversed(Cols);

for i in Rows do
for j in Cols do
if IsFree(i,j) then Toggle:=true; fi;;
od;od;
for i in RRows do
for j in RCols do
if IsFree(i,j) then Toggle:=true; fi;;
od;od;
od;
#####################

A:=UnframeArray(A);
K!.binaryArray:=A;

#####################
Rewrite:=function(wrd)
local  rewrt, a, A, Bnd;

####
Bnd:=function(v)
local i;

if IsEvenInt(v[1]) and IsEvenInt(v[2]) then
return [[v[1]-1,v[2]],[v[1]+1,v[2]],[v[1],v[2]-1],[v[1],v[2]+1]];
fi;
if IsEvenInt(v[1]) and IsOddInt(v[2]) then
return [[v[1]-1,v[2]],[v[1]+1,v[2]]];
fi;
if IsOddInt(v[1]) and IsEvenInt(v[2]) then
return [[v[1],v[2]-1],[v[1],v[2]+1]];
fi;

end;
####

A:=K!.binaryArray;
rewrt:=[];
for a in wrd do
if A[a[1]][a[2]]=1 then Add(rewrt,a);  fi;
if A[a[1]][a[2]]=0 and IsList(K!.vectors[a[1]][a[2]]) then
 v:=Vectors;
 v:=v[a[1]][a[2]];
 v:=Difference(Bnd(v),[a]); 
 Append(rewrt,Rewrite(v));
fi;
od;

return rewrt;
end;
#####################

K!.rewrite:=Rewrite;

end);
####################################################
####################################################

######################################################
######################################################
InstallGlobalFunction(ContractCubicalComplex_dim3,
function(K)
local A,IsFree,i,j,k,Rows,RRows,Cols,RCols,Pages,RPages,Toggle,Vectors,v,pos,L,Rewrite;

if not Dimension(K)=3 then
Print("This function works only for 3-dimensional cubical complexes.\n");
return fail;
fi;


A:=K!.binaryArray;
A:=FrameArray(A);
Vectors:=K!.vectors;

#####################
#####################
IsFree:=function(i,j,k)
local S;

###
if A[i][j][k]=0 then return false; fi;
###

### 3-cells
if IsOddInt(i) and IsOddInt(j) and IsOddInt(k) then
return false;
fi;
###

### 0-cells
if IsEvenInt(i) and IsEvenInt(j) and IsEvenInt(k) then
S:=A[i-1][j][k]+A[i+1][j][k]+A[i][j-1][k]+A[i][j+1][k] + A[i][j][k-1]+A[i][j][k+1] ;
if S=1 then
   pos:=Position([A[i-1][j][k],A[i+1][j][k],A[i][j-1][k],A[i][j+1][k],
                  A[i][j][k-1],A[i][j][k+1]],1);
   L:=[[i-1,j,k],[i+1,j,k],[i,j-1,k],[i,j+1,k],[i,j,k-1],[i,j,k+1]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

##
#if S=1 then
#A[i][j][k]:=0;
#A[i-1][j][k]:=0;A[i+1][j][k]:=0;A[i][j-1][k]:=0;A[i][j+1][k]:=0; A[i][j][k-1]:=0;A[i][j][k+1]:=0;
#return true;
#else return false;
#fi;
##
fi;
###

### 1-cells
if IsEvenInt(i) and IsOddInt(j) and IsEvenInt(k) then
S:=A[i-1][j][k]+A[i+1][j][k]+A[i][j][k-1]+A[i][j][k+1];

if S=1 then
   pos:=Position([A[i-1][j][k],A[i+1][j][k],
                  A[i][j][k-1],A[i][j][k+1]],1);
   L:=[[i-1,j,k],[i+1,j,k],[i,j,k-1],[i,j,k+1]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

#if S=1 then
#A[i][j][k]:=0;
#A[i-1][j][k]:=0;A[i+1][j][k]:=0;A[i][j][k-1]:=0;A[i][j][k+1]:=0;
#return true;
#else return false;
#fi;
fi;
###

### 1-cells
if IsEvenInt(i) and IsEvenInt(j) and IsOddInt(k) then
S:=A[i-1][j][k]+A[i+1][j][k]+A[i][j-1][k]+A[i][j+1][k];

if S=1 then
   pos:=Position([A[i-1][j][k],A[i+1][j][k],A[i][j-1][k],A[i][j+1][k],
                  ],1);
   L:=[[i-1,j,k],[i+1,j,k],[i,j-1,k],[i,j+1,k]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;



#if S=1 then
#A[i][j][k]:=0;
#A[i-1][j][k]:=0;A[i+1][j][k]:=0;A[i][j-1][k]:=0;A[i][j+1][k]:=0;
#return true;
#else return false;
#fi;
fi;
###

### 1-cells
if IsOddInt(i) and IsEvenInt(j) and IsEvenInt(k) then
S:=A[i][j-1][k]+A[i][j+1][k]+A[i][j][k-1]+A[i][j][k+1];

if S=1 then
   pos:=Position([A[i][j-1][k],A[i][j+1][k],
                  A[i][j][k-1],A[i][j][k+1]],1);
   L:=[[i,j-1,k],[i,j+1,k],[i,j,k-1],[i,j,k+1]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

#if S=1 then
#A[i][j][k]:=0;
#A[i][j-1][k]:=0;A[i][j+1][k]:=0;A[i][j][k-1]:=0;A[i][j][k+1]:=0;
#return true;
#else return false;
#fi;
fi;
###

### 2-cells
if IsOddInt(i) and IsOddInt(j) and IsEvenInt(k) then
S:=A[i][j][k-1]+A[i][j][k+1];

if S=1 then
   pos:=Position([ A[i][j][k-1],A[i][j][k+1]],1);
   L:=[[i,j,k-1],[i,j,k+1]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

#if S=1 then
#A[i][j][k]:=0;
#A[i][j][k-1]:=0;A[i][j][k+1]:=0;
#return true;
#else return false;
#fi;

fi;
###

### 2-cells
if IsOddInt(i) and IsEvenInt(j) and IsOddInt(k) then
S:=A[i][j-1][k]+A[i][j+1][k];

if S=1 then
   pos:=Position([A[i][j-1][k],A[i][j+1][k] ],1);
   L:=[[i,j-1,k],[i,j+1,k]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

#if S=1 then
#A[i][j][k]:=0;
#A[i][j-1][k]:=0;A[i][j+1][k]:=0;
#return true;
#else return false;
#fi;
fi;
###

### 2-cells
if IsEvenInt(i) and IsOddInt(j) and IsOddInt(k) then
S:=A[i-1][j][k]+A[i+1][j][k];

if S=1 then
   pos:=Position([A[i-1][j][k],A[i+1][j][k]],1);
   L:=[[i-1,j,k],[i+1,j,k]];
   v:=L[pos];
   Vectors[i-1][j-1][k-1]:=[v[1]-1,v[2]-1,v[3]-1];
   A[i][j][k]:=0;
   A[v[1]][v[2]][v[3]]:=0;
   return true;
else return false;
fi;

#if S=1 then
#A[i][j][k]:=0;
#A[i-1][j][k]:=0;A[i+1][j][k]:=0;
#return true;
#else return false;
#fi;

fi;
###

end;
#####################
#####################


#####################
Toggle:=true;
while Toggle do
Toggle:=false;

Rows:=Filtered([1..Length(A)],i->1 in Flat(A[i]));
RRows:=Reversed(Rows);
Cols:=[2..Length(A[1])-1];
RCols:=Reversed(Cols);
Pages:=[2..Length(A[1][1])-1];
RPages:=Reversed(Pages);


for i in Rows do
for j in Cols do
for k in Pages do
if IsFree(i,j,k) then Toggle:=true; fi;;
od;od;od;
for i in RRows do
for j in RCols do
for k in RPages do
if IsFree(i,j,k) then Toggle:=true; fi;;
od;od;od;
od;
#####################

A:=UnframeArray(A);
K!.binaryArray:=A;

#####################
Rewrite:=function(wrd)
local  rewrt, a, A, Bnd;

####
Bnd:=function(v)
local i;

if IsEvenInt(v[1]) and IsEvenInt(v[2]) and IsEvenInt(v[3]) then
return [[v[1]-1,v[2],v[3]],[v[1]+1,v[2],v[3]],[v[1],v[2]-1,v[3]],
        [v[1],v[2]+1,v[3]],[v[1],v[2],v[3]-1],[v[1],v[2],v[3]+1]];
fi;
if IsEvenInt(v[1]) and IsOddInt(v[2]) and IsEvenInt(v[3]) then
return [[v[1]-1,v[2],v[3]],[v[1]+1,v[2],v[3]],[v[1],v[2],v[3]-1],
        [v[1],v[2],v[3]+1]];
fi;
if IsOddInt(v[1]) and IsEvenInt(v[2]) and IsEvenInt(v[3]) then
return [[v[1],v[2]-1,v[3]],[v[1],v[2]+1,v[3]], [v[1],v[2],v[3]-1],
        [v[1],v[2],v[3]+1]];
fi;

if IsEvenInt(v[1]) and IsEvenInt(v[2]) and IsOddInt(v[3]) then
return [[v[1],v[2]-1,v[3]],[v[1],v[2]+1,v[3]], [v[1]-1,v[2],v[3]],
        [v[1]+1,v[2],v[3]]];
fi;

if IsEvenInt(v[1]) and IsOddInt(v[2]) and IsOddInt(v[3]) then
return [[v[1]-1,v[2],v[3]], [v[1]+1,v[2],v[3]]];
fi;

if IsOddInt(v[1]) and IsEvenInt(v[2]) and IsOddInt(v[3]) then
return [[v[1],v[2]-1,v[3]], [v[1],v[2]+1,v[3]]];
fi;

if IsOddInt(v[1]) and IsOddInt(v[2]) and IsEvenInt(v[3]) then
return [[v[1],v[2],v[3]-1], [v[1],v[2],v[3]+1]];
fi;

end;
####

A:=K!.binaryArray;
rewrt:=[];
for a in wrd do
if A[a[1]][a[2]][a[3]]=1 then Add(rewrt,a);  fi;
if A[a[1]][a[2]][a[3]]=0 and IsList(K!.vectors[a[1]][a[2]][a[3]]) then
 v:=Vectors;
 v:=v[a[1]][a[2]][a[3]];
 v:=Difference(Bnd(v),[a]);
 Append(rewrt,Rewrite(v));
fi;
od;

return rewrt;
end;
#####################

K!.rewrite:=Rewrite;

end);
####################################################
####################################################


####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(DVFReducedCubicalComplex,
function(YY)
local
        ChooseCriticalCell, B, N, S, SS, Y;

if not IsHapCubicalComplex(YY) then
Print("This function must be applied to a cubical complex.\n");
return fail;
fi;

Y:=Objectify(HapCubicalComplex,rec( ));;
Y!.binaryArray:=StructuralCopy(YY!.binaryArray);
Y!.properties:=StructuralCopy(YY!.properties);
B:=0*Y!.binaryArray;


######################################
ChooseCriticalCell:=function(N)
local i,j,k,A,toggle;

A:=Y!.binaryArray;
toggle:=false;

if Dimension(Y)=2 then
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
if A[i][j]=1 then
if Length(Filtered([i,j],x->IsEvenInt(x)))=N then
toggle:=true; A[i][j]:=0;B[i][j]:=1;break; fi;
fi;
od; if toggle then break; fi;
od;
fi;


if Dimension(Y)=3 then
for i in [1..Length(A)] do
for j in [1..Length(A[1])] do
for k in [1..Length(A[1][1])] do
if A[i][j][k]=1 then
if Length(Filtered([i,j,k],x->IsEvenInt(x)))=N then
toggle:=true; A[i][j][k]:=0;B[i][j][k]:=1;break; fi;
fi;
od; if toggle then break; fi;
od; if toggle then break; fi;
od;
fi;

Y!.binaryArray:=A;
end;
#######################################

ContractCubicalComplex(Y);
N:=Dimension(Y);
S:=ArraySum(Y!.binaryArray);
while S>0 do

ChooseCriticalCell(N);
ContractCubicalComplex(Y);
SS:=ArraySum(Y!.binaryArray);
if S=SS  then if N>0 then N:=N-1; else break; fi;fi; #N:=Dimension(Y);fi; fi;
S:=SS;
od;

Y!.binaryArray:=B;
Add(Y!.properties,["nonregular",true]);
return Y;
end);
######################################################
######################################################





