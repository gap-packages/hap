####################################################
####################################################
InstallGlobalFunction(SparseMat,
function(M)
local S, Rows, Cols, char, i, j;

Rows:=Length(M);
if Rows=0 then Cols:=0;
else
Cols:=Length(M[1]);;
fi;

if Rows=0 then char:=0;
else
char:=Characteristic(M[1][1]);
fi;

S:=List([1..Rows],i->[]);

for i in [1..Rows] do
for j in [1..Cols] do
if not IsZero(M[i][j]) then
Add(S[i],[j,M[i][j]]);
fi;
od;
od;


return Objectify(HapSparseMat,
               rec( rows:=Rows,
                    cols:=Cols,
                    characteristic:=char,
                    mat:=S));
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseRowMult,
function(M,i,k)
local x;

for x in M!.mat[i] do
x[2]:=k*x[2];
od;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseRowInterchange,
function(M,i,j)
local x;

x:=M!.mat[i];
M!.mat[i]:=M!.mat[j];
M!.mat[j]:=x;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseRowAdd,
function(M,i,j,k)
local x,y,z,r;

y:=List(M!.mat[i],a->a[1]);
z:=[];
for r in [1..Length(y)] do
z[y[r]]:=r;
od;

for x in M!.mat[j] do
if IsBound(z[x[1]]) then 
M!.mat[i][z[x[1]]][2]:= M!.mat[i][z[x[1]]][2]+k*x[2];
else
Add(M!.mat[i],[x[1],k*x[2]]);
fi;
od;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparsePivot,
function(A,p)
local i;

if not IsBound(A!.heads) then A!.heads:=[]; fi;

if not IsBound(A!.heads[p]) then
for i in [1..A!.rows] do
A!.mat[i]:=SSortedList(Filtered(A!.mat[i],x->not x[2]=0));
if IsBound(A!.mat[i][1]) then
if A!.mat[i][1][1]=p then A!.heads[p]:=i; return i; fi;
fi;
od;
A!.heads[p]:=fail;
fi;

return A!.heads[p];;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparsePivots,
function(A)
local i,x;

if not IsBound(A!.heads) then A!.heads:=[]; fi;

for i in [1..A!.rows] do
A!.mat[i]:=SSortedList(Filtered(A!.mat[i],x->not x[2]=0));
if Length(A!.mat[i])>0 then
 A!.heads[A!.mat[i][1][1]]:=i;  
fi;
od;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SemiEchelonSMat,
function(A)
local remrows,i,tmp;

SparsePivots(A);

###
remrows:=Difference([1..A!.rows],A!.heads);
tmp:=[];
for i in remrows do
if Length(A!.mat[i])=0 then Add(tmp,i);fi; 
od;
remrows:=Difference(remrows,tmp);
###

return remrows;
end);
####################################################
####################################################

