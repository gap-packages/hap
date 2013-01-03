#(C) Graham Ellis

####################################################
####################################################
InstallGlobalFunction(SparseMat,
function(M)
local S, Rows, Cols, char, i, j;

#Inputs a matrix M and returns a sparse matrix

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

#Multiplies the i-th row of a sparse matrix by k.

if IsZero(k) then return false; fi;

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

#Interchanges the i-th and j-th rows of a sparse matrix M.

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

#Adds k times the j-th row to the i-th row of a sparse matrix M.

y:=StructuralCopy(M!.mat[i]);
Apply(y,a->a[1]);
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

M!.mat[i]:=Filtered(M!.mat[i],a->not IsZero(a[2]));

Sort(M!.mat[i]);

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseSemiEchelon,
function(M)
local i, V, fn;

if not IsHapSparseMat(M) then
Print("This function must be applied to a HAP sparse matrix.\n");
return fail; fi;

if not IsBound(M!.heads) then M!.heads:=[]; fi;

# i = M!.heads[j] is the number of the row whose first non-zero entry
# is in column j.

V:=[1..M!.rows];
##########
fn:=function(i,j) 
local v; 
if IsBound(M!.mat[i]) and IsBound(M!.mat[j]) then return 
Maximum(List(M!.mat[i],x->x[1])) 
>
Maximum(List(M!.mat[j],x->x[1]));
fi;
return false;
end; 
#########
Sort(V,fn);

for i in V do
if IsBound(M!.mat[i]) then
   if not IsSSortedList(M!.mat[i]) then
   Sort(M!.mat[i]);
   fi;
   if IsBound(M!.mat[i][1]) then
   if not IsBound(M!.heads[M!.mat[i][1][1]]) and not IsZero(M!.mat[i][1][2]) then
   M!.heads[M!.mat[i][1][1]]:=i;
      if not IsOne(M!.mat[i][1][2])then 
      SparseRowMult(M,i,M!.mat[i][1][2]^-1);
      fi;
   fi;
   fi;
fi;
od;


for i in Difference([1..M!.rows],M!.heads) do
SparseRowReduce(M,i);
od;

#M!.mat:=Filtered(M!.mat,r->not Length(r)=0);
for i in M!.mat do
if Length(i)=0 then Unbind(i); fi;
od;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseRowReduce,
function(M,i)
local x, first, r,k;

#Row reduce the i-th row of a sparse matrix if it is non-empty. 
#The first non-zero entry of the reduced row will be 1.
#Return true if the matrix is modified, and false otherwise.

if not IsBound(M!.mat[i]) then return false; fi;
if Length(M!.mat[i])=0 then return false; fi;

first:=M!.mat[i][1][1];

   while IsBound(M!.heads[first]) do
   r:=M!.heads[first]; 
   k:=-M!.mat[i][1][2];
   SparseRowAdd(M,i,r,k);
   if Length(M!.mat[i])=0 then return true;  fi;
   first:=M!.mat[i][1][1];
   od;

k:=M!.mat[i][1][2];
if not IsOne(k) then
SparseRowMult(M,i,k^-1);
fi;
M!.heads[first]:=i;
return true;


end);
####################################################
####################################################

##########################################################
##########################################################
InstallGlobalFunction(SparseBoundaryMatrix,
function(C,n)
local
        M,i;

M:=[1..C!.dimension(n)];
for i in [1..C!.dimension(n)] do
M[i]:=C!.boundary(n,i);
od;

return Objectify(HapSparseMat,
               rec( rows:=C!.dimension(n),
                    cols:=C!.dimension(n-1),
                    characteristic:=EvaluateProperty(C,"characteristic"),
                    mat:=M));

end);
##########################################################
##########################################################

#####################################################################
#####################################################################
InstallOtherMethod(RankMatDestructive,
"Rank of a sparse matrix",
[IsHapSparseMat],
function(M) 
SparseSemiEchelon(M);
return Length(Flat(M!.heads));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Rank,
"Rank of a sparse matrix",
[IsHapSparseMat],
function(M)
local mat, rk, tog, heads;
mat:=StructuralCopy(M!.mat);
if IsBound(M!.heads) then tog:=true; heads:=StructuralCopy(M!.heads); else tog:=false; fi;
SparseSemiEchelon(M);
rk:= Length(Flat(M!.heads));
if tog then M!.heads:=heads; else Unbind(M!.heads); fi;
M!.mat:=mat;
return rk;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Bettinumbers,
"Betti numbers of a sparse chain complex",
[IsHapSparseChainComplex,IsInt],
function(C,n)
local B,betti;

if n<0 then return 0; fi;

if n=0 then betti:=C!.dimension(n);
else
B:=SparseBoundaryMatrix(C,n);
betti:=C!.dimension(n)-RankMatDestructive(B);
fi;

B:=SparseBoundaryMatrix(C,n+1);
betti:=betti-RankMatDestructive(B);
return betti;

end);
#####################################################################
#####################################################################



