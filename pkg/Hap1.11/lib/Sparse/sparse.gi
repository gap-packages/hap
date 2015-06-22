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

if IsZero(k) then M!.mat[i]:=[]; fi;

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
local x,z,r, pos;

#Adds k times the j-th row to the i-th row of a sparse matrix M.

Sort(M!.mat[i]); #Leave this in for safety. It should not be needed.

z:=List(M!.mat[i],a->a[1]);

for x in M!.mat[j] do
pos:=PositionSet(z,x[1]);
if  pos=fail then
Add(M!.mat[i],[x[1],k*x[2]]);
else
M!.mat[i][pos][2]:=M!.mat[i][pos][2]+k*x[2];
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

if not IsBound(M!.heads) then 

M!.heads:=[]; ;

# i = M!.heads[j] is the number of the row whose first non-zero entry
# is in column j.

V:=[1..M!.rows];
##########
fn:=function(i,j) 
local v; 
if IsBound(M!.mat[i]) and IsBound(M!.mat[j]) then 

if Length(M!.mat[i])=0 then return true; fi;
if Length(M!.mat[j])=0 then return false; fi;

return 
Maximum(List(M!.mat[i],x->x[1])) 
>
Maximum(List(M!.mat[j],x->x[1]));
fi;
return false;
end; 
#########
#Sort(V,fn);

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

for i in [1..M!.rows] do
if IsBound(M!.mat[i]) then
   M!.mat[i]:=Filtered(M!.mat[i],x->not x[2]=0);
   if Length(M!.mat[i])=0 then Unbind(M!.mat[i]); fi;
fi;
od;

fi;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(SparseRowReduce,
function(M,i)
local x, first, r,k, zz, xx, pos;

#Row reduce the i-th row of a sparse matrix if it is non-empty. 
#The first non-zero entry of the reduced row will be 1.
#Return true if the matrix is modified, and false otherwise.

if not IsBound(M!.mat[i]) then return false; fi;
if Length(M!.mat[i])=0 then return false; fi;

first:=M!.mat[i][1][1];

   while IsBound(M!.heads[first]) do
   r:=M!.heads[first]; 
   k:=-M!.mat[i][1][2];
###SparseRowAdd(M,i,r,k);
#############################################
zz:=List(M!.mat[i],a->a[1]);
for xx in M!.mat[r] do
pos:=PositionSet(zz,xx[1]);
if  pos=fail then
Add(M!.mat[i],[xx[1],k*xx[2]]);
else
M!.mat[i][pos][2]:=M!.mat[i][pos][2]+k*xx[2];
fi;
od;
M!.mat[i]:=Filtered(M!.mat[i],a->not IsZero(a[2]));
Sort(M!.mat[i]);
#############################################
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

##########################################################
##########################################################
InstallGlobalFunction(TransposeOfSparseMat,
function(A)
local
        M, entry, r;

M:=List([1..A!.cols],i->[]);
for r in [1..A!.rows] do
for entry in A!.mat[r] do
Add(M[entry[1]],[r,entry[2]]);
od;od;

return Objectify(HapSparseMat,
               rec( rows:=A!.cols,
                    cols:=A!.rows,
                    characteristic:=EvaluateProperty(A,"characteristic"),
                    mat:=M));

end);
##########################################################
##########################################################

#####################################################################
#####################################################################
InstallOtherMethod(RankMatDestructive,
"Rank of a sparse matrix",
[IsHapSparseMat],
function(M) local T;
if M!.cols>=M!.rows then 
SparseSemiEchelon(M);
return Length(Flat(M!.heads));
else
T:=TransposeOfSparseMat(M);
SparseSemiEchelon(T);
return Length(Flat(T!.heads));
fi;
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

#####################################################################
#####################################################################
InstallGlobalFunction(NullspaceSparseMatDestructive,
function(B)
local NS, heads, free, vec, revcol,  pos, i,v,w, x,y,z,h,s;

NS:=[];

if not IsBound(B!.heads) then
SparseSemiEchelon(B);
fi;

heads:=Filtered([1..B!.cols],x->IsBound(B!.heads[x]));

free:=Difference([1..B!.cols],heads);
revcol:=Reversed([1..B!.cols]);

for x in free do
v:=[];
Add(v,[x,1]);

for y in revcol do
if y in heads then
h:=B!.mat[B!.heads[y]];

s:=0;
for z in h do
pos:=PositionProperty(v,a->a[1]=z[1]);
if IsInt(pos) then s:=s+z[2]*v[pos][2]; fi;
od;
if not s=0 then Add(v,[y,-s]); fi;
fi;
od;

Add(NS,v);

od;
return NS;

end);
#####################################################################
#####################################################################

###############################################################
###############################################################
InstallGlobalFunction(ReverseSparseMat,
function(B)
local i,x,d;
d:=B!.cols+1;

for i in [1..B!.rows] do
if IsBound(B!.mat[i]) then
for x in B!.mat[i] do
x[1]:=d-x[1];
od;
fi;
od;

if IsBound(B!.heads) then Unbind(B!.heads); fi;

end);
###############################################################
###############################################################




