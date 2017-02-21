##########################################
InstallMethod(KendallMetric,
"for permutations g,h and the degree N of a group containing them",
[IsPerm,IsPerm,IsInt],
function(g,h,N)
local S, i,j;
S:=0;

for i in [1..N] do
for j in [i+1..N] do
if (i^g-j^g)*(i^h-j^h)<0 then S:=S+1; fi;
od;od;

return S;
end);
##########################################

##########################################
InstallMethod(KendallMetric,
"for permutations g,h",
[IsPerm,IsPerm],
function(g,h)
local N,M,S, i,j;

if IsOne(g) then N:=0;
else
N:=Maximum(ListPerm(g));
fi;
if IsOne(h) then M:=0;
else
M:=Maximum(ListPerm(h));
fi;

N:=Maximum(N,M);

S:=0;

for i in [1..N] do
for j in [i+1..N] do
if (i^g-j^g)*(i^h-j^h)<0 then S:=S+1; fi;
od;od;

return S;
end);
##########################################


##########################################
InstallMethod(CayleyMetric,
"for permutations g,h and the degree N of a group containing them",
[IsPerm,IsPerm,IsInt],
function(g,h,N)
local S,T,R,i,a;

a:=h*g^-1;
S:=Sum(Flat(CycleStructurePerm(a)));
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-S-T;

end);
##########################################

##########################################
InstallMethod(CayleyMetric,
"for permutations g,h ",
[IsPerm,IsPerm],
function(g,h)
local N,S,T,R,i,a;

N:=Maximum(ListPerm(g));
N:=Maximum(N,Maximum(ListPerm(h)));


a:=h*g^-1;
S:=Sum(Flat(CycleStructurePerm(a)));
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-S-T;

end);
##########################################


##########################################
InstallMethod(HammingMetric,
"for permutations g,h and the degree N of a group containing them",
[IsPerm,IsPerm,IsInt],

function(g,h,N)
local T,R,i,a;

a:=h*g^-1;
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-T;

end);
##########################################

##########################################
InstallMethod(HammingMetric,
"for permutations g,h",
[IsPerm,IsPerm],

function(g,h)
local N,T,R,i,a;

N:=Maximum(ListPerm(g));
N:=Maximum(N,Maximum(ListPerm(h)));

a:=h*g^-1;
T:=0;;
for i in [1..N] do
if i=i^a then T:=T+1; fi;
od;

return
N-T;

end);
##########################################

##########################################
InstallMethod(ManhattanMetric,
"for two vectors u and v",
[IsVector,IsVector],

function(u,v);

return Sum(List(u-v,x->AbsInt(x)));

end);
##########################################

##########################################
InstallMethod(EuclideanSquaredMetric,
"for two vectors u and v",
[IsVector,IsVector],

function(u,v);

return Sum(List(u-v,x->x^2));

end);
##########################################

##########################################
InstallMethod(EuclideanApproximatedMetric,
"for two vectors u and v",
[IsVector,IsVector],

function(u,v);

return Rat(Sqrt(Float(Sum(List(u-v,x->x^2)))));

end);
##########################################

####################################################
####################################################
InstallGlobalFunction(IsMetricMatrix,
function(A)
local x,y,z,N;

N:=Length(A);

#Is square matrix?
if not Length(A)=Length(A[1]) then return false; fi;

#d(x,x)=0?
for x in [1..N] do
if not A[x][x]=0 then return false;fi;
od;

#Symmetric? And d(x,y)>=0 and =0 only when x=y
for x in [1..N] do
for y in [x+1..N] do
if A[x][y]<=0 then return false; fi;
if not A[x][y]=A[y][x] then return false; fi;
od;
od;

#Triangle inequality?
for x in [1..N] do
for y in [x+1..N] do
for z in [y+1..N] do
if A[x][y]+A[y][z] <A[x][z] then return false; fi;
if A[x][z]+A[z][y] <A[x][y] then return false; fi;
if A[y][x]+A[x][z] <A[y][z] then return false; fi;
od;
od;
od;

return true;
end);
####################################################
####################################################

