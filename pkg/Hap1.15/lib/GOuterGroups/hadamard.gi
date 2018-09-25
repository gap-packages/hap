
###############################################################
InstallGlobalFunction(CocyclicMatrices,
function(arg)
local G,q,a, A, R, C, H, h, L, n, c, M, g, row, g1, x ;
G:=arg[1];
#if Length(arg)=2 then q:=arg[2]; else q:=2; fi;
#a:=Group([[E(q)]]);;
a:=Group((1,2));
A:=TrivialGModuleAsGOuterGroup(G,a);;
R:=ResolutionFiniteGroup(G,3);;  
C:=HomToGModule(R,A);; 
H:=CohomologyModule(C,2);; 
h:=ActedGroup(H);; 

L:=[];
for n in [1..Length(Elements(h))] do
c:=H!.representativeCocycle(Elements(h)[n]);
 
M:=[];
row:=0;
for g in G do
row:=row+1;
M[row]:=[];
for g1 in G do
x:=c!.Mapping(g,g1);
if x=() then x:=1; else x:=-1; fi;
#x:=x[1][1];
Add(M[row],x);
od;od;

Add(L,M);
od;

return L;
end);
###############################################################

###############################################################
InstallGlobalFunction(IsHadamardMatrix,
function(H) local N,T;
if not IsMatrix(H) then return false; fi;
if not SSortedList(Flat(H))=[-1,1] then return false; fi;
N:=Length(H);
T:=1*TransposedMat(H);
#T:=List(T,r->List(r,x->ComplexConjugate(x)));
if not H*T=N*IdentityMat(N) then return false; fi;
return true;
end);
###############################################################

###############################################################
InstallGlobalFunction(HadamardGraph,
function(H)
local G, N, A, i, j;
N:=Length(H);
if not IsHadamardMatrix(H) then
Print("This function must be applied to a Hadamard matrix\n");
return fail;
fi;

A:=NullMat(4*N,4*N);
for i in [1..N] do
for j in [1..N] do
if H[i][j]=1 then
A[i][2*N+j]:=1;
A[N+i][3*N+j]:=1;
else
A[i][3*N+j]:=1;
A[N+i][2*N+j]:=1;
fi;
od;od;

return IncidenceMatrixToGraph(A);
end);
##############################################################

###############################################################
InstallGlobalFunction(CocyclicHadamardMatrices,
function(arg)
local G,q,L, C2, K, elts, M, x, g, h, gh, N, CHM;
G:=arg[1];
#if Length(arg)=2 then q:=arg[2]; else q:=2; fi;
elts:=Elements(G);
C2:=[-1,1];
#C2:=List([1..q],i->E(q)^i);
K:=Cartesian(List([1..Size(G)],i->C2));;
K:=Filtered(K,x->x[1]=1);
#L:=CocyclicMatrices(G,q);
L:=CocyclicMatrices(G);


CHM:=[];
for M in L do
N:=1*M;
for x in K do
for g in [1..Length(elts)] do
for h in [1..Length(elts)] do
gh:=Position(elts,elts[g]*elts[h]);
N[g][h]:=M[g][h]*x[g]*x[h]*x[gh]^-1;
od;
od;
if IsHadamardMatrix(N) then Add(CHM,1*N);
fi;
od;
od;

return SSortedList(CHM);
end);
###############################################################

