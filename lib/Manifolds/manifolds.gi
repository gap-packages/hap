#####################################################
#####################################################
InstallGlobalFunction(HAP_WedgeSumOfSimplicialComplexes,
function(K,L)
local W;

W:=MaximalSimplicesOfSimplicialComplex(K);
Append(W,  (K!.nrSimplices(0)-1)+MaximalSimplicesOfSimplicialComplex(L));

return MaximalSimplicesToSimplicialComplex(W);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallMethod(WedgeSum,
"wedge of simplicial complexes",
[IsHapSimplicialComplex,IsHapSimplicialComplex],
function(K,L)
return HAP_WedgeSumOfSimplicialComplexes(K,L);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(WedgeSum,
"wedge of simplicial complexes",
[IsHapSimplicialComplex,IsHapSimplicialComplex,IsHapSimplicialComplex],
function(K,L,M)
return HAP_WedgeSumOfSimplicialComplexes(HAP_WedgeSumOfSimplicialComplexes(K,L),M);
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallOtherMethod(WedgeSum,
"wedge of regulat CW complexes",
[IsHapRegularCWComplex,IsHapRegularCWComplex],
function(K,L)
return RegularCWComplex_WedgeSum(K,L,1,1);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(WedgeSum,
"wedge of regulat CW complexes",
[IsHapRegularCWComplex,IsHapRegularCWComplex,IsHapRegularCWComplex],
function(K,L,M)
return RegularCWComplex_WedgeSum(RegularCWComplex_WedgeSum(K,L,1,1),M,1,1);
end);
#####################################################
#####################################################



#####################################################
#####################################################
InstallOtherMethod(ClosedSurface,
"simplicial surface or genus +/- g (where - means nonorientable)",
[IsInt],
function(n)
local M, i;
if n=0 then return Sphere(n); fi;

if n>0 then M:=HAP_SimplicialTorus();;
for i in [2..n] do
M:=ConnectedSum(M,HAP_SimplicialTorus());
od;
return M;
fi;

M:=HAP_SimplicialProjectivePlane();;
for i in [2..-n] do
M:=ConnectedSum(M,HAP_SimplicialProjectivePlane());
od;
return M;

end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(ClosedSurface,
"simplicial surface or genus +/- g (where - means nonorientable)",
[IsInt,IsString],
function(n,str)
local Y;
if str="CW" then
Y:=ClosedSurface(n);
Y:=RegularCWComplex(Y);;
Y:=BarycentricallySimplifiedComplex(Y);
return Y;
fi;
if str="Simplicial" then
Y:=ClosedSurface(n);
return Y;
fi;
Print("The second argument should be \"CW\" or \"Simplicial\".\n");
return fail;
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallGlobalFunction(HAP_SimplicialTorus, 
function()
local S;

#S:=[ [1,2,4], [2,4,6], [2,3,6], [3,6,7], [1,3,7], [1,4,7],
#     [4,5,6], [5,6,8], [6,7,8], [7,8,9], [4,7,9], [4,5,9],
#     [1,5,8], [1,2,8], [2,8,9], [2,3,9], [3,5,9], [1,3,5]];; 
S:=[[1,2,6],[2,6,7],[2,3,7],[1,3,7],[4,6,7],[4,5,7],[1,5,7],
    [1,5,6],[1,2,4],[2,4,5],[2,3,5],[3,5,6],[3,4,6],[1,3,4]];;
return MaximalSimplicesToSimplicialComplex(S);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(HAP_SimplicialProjectivePlane,
function()
local S;

S:= [[1,5,8],[1,2,8],[2,7,8],[2,3,7],[3,4,7],
     [4,5,8],[4,8,9],[7,8,9],[6,7,9],[4,6,7],[4,5,6],
     [3,4,9],[2,3,9],[2,6,9],[1,2,6],[1,5,6]];;
return MaximalSimplicesToSimplicialComplex(S);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(SimplicialK3Surface,
function()
local S;

S:=[[2, 10, 13, 15, 16], [3, 5, 10, 13, 15], [2, 5, 7, 8, 10], [1, 9, 11, 13, 14], [1, 2, 8, 10, 12], [1, 3, 5, 6, 11], [1, 5, 6, 9, 12], [1, 4, 5, 6, 14], [1, 4, 10, 13, 14], [1, 9, 10, 14, 15], [2, 4, 7, 8, 12], [3, 4, 6, 10, 12], [1, 6, 7, 8, 9], [3, 4, 5, 7, 15], [1, 7, 12, 15, 16], [4, 5, 7, 13, 16], [5, 8, 11, 12, 15], [2, 4, 7, 12, 14], [1, 4, 5, 14, 16], [2, 5, 6, 10, 11], [1, 6, 8, 12, 14], [5, 8, 9, 14, 16], [5, 10, 11, 12, 13], [2, 4, 8, 9, 12], [7, 9, 12, 15, 16], [1, 2, 6, 9, 15], [1, 5, 14, 15, 16], [2, 3, 4, 5, 9], [6, 8, 10, 11, 15], [1, 5, 8, 10, 12], [1, 3, 7, 9, 10], [6, 7, 8, 9, 13], [1, 2, 9, 11, 15], [2, 8, 11, 14, 16], [2, 4, 5, 13, 16], [1, 4, 8, 13, 15], [4, 7, 8, 10, 11], [2, 4, 8, 9, 10], [2, 3, 4, 9, 13], [2, 8, 10, 12, 13], [1, 2, 4, 11, 15], [2, 3, 9, 11, 15], [2, 8, 11, 15, 16], [3, 4, 5, 9, 11], [6, 10, 13, 15, 16], [8, 10, 11, 15, 16], [6, 7, 11, 13, 15], [1, 5, 7, 15, 16], [4, 5, 7, 9, 15], [3, 4, 6, 7, 16], [2, 3, 11, 14, 16], [3, 4, 9, 11, 13], [1, 2, 5, 14, 15], [2, 3, 9, 13, 14], [1, 7, 12, 13, 16], [1, 2, 5, 13, 16], [2, 3, 7, 8, 12], [2, 9, 11, 12, 14], [1, 9, 11, 15, 16], [2, 4, 7, 8, 10], [1, 4, 9, 13, 14], [1, 2, 3, 12, 16], [8, 11, 12, 14, 15], [1, 2, 6, 13, 16], [1, 4, 10, 12, 13], [3, 6, 8, 10, 16], [2, 7, 8, 14, 16], [1, 6, 8, 9, 12], [6, 9, 10, 14, 16], [5, 8, 11, 12, 16], [5, 9, 10, 14, 15], [3, 9, 12, 15, 16], [4, 6, 10, 11, 15], [2, 4, 9, 10, 16], [5, 8, 9, 13, 15], [2, 3, 6, 9, 15], [6, 11, 12, 14, 16], [2, 3, 10, 13, 15], [3, 5, 9, 11, 16], [3, 4, 8, 11, 13], [3, 4, 5, 7, 13], [5, 7, 8, 10, 14], [4, 12, 13, 14, 15], [6, 7, 10, 14, 16], [5, 10, 11, 13, 14], [3, 4, 7, 13, 16], [6, 8, 9, 12, 13], [1, 3, 4, 10, 14], [2, 4, 6, 11, 12], [1, 7, 9, 10, 14], [4, 6, 8, 13, 14], [4, 9, 10, 11, 16], [3, 7, 8, 10, 16], [5, 7, 9, 15, 16], [1, 7, 9, 11, 14], [6, 8, 10, 15, 16], [4, 6, 9, 14, 16], [5, 8, 9, 10, 14], [7, 8, 10, 14, 16], [2, 6, 7, 9, 11], [7, 9, 10, 13, 15], [3, 6, 7, 10, 12], [2, 4, 6, 10, 11], [4, 5, 8, 9, 11], [1, 2, 3, 8, 16], [1, 3, 8, 9, 12], [1, 2, 6, 8, 14], [3, 5, 6, 13, 15], [1, 5, 6, 12, 14], [2, 5, 7, 14, 15], [1, 5, 10, 11, 12], [3, 7, 8, 10, 11], [1, 2, 6, 14, 15], [1, 2, 6, 8, 16], [7, 9, 10, 12, 15], [3, 4, 6, 8, 14], [3, 7, 13, 14, 16], [2, 5, 7, 8, 14], [6, 7, 9, 10, 14], [2, 3, 7, 12, 14], [4, 10, 12, 13, 14], [2, 5, 6, 11, 13], [4, 5, 6, 7, 16], [1, 3, 12, 13, 16], [1, 4, 11, 15, 16], [1, 3, 4, 6, 10], [1, 10, 11, 12, 13], [3, 9, 11, 15, 16], [3, 5, 10, 14, 15], [5, 8, 9, 10, 13], [1, 2, 5, 7, 15], [2, 4, 11, 12, 14], [3, 11, 13, 14, 16], [1, 2, 5, 7, 13], [4, 7, 8, 9, 15], [1, 5, 6, 10, 11], [6, 7, 10, 13, 15], [3, 4, 7, 14, 15], [3, 4, 10, 12, 14], [1, 2, 6, 7, 13], [2, 3, 4, 5, 13], [5, 8, 12, 13, 15], [4, 6, 9, 13, 14], [2, 4, 5, 6, 12], [2, 9, 10, 13, 16], [8, 11, 12, 14, 16], [1, 7, 12, 13, 15], [8, 12, 13, 14, 15], [2, 8, 9, 12, 13], [4, 6, 10, 12, 15], [2, 8, 11, 14, 15], [2, 6, 9, 11, 12], [8, 9, 10, 11, 16], [2, 3, 6, 13, 15], [2, 3, 12, 15, 16], [1, 3, 5, 9, 12], [2, 5, 6, 9, 12], [5, 6, 8, 11, 15], [2, 6, 13, 15, 16], [2, 3, 11, 15, 16], [3, 5, 6, 8, 15], [2, 4, 5, 9, 12], [2, 10, 12, 13, 14], [6, 8, 12, 13, 14], [1, 2, 3, 8, 12], [1, 4, 7, 8, 11], [5, 6, 7, 12, 16], [3, 5, 7, 13, 14], [3, 4, 5, 8, 11], [6, 7, 11, 12, 15], [3, 4, 6, 7, 12], [1, 2, 4, 7, 11], [3, 9, 10, 12, 15], [4, 10, 12, 15, 16], [3, 5, 7, 14, 15], [3, 9, 11, 13, 14], [5, 9, 14, 15, 16], [4, 5, 6, 7, 12], [1, 3, 6, 10, 11], [1, 3, 9, 10, 15], [4, 7, 8, 9, 12], [5, 9, 10, 13, 15], [1, 3, 8, 13, 16], [2, 9, 12, 13, 14], [6, 7, 10, 12, 15], [2, 6, 8, 14, 15], [3, 5, 6, 8, 11], [3, 4, 7, 12, 14], [1, 3, 10, 14, 15], [7, 11, 12, 13, 16], [3, 11, 12, 13, 16], [3, 4, 5, 8, 15], [7, 11, 13, 14, 16], [2, 4, 7, 14, 15], [1, 2, 10, 12, 16], [1, 6, 8, 13, 16], [1, 7, 8, 13, 15], [6, 9, 11, 12, 14], [3, 6, 8, 14, 15], [2, 4, 11, 14, 15], [3, 7, 9, 10, 12], [1, 3, 6, 14, 15], [2, 4, 5, 6, 10], [1, 4, 9, 14, 16], [1, 3, 7, 8, 9], [5, 7, 9, 12, 16], [1, 3, 7, 10, 11], [7, 8, 9, 13, 15], [1, 4, 7, 8, 15], [1, 4, 10, 12, 16], [1, 7, 10, 11, 14], [1, 2, 6, 7, 9], [1, 3, 11, 12, 13], [1, 5, 7, 13, 16], [5, 7, 10, 11, 14], [2, 10, 12, 15, 16], [3, 6, 7, 10, 16], [1, 2, 5, 8, 10], [4, 10, 11, 15, 16], [5, 8, 10, 12, 13], [3, 6, 8, 10, 11], [4, 5, 7, 9, 12], [6, 7, 11, 12, 16], [2, 8, 9, 10, 13], [8, 9, 10, 14, 16], [3, 4, 6, 8, 16], [1, 10, 11, 13, 14], [1, 2, 5, 8, 14], [2, 4, 5, 10, 16], [1, 2, 7, 9, 11], [1, 3, 5, 6, 9], [5, 7, 11, 13, 14], [3, 5, 10, 13, 14], [2, 3, 9, 11, 14], [4, 11, 12, 14, 15], [2, 3, 7, 14, 16], [3, 4, 8, 13, 16], [6, 7, 9, 11, 14], [5, 6, 11, 13, 15], [4, 5, 6, 14, 16], [3, 4, 8, 14, 15], [4, 5, 8, 9, 15], [1, 4, 8, 11, 13], [5, 6, 12, 14, 16], [2, 3, 10, 12, 14], [1, 2, 5, 10, 16], [2, 5, 7, 10, 11], [2, 6, 7, 11, 13], [1, 4, 5, 10, 16], [2, 6, 8, 15, 16], [2, 3, 10, 12, 15], [7, 11, 12, 13, 15], [1, 3, 8, 11, 13], [4, 8, 9, 10, 11], [1, 9, 14, 15, 16], [1, 3, 6, 9, 15], [6, 9, 12, 13, 14], [2, 3, 10, 13, 14], [2, 5, 7, 11, 13], [2, 3, 5, 6, 13], [4, 6, 8, 13, 16], [6, 7, 9, 10, 13], [5, 8, 12, 14, 16], [4, 6, 9, 13, 16], [5, 8, 9, 11, 16], [2, 3, 5, 6, 9], [1, 3, 5, 11, 12], [3, 7, 8, 9, 12], [4, 6, 11, 12, 15], [3, 5, 9, 12, 16], [5, 11, 12, 13, 15], [1, 3, 4, 6, 14], [3, 5, 11, 12, 16], [1, 5, 8, 12, 14], [4, 8, 13, 14, 15], [1, 3, 7, 8, 11], [6, 9, 10, 13, 16], [2, 4, 9, 13, 16], [1, 6, 7, 8, 13], [1, 4, 12, 13, 15], [2, 4, 7, 10, 11], [1, 4, 9, 11, 13], [6, 7, 11, 14, 16], [1, 4, 9, 11, 16], [1, 4, 12, 15, 16], [1, 2, 4, 7, 15], [2, 3, 7, 8, 16], [1, 4, 5, 6, 10]];;
return MaximalSimplicesToSimplicialComplex(S);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(Sphere,
"n-sphere as a simplicial complex",
[IsInt],
function(n)
local boundaries, k;

return MaximalSimplicesToSimplicialComplex( Combinations([1..n+2],n+1) );
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(ComplexProjectiveSpace,
"Complex projective plane as simplicial complex",
[IsInt],
function(n)
local B;

if n>2 or n<1 then 
Print("This is currently implemented only for n=1, 2\n");
return fail;
fi;

if n=1 then return Sphere(2); fi;

B:=
[ [ 1, 2, 4, 5, 6 ], [ 1, 2, 4, 5, 9 ], [ 1, 2, 5, 6, 8 ], 
  [ 1, 2, 6, 4, 7 ], [ 2, 3, 4, 5, 8 ], [ 2, 3, 5, 6, 4 ], 
  [ 2, 3, 5, 6, 7 ], [ 2, 3, 6, 4, 9 ], [ 3, 1, 4, 5, 7 ],
  [ 3, 1, 5, 6, 9 ], [ 3, 1, 6, 4, 5 ], [ 3, 1, 6, 4, 8 ], 
  [ 4, 5, 7, 8, 3 ], [ 4, 5, 7, 8, 9 ], [ 4, 5, 8, 9, 2 ], 
  [ 4, 5, 9, 7, 1 ], [ 5, 6, 7, 8, 2 ], [ 5, 6, 8, 9, 1 ],
  [ 5, 6, 8, 9, 7 ], [ 5, 6, 9, 7, 3 ], [ 6, 4, 7, 8, 1 ], 
  [ 6, 4, 8, 9, 3 ], [ 6, 4, 9, 7, 2 ], [ 6, 4, 9, 7, 8 ], 
  [ 7, 8, 1, 2, 3 ], [ 7, 8, 1, 2, 6 ], [ 7, 8, 2, 3, 5 ],
  [ 7, 8, 3, 1, 4 ], [ 8, 9, 1, 2, 5 ], [ 8, 9, 2, 3, 1 ], 
  [ 8, 9, 2, 3, 4 ], [ 8, 9, 3, 1, 6 ], [ 9, 7, 1, 2, 4 ], 
  [ 9, 7, 2, 3, 6 ], [ 9, 7, 3, 1, 2 ], [ 9, 7, 3, 1, 5 ] ];;

return MaximalSimplicesToSimplicialComplex(B);
end);
#####################################################
#####################################################

#####################################################
#####################################################
## Inputs simplicial complexes K, L representing 
## n-manifolds and integer e=+ or e=-1. It returns the 
## connected sum K#L or K#-L depending on e. 
## The implementation assumes K!.vertices=[1..m] and
## L!.vertices=[1..n]
##
SimplicialComplexConnectedSum:=function(arg)
local K, L, e, V, W, maxK, maxL,newmaxL, x,y,w,N,reindex, pos,a, b;

K:=arg[1];
L:=arg[2];

if not Dimension(K)=Dimension(L) then
Print("The simplicial manifolds must have the same dimension.\n");
return fail;
fi;

e:=-arg[3]; 
maxK:=1*K!.simplicesLst[Dimension(K)+1];
maxL:=1*L!.simplicesLst[Dimension(L)+1];

if (not Size(MaximalSimplicesToSimplicialComplex(maxK)) = Size(K))
or (not Size(MaximalSimplicesToSimplicialComplex(maxL)) = Size(L))
then
Print("The simplicial complexes are not pure.\n");
return fail;
fi;

V:=maxK[1];
W:=maxL[1];
maxK:=maxK{[2..Length(maxK)]};
maxK:=Concatenation(maxK,Combinations(V,Length(V)-1));
maxL:=maxL{[2..Length(maxL)]};
maxL:=Concatenation(maxL,Combinations(W,Length(W)-1));
N:=K!.nrSimplices(0)+1;
reindex:=[];
for w in L!.vertices do 
pos:=Position(W,w);
if not pos=fail then reindex[w]:=V[pos];
else reindex[w]:=N; N:=N+1; fi;
od;

if e<0 then
W:=Flat(W);
a:=reindex[W[1]];
b:=reindex[W[2]];
reindex[W[1]]:=b;
reindex[W[2]]:=a;
fi;

newmaxL:=[];
for x in maxL do
Add(newmaxL, List(x,i->reindex[i]) );
od;

return MaximalSimplicesToSimplicialComplex( Concatenation(maxK,newmaxL) );;
end;
#####################################################
#####################################################

#####################################################
#####################################################
InstallMethod(ConnectedSum,
"Connected sum of simplicial manifolds",
[IsHapSimplicialComplex,IsHapSimplicialComplex,IsInt],
function(K,L,e);
return SimplicialComplexConnectedSum(K,L,e);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(ConnectedSum,
"Connected sum of simplicial manifolds",
[IsHapSimplicialComplex,IsHapSimplicialComplex],
function(K,L);
return SimplicialComplexConnectedSum(K,L,1);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(ConnectedSum,
"Connected sum of regular CW-manifolds  [temporary method -- needs replacing]",
[IsHapRegularCWComplex,IsHapRegularCWComplex,IsInt],
function(KK,LL,e)
local K, L;
K:=BarycentricSubdivision(KK);
L:=BarycentricSubdivision(LL);
return SimplifiedComplex(RegularCWComplex(ConnectedSum(K,L,e)));
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(ConnectedSum,
"Connected sum of regular CW-manifolds  [temporary method -- needs replacing]",
[IsHapRegularCWComplex,IsHapRegularCWComplex],
function(KK,LL)
local K, L;
K:=BarycentricSubdivision(KK);
L:=BarycentricSubdivision(LL);
return SimplifiedComplex(RegularCWComplex(ConnectedSum(K,L)));
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(Cohomology,
"Integral cohomology of a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(K,n);
return Cohomology(RegularCWComplex(K),n);
end);
#####################################################
#####################################################

#####################################################
#####################################################
IntersectionForm:=function(arg)
local Y, H, gens, A, i, j, cup, c, x, y, n,B;

Y:=arg[1];
#if IsHapSimplicialComplex(YY) then
#Y:=RegularCWComplex(YY); 
#else
#Y:=YY;
#fi;

if not ( IsEvenInt(Dimension(Y)) and  Homology(Y,Dimension(Y))=[0]   )
then Print("The argument must have even dimension and top homology equal to Z.\n"); return fail; fi;

n:=Dimension(Y)/2;

H:=Cohomology(Y,n);
gens:=Filtered([1..Length(H)],i->H[i]=0);
A:=NullMat(Length(gens),Length(gens));

cup:=CupProduct(Y);
if Length(arg)=1 then
B:=IdentityMat(Length(H));
B:=B{gens};
else B:=arg[2];
fi;

for i in [1..Length(B)] do
for j in [i..Length(B)] do
c:=cup(n,n,B[i],B[j]);
A[i][j]:=c[1];
A[j][i]:=(-1)^n*c[1];
od;
od;
return A;
end;
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(CohomologyRingOfSimplicialComplex,
function(K,prime)
local R, dims, Y, C, D, n, k, i, j, ii, jj, dim, SCT, Degree, 
      cup, V, pair2int, int2pair, cnt, x, y, c, s, dims1;

Y:=RegularCWComplex(K);
C:=ChainComplex(Y);
D:=HomToIntegersModP(C,prime);
dims:=List([0..Dimension(K)], i->Cohomology(D,i));
dim:=Sum(dims);
cnt:=0;
dims1:=[];
dims1[1]:=0;
for k in [1..Length(dims)] do
cnt:=cnt+dims[k];
dims1[k+1]:=cnt;
od;
cup:=CupProduct(K,prime);

###############################
cnt:=0;
pair2int:=[];
int2pair:=[];
for n in [1..Length(dims)] do
pair2int[n]:=[];
for k in [1..dims[n]] do
cnt:=cnt+1;
pair2int[n][k]:=cnt;
int2pair[cnt]:=[n,k];
od;
od;
###############################

SCT:=EmptySCTable(dim,Zero(GF(prime)));
cnt:=0;

for i in [1..dim] do
for j in [i..dim] do
ii:=int2pair[i];
jj:=int2pair[j];
x:=[1..dims[ii[1]]]*0; x[ii[2]]:=1;
y:=[1..dims[jj[1]]]*0; y[jj[2]]:=1;
if ii[1]-1+jj[1]-1<=Length(dims)-1 then
c:=cup(ii[1]-1,jj[1]-1, x, y);
s:=dims1[ii[1]-1+jj[1]];
V:=[];
for k in [1..Length(c)] do
if not c[k]=0 then
Add(V, c[k]);
Add(V, s+k);
fi;
od;
SetEntrySCTable(SCT,i,j,V);
fi;
od;
od;

for i in [1..dim] do
for j in [1..i-1] do
ii:=int2pair[i];
jj:=int2pair[j];
x:=[];
x[1]:=List(SCT[j][i][1],a->1*a);
x[2]:=List(SCT[j][i][2],a->1*a);
if IsOddInt((ii[1]-1)*(jj[1]-1)) then
x[2]:=-1*x[2];
fi;
SCT[i][j]:=x;
od;od;

R:=AlgebraByStructureConstants(GF(prime),SCT);
R!.chainComplex:=C;

#####################################################################
Degree:=function(x)
local i, w, bas;
# returns the highest degree of a non-zero coefficient of x

if IsZero(x) then return 0; fi;
i:=Position(GeneratorsOfAlgebra(R),x);

if i=1 then return 0; fi;

if not i=fail then return int2pair[i][1]-1; fi;

bas:=Basis(R);
w:=Coefficients(bas,x);
w:=Filtered([2..Length(bas)],i->not IsZero(w[i]));
w:=List(w,i->int2pair[i][1]-1);
return Maximum(w);
end;
#####################################################################

R!.degree:=Degree;
R!.int2pair:=int2pair;
R!.pair2int:=pair2int;
#R!.bockstein:=HAP_bockstein(R); 
return R;
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(CohomologyRing,
"cohomology of a simplicial complex over a field of p elements",
[IsHapSimplicialComplex,IsInt],
function(K,prime);
return CohomologyRingOfSimplicialComplex(K,prime);
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallOtherMethod(CohomologyRing,
"cohomology of a regular CW-complex complex over a field of p elements",
[IsHapRegularCWComplex,IsInt],
function(W,prime)
local K;
#K:=BarycentricSubdivision(W);     #WE SHOULD USE THE DIAGONAL FOR W HERE!!!
                                  #THIS IS A TEMPORARY METHOD
return CohomologyRingOfSimplicialComplex(W,prime);
end);
#####################################################
#####################################################

