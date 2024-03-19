#  APPENDIX TO: CELLULAR APPROXIMATIONS TO THE DIAGONAL MAP
#  by GRAHAM ELLIS

LoadPackage("HAP");
##To time a piece of code do
#gap> t:=Runtime();
#gap>  ---the code
#gap> Print(Runtime()-t);

##To obtain an upper bound on the memory needed for a piece of code do
#gap> t:=TotalMemoryAllocated();
#gap> ---the code
#gap> Print(TotalMemoryAllocated()-t);

###########################################################
###########################################################
#Algorithm 1 applied to 18-dimensional simplex (Table 1)
Y:=RegularCWSimplex(18,false);;
Print("Number of cells = ", Size(Y),"\n");
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..100] do
Y:=RegularCWComplexReordered(Y);
c:=CriticalCells(Y);;
Add(C,Length(c));
od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Algorithm 1 applied to 12-dimensional cube (Table 1)
Y:=RegularCWCube(12);;
Print("Number of cells = ", Size(Y),"\n");
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..100] do
Y:=RegularCWComplexReordered(Y);
c:=CriticalCells(Y);;
Add(C,Length(c));
od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Algorithm 1 applied to 9-dimensional associahedron (Table 1)
Y:=RegularCWAssociahedron(11,false);; 
Print("Number of cells = ", Size(Y),"\n");
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..100] do
Y:=RegularCWComplexReordered(Y);
c:=CriticalCells(Y);;
Add(C,Length(c));
od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Algorithm 1 applied to 7-dimensional permutahedron (Table 1)
Y:=RegularCWPermutahedron(7);;
Print("Number of cells = ", Size(Y),"\n");
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..100] do
Y:=RegularCWComplexReordered(Y);
c:=CriticalCells(Y);;
Add(C,Length(c));
od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Algorithm 1 applied to random 7-dimensional polytope (Table 1)
P:=[];;
for i in [1..400] do
Add(P , List([1..7],i->Random([0..20])/Random([1..20])) );
od;
Y:=RegularCWPolytope(P);;
Print("Number of cells = ", Size(Y),"\n");
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..100] do
Y:=RegularCWComplexReordered(Y);
c:=CriticalCells(Y);;
Add(C,Length(c));
od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Algorithm 1 applied to the CW-complex [0,2]x[0,3]x[0,5]
A:=List([1..2],k->List([1..3],j->List([1..5],i->1)));;
Y:=PureCubicalComplex(A);;
Y:=RegularCWComplex(Y);;
c:=CriticalCells(Y);;
C:=[Length(c)];;
for i in [1..10^6] do
 Y:=RegularCWComplexReordered(Y);
 c:=CriticalCells(Y);;
 Add(C,Length(c));
 od;
Print("Maximum number critical cells = ", Maximum(C), "\n");

###########################################################
###########################################################
#Computing "size" of the diagonal map using Algorithms 1 & 2 (Table 2)
for n in [1..5] do
K:=RegularCWAssociahedron(n+2);;
Unbind(K!.directed); #To ensur Algs 1 and 2 are used
D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWPermutahedron(n);;
Unbind(K!.directed); #To ensur Algs 1 and 2 are used
D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWCube(n);;
Unbind(K!.directed); #To ensur Algs 1 and 2 are used
D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWSimplex(n);;
Unbind(K!.directed); #To ensur Algs 1 and 2 are used
D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;

###########################################################
###########################################################
#Enumerating the planar rooted trees with 5 leaves in Exa,ple 3.9
#The first tree is displayed.
K:=RegularCWAssociahedron(5);        
trees:=K!.trees;
t:=trees[1][1];  #first tree       
HAP_DisplayPlanarTree(t);

###########################################################
###########################################################
#Computing "size" of the diagonal map using "classical formulae" (Table 2)
for n in [1..5] do
K:=RegularCWAssociahedron(n+2);;D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWPermutahedron(n);;D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWCube(n);;D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;
for n in [1..5] do
K:=RegularCWSimplex(n);;D:=DiagonalChainMap(K);;
s:=Size(Filtered(D!.mapping([1],n), x -> not IsZero(x)));
Print(s,"\n");
od;

###########################################################
###########################################################
#Commands for Example 7.1
K:=SimplicialK3Surface();;
KpK:=ConnectedSum(K,K,+1);;
KmK:=ConnectedSum(K,K,-1);;
SignatureOfSymmetricMatrix(IntersectionForm(KpK));
SignatureOfSymmetricMatrix(IntersectionForm(KmK));

###########################################################
###########################################################
#Commands for Example 7.2
ap:=[[6,10],[5,9],[10,7],[1,6],[8,5],[9,4],[7,3],[4,2],[3,1],[2,8]];;
Y:=SphericalKnotComplement(ap);;
W:=ContractedComplex(Y);;
List([0..2],n->Cohomology(W,n));
Display(IntersectionForm(W));

###########################################################
###########################################################
#Commands for Example 7.3
Y:=ContractedComplex(SpunLinkComplement([[1,3],[2,4],[1,3],[2,4]]));;
Size(Y);
W:=SimplifiedComplex(Y);;
List([0..3],n->Cohomology(W,n));
Display(CupProductMatrix(W,1,1));;
Display(CupProductMatrix(W,1,2));;

###########################################################
###########################################################
#Commands for Example 7.4
T:=ClosedSurface(2,"CW");;
W:=DirectProduct(T,T);;
A:=IntersectionForm(W);;time;

K:=BarycentricSubdivision(W);;
B:=IntersectionForm(K);;time;

T:=ClosedSurface(2,"Simplicial");;
W:=DirectProduct(T,T);;
A:=IntersectionForm(W);;time;

###########################################################
###########################################################
#Commands for Example 7.5
fl:=HapFile("data247.txt");;Read(fl);;
F:=ThickeningFiltration(T,19);;
P:=PersistentBettiNumbersAlt(F,[0,1,2]);;time; 
BarCodeCompactDisplay(P); 
X12:=RegularCWComplex(FiltrationTerm(F,12));;
cup:=LowDimensionalCupProduct(X12);;
cup([1,0],[0,1]);

##########################################################
###########################################################
#Construction of Poincare's third manifold (page 31)
P1:=[[1,4,3,2],[8,7,6,5]];;
P2:=[[1,4,8,5],[6,2,3,7]];;
P3:=[[1,2,6,5],[3,7,8,4]];;
M:=PoincareCubeCWComplex(P1,P2,P3);
IsClosedManifold(M);
StructureDescription(FundamentalGroup(M));

##########################################################
##########################################################
#Construction of a manifold using Dehn surgery (page33)
ap:=ArcPresentation(PureCubicalKnot(3,1));;
W:=ThreeManifoldViaDehnSurgery(ap,17,16);
Size(W);

##########################################################
##########################################################
#Computation of the Dijkgraaf-Witten invariant on page 34.
ap:=[[2,1],[2,1]];; #Arc presentation for the trivial knot
L51:=ThreeManifoldViaDehnSurgery(ap,5,1);;
D:=DijkgraafWittenInvariant(L51,CyclicGroup(5));
L52:=ThreeManifoldViaDehnSurgery(ap,5,2);;
D:=DijkgraafWittenInvariant(L52,CyclicGroup(5));

##########################################################
##########################################################
#Computation of the linking form invariants on page 35.
LensSpaces:=[];;
for q in [1..12] do
Add(LensSpaces,ThreeManifoldViaDehnSurgery([[1,2],[1,2]],13,q));
od;
c:=Classify([1..12],q->LinkingFormHomotopyInvariant(LensSpaces[q]));

##########################################################
##########################################################
#Computation of the cohomology ring of a lens space (page 35).
ap:=[[2,1],[2,1]];; #Arc presentation for the trivial knot
L51:=ThreeManifoldViaDehnSurgery(ap,5,1);;
v:=CohomologyRing(L51,5);
v.2*v.3;

##########################################################
##########################################################
#Computation of the cohomology ring of a covering space (page 35)
ap:=ArcPresentation(PureCubicalKnot(3,1));;
W:=ThreeManifoldViaDehnSurgery(ap,5,2);;
U:=UniversalCover(W);;F:=U!.group;;
H:=LowIndexSubgroupsFpGroup(F,6)[7];;
Index(F,H);
V:=EquivariantCWComplexToRegularCWComplex(U,H);;
V:=BarycentricallySimplifiedComplex(V);;
R:=CohomologyRing(V,7);

##########################################################
##########################################################
#Computation of the cohomology ring of a 6-fold cover of 
#Seifert-Weber space (page 36) 
W:=PoincareDodecahedronCWComplex(
 [[1,2,3,4,5],[7,8,9,10,6]],
 [[1,11,16,12,2],[9,8,18,14,19]],
 [[2,12,17,13,3],[10,9,19,15,20]],
 [[3,13,18,14,4],[6,10,20,11,16]],
 [[4,14,19,15,5],[7,6,16,12,17]],
 [[5,15,20,11,1],[8,7,17,13,18]]);
U:=UniversalCover(W);
G:=U!.group;;
L:=LowIndexSubgroupsFpGroup(G,6);;
H:=L[86];;
WH:=EquivariantCWComplexToRegularCWComplex(U,H);
v:=CohomologyRing(WH,2);
v.4^3;

###########################################################
###########################################################
#Construction of the iterated suspension of the real projective plane
#and its homology
RP2:=ClosedSurface(-1);;   
RP2:=RegularCWComplex(RP2);;
RP2:=SimplifiedComplex(RP2);;
S10:=Suspension(RP2,10);;
Print(Size(S10),"\n");

Y:=ClosedSurface(-1);;      
Y:=RegularCWComplex(Y);;
Y:=SimplifiedComplex(Y);;
S100:=Suspension(Y,100);;
Print(List([0..102],i->Homology(S100,i),"\n"));

Y:=ClosedSurface(-1);;
Y:=RegularCWComplex(Y);;
Y:=SimplifiedComplex(Y);;
S1000:=Suspension(Y,1000);;
CriticalCells(S1000);;
Print(List([0..1002],i->Homology(S1000,i)));
