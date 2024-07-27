##########################################################
##########################################################
InstallOtherMethod(Radical,
"Radical of FpG-module",
[IsHapFPGModule],
function(M)
return RadicalOfFpGModule(M);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(RadicalSeries,
"Radical series of FpG-module",
[IsHapFPGModule],
function(M)
return RadicalSeriesOfFpGModule(M);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RadicalSeries,
"Radical series of FpG-resolution",
[IsHapResolution],
function(R)
return RadicalSeriesOfResolution(R);
end);
##########################################################
##########################################################




##########################################################
##########################################################
InstallOtherMethod(SimplifiedComplex,
"simplify a free ZG-resolution for a group G",
[IsHapResolution],
function(R)
return TietzeReducedResolution(R);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(SimplifiedComplex,
"simplify a chain complex of free abelian groups",
[IsHapChainComplex],
function(R)
return CoreducedChainComplex(R);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(Resolution,
"free ZG-resolution for a group G",
[IsGroup,IsInt],
function(G,n)
return ResolutionGenericGroup(G,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(CupProduct,
"integral cohomology cup product for regular CW-spaces",
[IsHapRegularCWComplex],
function(Y)
return CupProductOfRegularCWComplex(Y);
end);
##########################################################
##########################################################
#
###########################################################
##########################################################
InstallOtherMethod(CupProduct,
"Mod p cohomology cup product for regular CW-spaces",
[IsHapRegularCWComplex,IsInt],
function(Y,p)
return CupProductOfRegularCWComplexModP(Y,p);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(CupProduct,
"integral cohomology cup product for 2-complex of group presentation",
[IsFpGroup],
function(G)
return HAP_CupProductOfPresentation(G);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(CupProduct,
"integral cohomology cup product for a group",
[IsHapResolution, IsInt, IsInt,IsList,IsList],
function(R,p,q,v,w)
return IntegralCupProduct(R,v,w,p,q);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(Size,
"size of a chain complex = number of generators in it",
[IsHapChainComplex],
function(C)
local n, sz;
sz:=0;
for n in [0..Length(C)] do
sz:=sz+C!.dimension(n);
od;
return sz;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(ChainMap,
"chain map of regular cw-map ",
[IsHapRegularCWMap],
function(f)
return CWMap2ChainMap(f);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ChainMap,
"chain map of simplicial map ",
[IsHapSimplicialMap],
function(f)
return ChainMapOfSimplicialMap(f);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ChainMap,
"chain map of a pure cubical pairs",
[IsHapPureCubicalComplex, IsHapPureCubicalComplex, IsHapPureCubicalComplex, IsHapPureCubicalComplex],
function(X,A,Y,B)
return ChainMapOfCubicalPairs(X,A,Y,B);
end);
##########################################################
##########################################################




##########################################################
##########################################################
InstallMethod(ExcisedPair,
"excision for pure cubical complexes ",
[IsHapPureCubicalComplex, IsHapPureCubicalComplex],
function(X,A)
return ExcisedPureCubicalPair(X,A);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(FiltrationTerm,
"Term of filtered pure cubical complex ",
[IsHapPureCubicalComplex, IsInt],
function(Y,n)
return FiltrationTermOfPureCubicalComplex(Y,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(FiltrationTerm,
"Term of filtered graph ",
[IsHapFilteredGraph, IsInt],
function(Y,n)
return FiltrationTermOfGraph(Y,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(FiltrationTerm,
"Term of filtered regular CW-complex ",
[IsHapFilteredRegularCWComplex, IsInt],
function(Y,n)
return FiltrationTermOfRegularCWComplex(Y,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(PathComponent,
"path component of simplicial complex ",
[IsHapSimplicialComplex, IsInt],
function(Y,n)

return PathComponentOfSimplicialComplex(Y,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(PathComponent,
"path component of pure cubical complex ",
[IsHapPureCubicalComplex, IsInt],
function(Y,n)
return PathComponentOfPureComplex(Y,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PathComponent,
"path component of pure permutahedral complex ",
[IsHapPurePermutahedralComplex, IsInt],
function(Y,n)

return PathComponentOfPureComplex(Y,n);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(FundamentalGroup,
"fundamental group of regular CW space without Tietze reductions ",
[IsHapRegularCWComplex, IsString],
function(Y,s)

return FundamentalGroupOfRegularCWComplex(Y,1,s);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(Pushout,
"amalgamated sum of fp groups ",
[IsGroupHomomorphism, IsGroupHomomorphism],
function(f,g)

return PushoutOfFpGroups(f,g);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(Dimensions,
"array dimensions of a pure cubical complex ",
[IsHapPureCubicalComplex],
function(M)

return EvaluateProperty(M,"arraySize");
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Dimensions,
"array dimensions of a pure permutahedral complex ",
[IsHapPurePermutahedralComplex],
function(M)

return EvaluateProperty(M,"arraySize");
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(RegularCWMap,
"inclusion of pure cubical complexes to regular CW map ",
[IsHapPureCubicalComplex, IsHapPureCubicalComplex],
function(M,A)

return HAP_PureCubicalPairToCWMap(M,A);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWMap,
"inclusion of simplicial complexes to regular CW map ",
[IsHapSimplicialComplex, IsHapSimplicialComplex],
function(M,A)

return HAP_SimplicialPairToCWMap(M,A);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(PureComplexSubcomplex,
"subcomplex of  pure cubical complex ",
[IsHapPureCubicalComplex, IsList],
function(M,L)

return HAP_PureComplexSubcomplex(M,L);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PureComplexSubcomplex,
"subcomplex of  pure permutahedral complex ",
[IsHapPurePermutahedralComplex, IsList],
function(M,L)

return HAP_PureComplexSubcomplex(M,L);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallGlobalFunction(HAP_PureComplexSubcomplex,
function(M,L)
local A, dim,  ArrayValueDim, ArrayAssignDim, 
ArrayIt,  Fun;

A:=0*M!.binaryArray;

dim:=Dimension(M);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);

#########################
Fun:=function(i);
if ArrayValueDim(M!.binaryArray,i)=1 then
ArrayAssignDim(A,i,1);
fi;
end;
#########################
ArrayIt(L,Fun);

return PureComplex(M,A);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(PureComplexMeet,
"meet of two pure cubical complexes ",
[IsHapPureCubicalComplex, IsHapPureCubicalComplex],
function(M,N)
local T;

if not 
EvaluateProperty(M,"arraySize") = EvaluateProperty(N,"arraySize")
then Print("Complexes should have the same array size.\n");
return fail;
fi;
 
T:=PureComplexThickened(N);       #This is very inefficient!
T:=PureComplexIntersection(T,M);
return T;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PureComplexMeet,
"meet of two pure permutahedral complexes ",
[IsHapPurePermutahedralComplex, IsHapPurePermutahedralComplex],
function(M,N)
local T;

if not 
EvaluateProperty(M,"arraySize") = EvaluateProperty(N,"arraySize")
then Print("Complexes should have the same array size.\n");
return fail;
fi;

T:=PureComplexThickened(N);       #This is very inefficient!
T:=PureComplexIntersection(T,M);
return T;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(PureComplexRandomCell,
"random cell of a pure cubical complex ",
[IsHapPureCubicalComplex],
function(K);
return RandomCellOfPureComplex(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PureComplexRandomCell,
"random cell of a pure cubical complex ",
[IsHapPurePermutahedralComplex],
function(K);
return RandomCellOfPureComplex(K);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(ContractibleSubcomplex,
"contractible subcomplex of pure cubical complex ",
[IsHapPureCubicalComplex],
function(K)
local A;
A:=HomotopyEquivalentLargerSubArray(K!.binaryArray,0*K!.binaryArray);

return PureCubicalComplex(A);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractibleSubcomplex,
"contractible subcomplex of pure permutahedral complex ",
[IsHapPurePermutahedralComplex],
function(K)
local A;

A:=HomotopyEquivalentLargerSubPermArray(K!.binaryArray,0*K!.binaryArray);

return PurePermutahedralComplex(A);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(ContractibleSubcomplex,
"contractible subcomplex of simplicial complex ",
[IsHapSimplicialComplex],
function(K);
return ContractibleSubcomplexOfSimplicialComplex(K);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(BoundaryMap,
"inclusion of the boundary of a regular CW complex ",
[IsHapRegularCWComplex],
function(K);
return BoundaryPairOfPureRegularCWComplex(K);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(CriticalCells,
"critical cells of regular CW complex",
[IsHapRegularCWComplex],
function(K);
if IsBound(K!.directed) then
return HAP_CriticalCellsDirected(K);
fi;

return CriticalCellsOfRegularCWComplex(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(DisplayArcPresentation,
"displays a pure cubical knot",
[IsHapPureCubicalComplex],
function(K);
ViewPureCubicalKnot(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(KnotReflection,
"reflects a pure cubical knot",
[IsHapPureCubicalComplex],
function(K);
return ReflectedCubicalKnot(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(DirectProductOp,
"direct product of  pure cubical complexes",
[IsList, IsGroupHomomorphism],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return DirectProductOfGroupHomomorphisms(L[1],L[2]);
else
D:=DirectProductOp(L{[2..Length(L)]},K);
return  DirectProductOfGroupHomomorphisms(L[1],D);
fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(DirectProductOp,
"direct product of  pure cubical complexes",
[IsList, IsHapPureCubicalComplex],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return DirectProductOfPureCubicalComplexes(L[1],L[2]);
else
D:=DirectProductOp(L{[2..Length(L)]},K);
return  DirectProductOfPureCubicalComplexes(L[1],D);
fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(DirectProductOp,
"direct product of  regular CW-complexes",
[IsList, IsHapRegularCWComplex],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return DirectProductOfRegularCWComplexes(L[1],L[2]);
else
D:=DirectProductOp(L{[2..Length(L)]},K);
return  DirectProductOfRegularCWComplexes(L[1],D);
fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(DirectProductOp,
"direct product of  simplicial complexes",
[IsList, IsHapSimplicialComplex],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return DirectProductOfSimplicialComplexes(L[1],L[2]);
else
D:=DirectProductOp(L{[2..Length(L)]},K);
return  DirectProductOfSimplicialComplexes(L[1],D);
fi;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(TensorProductOp,
"tensor product of free resolutions",
[IsList, IsHapResolution],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return ResolutionDirectProduct(L[1],L[2]);
else
D:=TensorProductOp(L{[2..Length(L)]},K);
return  ResolutionDirectProduct(L[1],D);
fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(TensorProductOp,
"tensor product of chain complex",
[IsList, IsHapChainComplex],
function(L,K)
local D;
if Length(L)=1 then return L[1]; fi;
if Length(L)=2 then
return TensorProductOfChainComplexes(L[1],L[2]);
else
D:=TensorProductOp(L{[2..Length(L)]},K);
return  TensorProductOfChainComplexes(L[1],D);
fi;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(ZigZagContractedComplex,
"zig zag contracted  pure cubical complex",
[IsHapPureCubicalComplex],
function(K);
return ZigZagContractedPureComplex(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ZigZagContractedComplex,
"zig zag contracted  pure cubical complex",
[IsHapPureCubicalComplex,IsInt],
function(K,n);
return ZigZagContractedPureComplex(K,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(ZigZagContractedComplex,
"zig zag contracted  pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(K);
return ZigZagContractedPureComplex(K);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ZigZagContractedComplex,
"zig zag contracted  pure permutahedral complex",
[IsHapPurePermutahedralComplex,IsInt],
function(K,n);
return ZigZagContractedPureComplex(K,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(ZigZagContractedComplex,
"zig zag contracted  filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex],
function(K); 
return ZigZagContractedFilteredPureCubicalComplex(K);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  graph",
[IsHapGraph],
function(G)
local K,v;
v:=EvaluateProperty(G,"numberofvertices");

K:=Objectify(HapGraph,
                rec(
                incidenceMatrix:=StructuralCopy(G!.incidenceMatrix),
                properties:=
                [
                ["numberofvertices",v]
                ]));

ContractGraph(K);

return K;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  simplicial  complex",
[IsHapSimplicialComplex],
function(Y)
local K, Vertices, NrSimplices,Simplices,SimplicesLst,EnumeratedSimplex, dim;

Vertices:=StructuralCopy(Y!.vertices);
SimplicesLst:=StructuralCopy(Y!.simplicesLst);
dim:=Dimension(Y);                       

#####################
NrSimplices:=function(n);
if n<0 or n>dim then return 0; fi;
return Length(SimplicesLst[n+1]);
end;
#####################

#####################
Simplices:=function(n,i);
return SimplicesLst[n+1][i];
end;
#####################


#####################
EnumeratedSimplex:=function(v);
return PositionSet(SimplicesLst[Length(v)],v);
end;
#####################

K:=Objectify(HapSimplicialComplex,
           rec(
           vertices:=Vertices,
           nrSimplices:=NrSimplices,
           simplices:=Simplices,
           simplicesLst:=SimplicesLst,
           enumeratedSimplex:=EnumeratedSimplex,
           properties:=[
           ["dimension",dim]]
           ));

ContractSimplicialComplex(K);

dim:=PositionProperty(SimplicesLst,x->Size(x)=0);
if dim=fail then dim:=Length(SimplicesLst)-1; else dim:=dim-2; fi;

return K;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  filtered regular CW complex",
[IsHapFilteredRegularCWComplex],
function(Y);
return ContractedFilteredRegularCWComplex(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  pure cubical complex",
[IsHapPureCubicalComplex,IsHapPureCubicalComplex],
function(M,S)
local A;

A:=HomotopyEquivalentSmallerSubArray(M!.binaryArray,S!.binaryArray);
return PureCubicalComplex(A);;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ExpandedComplex,
"Homotopy expanded  pure cubical complex",
[IsHapPureCubicalComplex,IsHapPureCubicalComplex],
function(M,S)
local A;

A:=HomotopyEquivalentLargerSubArray(S!.binaryArray,M!.binaryArray);
return PureCubicalComplex(A);;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ExpandedComplex,
"Homotopy expanded  pure peermutahedral complex",
[IsHapPurePermutahedralComplex,IsHapPurePermutahedralComplex],
function(M,S)
local A;

A:=HomotopyEquivalentLargerSubPermArray(S!.binaryArray,M!.binaryArray);
return PurePermutahedralComplex(A);;
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(ExpandedComplex,
"Homotopy expanded  pure cubical complex",
[IsHapPureCubicalComplex],
function(M)
local A,S;

S:=0*M!.binaryArray;
S:=S+1;
A:=HomotopyEquivalentLargerSubArray(S,M!.binaryArray);
return PureCubicalComplex(A);;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ExpandedComplex,
"Homotopy expanded  pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(M)
local A,S;

S:=0*M!.binaryArray;
S:=S+1;
A:=HomotopyEquivalentLargerSubPermArray(S,M!.binaryArray);
return PurePermutahedralComplex(A);;
end);
##########################################################
##########################################################




##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  pure permutahedral complex",
[IsHapPurePermutahedralComplex,IsHapPurePermutahedralComplex],
function(M,S);
return HomotopyEquivalentMinimalPureSubcomplex(M,S);;
end);
##########################################################
##########################################################





##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex],
function(Y);
return
ContractedFilteredPureCubicalComplex(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(SimplifiedComplex,
"simplify a regular CW-complex while retaining the homeomorphism type",
[IsHapRegularCWComplex],
function(Y);
return
SimplifiedRegularCWComplex(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(SimplifiedComplex,
"simplify a pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(Y);
return
ZigZagContractedComplex(Y);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted sparse chain complex",
[IsHapSparseChainComplex],
function(C);
return
SimplifiedSparseChainComplex(C);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  chain complex",
[IsHapChainComplex],
function(C)
local D;
D:=ChainComplexToSparseChainComplex(C);
D:=ContractedComplex(D);
return
SparseChainComplexToChainComplex(D);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  free ZG-resolution",
[IsHapResolution],
function(R)
return
TietzeReducedResolution(R);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(ContractedComplex,
"contracted  chain complex",
[IsHapCochainComplex],
function(C)
local D;
D:=HomToIntegers(C);
D:=ContractedComplex(D);
D:=HomToIntegers(D);
return D;
end);
##########################################################
##########################################################





##########################################################
##########################################################
InstallMethod(ConcentricFiltration,
"Concentric filtration on pure cubical complex",
[IsHapPureCubicalComplex,IsInt],
function(Y,n);
return 
ConcentricallyFilteredPureCubicalComplex(Y,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsInt],
function(Y,n)
local W;

if n=0 and EvaluateProperty(Y,"monotonic")=true then
return HAP_BettiZeroMonotonic(Y);
fi;

W:=ContractedFilteredPureCubicalComplex(Y);
W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
if Dimension(W)<=3 then 
return PersistentBettiNumbers(W,n,2);
else
return PersistentBettiNumbers(W,n);
fi;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsInt,IsInt],
function(Y,n,prime)
local W;
if n=0 and EvaluateProperty(Y,"monotonic")=true then
return HAP_BettiZeroMonotonic(Y);
fi;

W:=ContractedFilteredPureCubicalComplex(Y);
W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbers(W,n,prime);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsInt,IsInt],
function(Y,n,prime)
local W;
if n=0 and EvaluateProperty(Y,"monotonic")=true then
return HAP_BettiZeroMonotonic(Y);
fi;

W:=ContractedFilteredPureCubicalComplex(Y);
W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbersAlt(W,n,prime);
end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsInt],
function(Y,n)
local W;

if n=0 and EvaluateProperty(Y,"monotonic")=true then
return HAP_BettiZeroMonotonic(Y);
fi;

W:=ContractedFilteredPureCubicalComplex(Y);
W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbersAlt(W,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsInt,IsBool],
function(Y,n,bool)
local W;

if n=0 and (not bool) and EvaluateProperty(Y,"monotonic")=true then
return HAP_BettiZeroMonotonic(Y);
fi;

W:=ContractedFilteredPureCubicalComplex(Y);
W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbersAlt(W,n,bool);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsList],
function(Y,N)
local W,P, NN;

W:=ContractedFilteredPureCubicalComplex(Y);

if 0 in N and EvaluateProperty(Y,"monotonic")=true then
P:=[ HAP_BettiZeroMonotonic(Y) ];
NN:=Filtered(N,i->not i=0);
Append(P,PersistentBettiNumbersAlt(W,NN));
return P;
fi;

W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbersAlt(W,N);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered pure cubical complex",
[IsHapFilteredPureCubicalComplex,IsList,IsBool],
function(Y,N,bool)
local W,P, NN;

W:=ContractedFilteredPureCubicalComplex(Y);

if 0 in N and EvaluateProperty(Y,"monotonic")=true and (not bool) then
P:=[ HAP_BettiZeroMonotonic(Y) ];
NN:=Filtered(N,i->not i=0);
Append(P,PersistentBettiNumbersAlt(W,NN));
return P;
fi;

W:= FilteredPureCubicalComplexToCubicalComplex(W);
W:=FilteredCubicalComplexToFilteredRegularCWComplex(W);
return PersistentBettiNumbersAlt(W,N,bool);
end);
##########################################################
##########################################################






##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered simplicial complex",
[IsHapFilteredSimplicialComplex,IsInt],
function(Y,n)
local W;
W:=FilteredRegularCWComplex(Y);
return PersistentBettiNumbers(W,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered simplicial complex",
[IsHapFilteredSimplicialComplex,IsInt],
function(Y,n)
local W;
W:=FilteredRegularCWComplex(Y);
return PersistentBettiNumbersAlt(W,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered simplicial complex",
[IsHapFilteredSimplicialComplex,IsInt,IsInt],
function(Y,n,prime)
local W;
W:=FilteredRegularCWComplex(Y);
return PersistentBettiNumbers(W,n,prime);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered simplicial complex",
[IsHapFilteredSimplicialComplex,IsInt,IsInt],
function(Y,n,prime)
local W;
W:=FilteredRegularCWComplex(Y);
return PersistentBettiNumbersAlt(W,n,prime);
end);
##########################################################
##########################################################





##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt],
function(Y,n)
local W;
#W:=ContractedFilteredRegularCWComplex(Y);
W:=Y;
W:=SparseChainComplexOfFilteredRegularCWComplex(W);
return PersistentBettiNumbers(W,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt,IsInt],
function(Y,n,prime)
local W;
#W:=ContractedFilteredRegularCWComplex(Y);
W:=Y;
W:=SparseChainComplexOfFilteredRegularCWComplex(W);
return PersistentBettiNumbers(W,n,prime);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt,IsInt],
function(Y,n,prime);
return PersistentBettiNumbersViaContractions(Y,n,prime);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt,IsInt,IsBool],
function(Y,n,prime,bool);
return PersistentBettiNumbersViaContractions(Y,n,prime,bool);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt],
function(Y,n);
Print("Using homology over GF(2).\n");
return PersistentBettiNumbersViaContractions(Y,n,2);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsInt,IsBool],
function(Y,n,bool);
Print("Using homology over GF(2).\n");
return PersistentBettiNumbersViaContractions(Y,n,2,bool);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsList],
function(Y,N);
Print("Using homology over GF(2).\n");
return PersistentBettiNumbersViaContractions(Y,N,2);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbersAlt,
"Betti number of a filtered regular CW-complex",
[IsHapFilteredRegularCWComplex,IsList,IsBool],
function(Y,N,bool);
Print("Using homology over GF(2).\n");
return PersistentBettiNumbersViaContractions(Y,N,2,bool);
end);
##########################################################
##########################################################









##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a sparse filtered chain complex",
[IsHapFilteredSparseChainComplex,IsInt],
function(C,n);
return PersistentHomologyOfFilteredSparseChainComplex(C,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a sparse filtered chain complex",
[IsHapFilteredSparseChainComplex,IsInt,IsInt],
function(C,n,prime);
return PersistentHomologyOfFilteredSparseChainComplex(
                   TensorWithIntegersModPSparse(C,prime),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a  filtered chain complex",
[IsHapFilteredChainComplex,IsInt],
function(C,n);
return PersistentBettiNumbers(
FilteredChainComplexToFilteredSparseChainComplex(C),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(PersistentBettiNumbers,
"Betti number of a  filtered chain complex",
[IsHapFilteredChainComplex,IsInt,IsInt],
function(C,n,prime);
return PersistentBettiNumbers(
FilteredChainComplexToFilteredSparseChainComplex(C),n,prime);
end);
##########################################################
##########################################################





##########################################################
##########################################################
InstallOtherMethod(Size,
"Size of an integer",
[IsInt],
function(n);
return n;
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a chain complex",
[IsHapChainComplex,IsInt],
function(C,n);
if IsPrimeInt(EvaluateProperty(C,"characteristic")) then
return Size(Homology(C,n));
fi;
if EvaluateProperty(C,"characteristic")=0 then
return Size(Homology(TensorWithRationals(C),n));
fi;
return Size(Homology(TensorWithRationals(C),n));

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a sparse chain complex",
[IsHapSparseChainComplex,IsInt],
function(C,n);
return Bettinumbers(C,n);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a pure cubical complex",
[IsHapPureCubicalComplex, IsInt],
function(S,n);
if n=0 then
return PathComponentOfPureComplex(S,0); 
fi;
if n<0 or n>Dimension(S) then return 0; fi;
if Dimension(S)<=3 then
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n,2);
fi;
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a pure cubical complex in characteristic p",
[IsHapPureCubicalComplex, IsInt, IsInt],
function(S,n,p);
if n=0 then
return PathComponentOfPureComplex(S,0);
fi;
if n<0 or n>Dimension(S) then return 0; fi;
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n,p);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a pure permutahedral complex",
[IsHapPurePermutahedralComplex,IsInt],
function(S,n);
if n=0 then
return PathComponentOfPureComplex(S,0);       
fi;
if n<0 or n>Dimension(S) then return 0; fi;
if Dimension(S)<=3 then
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n,2);
fi;
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a pure permutahedral complex in characteristic p",
[IsHapPurePermutahedralComplex,IsInt,IsInt],
function(S,n,p);
if n=0 then
return PathComponentOfPureComplex(S,0);      
fi;
if n<0 or n>Dimension(S) then return 0; fi;
return BettiNumber(RegularCWComplex(ContractedComplex(S)),n,p);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a cubical complex",
[IsHapCubicalComplex,IsInt],
function(S,n);
return BettiNumber(RegularCWComplex(S),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a cubical complex in characteristic p",
[IsHapCubicalComplex,IsInt,IsInt],
function(S,n,p);
return BettiNumber(RegularCWComplex(S),n,p);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(S,n)
local C;
if n=0 then
return PathComponentOfSimplicialComplex(S,0);
fi;
return BettiNumber(RegularCWComplex(S),n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a simplicial complex in characteristic p",
[IsHapSimplicialComplex,IsInt,IsInt],
function(S,n,p)
local C;
if n=0 then
return PathComponentOfSimplicialComplex(S,0);
fi;
return BettiNumber(RegularCWComplex(S),n,p);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallMethod(BettiNumber,
"Betti number of a regular complex",
[IsHapRegularCWComplex,IsInt],
function(S,n)
local C;
if n=0 then
return Size(PiZero(S)[1]);
fi;
if n<0 or n>Dimension(S) then return 0; fi;

C:=TruncatedRegularCWComplex(S,n+1);
C:=ChainComplex(C);
C:=TensorWithRationals(C);
return Size(Homology(C,n));

end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(BettiNumber,
"Betti number of a regular complex in characteristic p",
[IsHapRegularCWComplex,IsInt,IsInt],
function(S,n,p)
local C;
if n=0 then
return Size(PiZero(S)[1]);
fi;
if n<0 or n>Dimension(S) then return 0; fi;

C:=TruncatedRegularCWComplex(S,n+1);
C:=ChainComplex(C);
if p>0 then 
C:=TensorWithIntegersModP(C,p);
return Size(Homology(C,n));
fi;

if p=0 then 
C:=TensorWithRationals(C);
return Size(Homology(C,n));
fi;

end);
##########################################################
##########################################################



##########################################################
##########################################################
InstallMethod(FilteredRegularCWComplex,
"filtered simplicial complex to filtered regular CW-complex",
[IsHapFilteredSimplicialComplex],
function(S)
return FilteredSimplicialComplexToFilteredCWComplex(S);
end);
##########################################################
##########################################################





##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"simplicial complex to regular CW-complex",
[IsHapSimplicialComplex],
function(S)
return SimplicialComplexToRegularCWComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(RegularCWComplex,
"List of boundaries to regular CW-omplex",
[IsList],
function(S)
return HAPRegularCWComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"List of boundaries to regular CW-omplex",
[IsList,IsList],
function(S,T)
return HAPRegularCWComplex(S,T);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"List of boundaries to regular CW-omplex",
[IsList,IsPseudoList],
function(S,T)
return HAPRegularCWComplex(S,T);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"pure cubical complex to regular CW-complex",
[IsHapPureCubicalComplex],
function(S)
return CubicalComplexToRegularCWComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"cubical complex to regular CW-complex",
[IsHapCubicalComplex],
function(S)
return CubicalComplexToRegularCWComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWComplex,
"pure permutahedral complex to regular CW-complex",
[IsHapPurePermutahedralComplex],
function(S)
return PermutahedralComplexToRegularCWComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(RegularCWPolytope,
"convex hull of set of points",
[IsList],
function(S)
return HAPRegularCWPolytope(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(RegularCWPolytope,
"convex hull of a group's orbit on a vector",
[IsPermGroup,IsList],
function(G,v)
return RegularCWOrbitPolytope(G,v);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(SimplicialComplex,
"simplicial complex generated by a set of simplices",
[IsList],
function(S)
return MaximalSimplicesToSimplicialComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Nerve,
"Nerve of a pure cubical complex",
[IsHapPureCubicalComplex],
function(S)
return PureComplexToSimplicialComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Nerve,
"n-skeleton of nerve of a pure cubical complex",
[IsHapPureCubicalComplex, IsInt],
function(S,n)
return PureComplexToSimplicialComplex(S,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Nerve,
"Nerve of a pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(S)
return PureComplexToSimplicialComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Nerve,
"n-skeleton of nerve of a pure permutahedral complex",
[IsHapPurePermutahedralComplex, IsInt],
function(S,n)
return PureComplexToSimplicialComplex(S,n);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallMethod(Graph,
"1-skeleton of a simplicial complex returned as a HAP graph",
[IsHapSimplicialComplex],
function(S)
return GraphOfSimplicialComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Graph,
"1-skeleton of a regular CW-complex returned as a HAP graph",
[IsHapRegularCWComplex],
function(S)
return GraphOfRegularCWComplex(S);
end);
##########################################################
##########################################################


##########################################################
##########################################################
InstallOtherMethod(Display,
"Display a graph",
[IsHapGraph],
function(S)
GraphDisplay(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Display,
"Display a 3-dimensional pure permutahedral complex",
[IsHapPurePermutahedralComplex],
function(S)
ViewPureComplex(S);
end);
##########################################################
##########################################################

##########################################################
##########################################################
InstallOtherMethod(Display,
"Display a 2- or 3-dimensional pure cubical complex",
[IsHapPureCubicalComplex],
function(S)
ViewPureComplex(S);
end);
##########################################################
##########################################################






