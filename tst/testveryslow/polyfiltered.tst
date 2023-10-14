gap> START_TEST("HAP library");
gap> file:=HapFile("data253a.txt");;
gap> Read(file);;
gap> G:=SymmetricMatrixToFilteredGraph(D,100);;
gap> K:=FilteredRegularCWComplex(CliqueComplex(G,2));;
gap> K:=ContractedComplex(K);;
gap> P:=PersistentBettiNumbers(K,0)[13];
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 74, 74, 74, 74, 74, 74, 74, 74\
, 74, 74, 
  74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 
  74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 71, 37, 37, 37, 37, 37, 
  34, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]
gap> P:=PersistentBettiNumbers(K,1)[13];
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> L:=[];;
gap> for p in [1..14] do
> K:=RegularCWComplex(RandomSimplicialTwoComplex(50,p/1000));;
> h1:=Length(Homology(K,1));;
> Add(L, [1.0*(p/1000),h1,"red"]);
> od;
gap> L[13]{[1,3]};
[ 0.013, "red" ]
gap> M:=HenonOrbit([0,0],14/10,3/10,5*10^5,50,30);;
gap> Size(M);
611
gap> Size(ContractedComplex(M));
139
gap> STOP_TEST( "tst.tst", 1000 );
