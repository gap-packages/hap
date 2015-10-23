gap> START_TEST("HAP library");
gap> Read(Concatenation(dir,"data245.txt"));
gap> G:=SymmetricMatrixToFilteredGraph(A,10,60);;
gap> K:=CliqueComplex(G,2);;
gap> P1:=PersistentBettiNumbers(K,1);;
gap> BarCodeCompactDisplay(P1);
gap> STOP_TEST( "tst.tst", 1000 );


