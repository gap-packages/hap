#D
gap> START_TEST("HAP library");
gap> M:=PurePermutahedralKnot(3,1);
Pure permutahedral complex of dimension 3.

gap> K:=Nerve(M);
Simplicial complex of dimension 1.

gap> #Display(M); #See Figure 1.10
gap> #Display(Graph(K)); #See Figure 1.12
gap> STOP_TEST( "tst.tst", 1000 );





