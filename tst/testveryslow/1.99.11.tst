gap> START_TEST("HAP library");
gap> K:=ComplexProjectiveSpace(2);;
gap> M:=WedgeSum(K,K);;
gap> Y:=RegularCWComplex(M);;
gap> cup:=CupProduct(Y);;
gap> cup(2,2,[1,0],[0,1]);
[ 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
