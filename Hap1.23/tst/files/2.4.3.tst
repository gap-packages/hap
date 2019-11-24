gap> START_TEST("HAP library");
gap> 2simplices:=
> [[1,2,5], [2,5,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9],
> [4,5,8], [4,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7],
> [1,4,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
gap> K:=SimplicialComplex(2simplices);;
gap> C:=ChainComplex(K);;
gap> BettiNumber(C,1);
1
gap> BettiNumber(C,2);
0
gap> D:=TensorWithIntegersModP(C,2);;
gap> BettiNumber(D,1);
2
gap> BettiNumber(D,2);
1
gap> STOP_TEST( "tst.tst", 1000 );


