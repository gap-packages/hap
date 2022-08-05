gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/25spheres.txt");;
gap> Size(M);
499874
gap> A:=ContractibleSubcomplex(M);;
gap> Size(A);
499849
gap> P:=ExcisedPair(M,A);;
gap> Size(P[1]);
353
gap> Size(P[2]);
328
gap> C:=ChainComplexOfPair(P[1],P[2]);;
gap> Homology(C,1);
[  ]
gap> Homology(C,2);
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


