#
gap> START_TEST("HAP library");
gap> Q:=QuaternionGroup(8);;
gap> B:=BarComplexOfMonoid(Q,4);;
gap> C:=ContractedComplex(B);;
gap> Homology(C,3);
[ 8 ]
gap> STOP_TEST( "tst.tst", 1000 );
