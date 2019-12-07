#E
gap> START_TEST("HAP library");
gap> G:=SignedPermutationGroup(5);;
gap> P:=PolytopalComplex(G,[1,2,3,4,5]);;
gap> R:=FreeGResolution(P,5);;
gap> Homology(TensorWithIntegers(R),4);
[ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


