gap> START_TEST("HAP library");
gap> G:=AbelianGroup([4,8,9,27]);;
gap> R:=ResolutionAbelianGroup(G,5);;
gap> C:=TensorWithIntegers(R);;
gap> HomologyPrimePart(C,3,2);
[ 4, 4, 8 ]
gap> HomologyPrimePart(C,3,3);
[ 9, 9, 27 ]
gap> STOP_TEST( "tst.tst", 1000 );
