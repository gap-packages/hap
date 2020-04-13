#D
gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(3,1);
prime knot 1 with 3 crossings

gap> L:=KnotReflection(K);
Reflected( prime knot 1 with 3 crossings )

gap> S:=KnotSum(K,L);
prime knot 1 with 3 crossings + Reflected( prime knot 1 with 3 crossings )

gap> #DisplayArcPresentation(S);
gap> STOP_TEST( "tst.tst", 1000 );


