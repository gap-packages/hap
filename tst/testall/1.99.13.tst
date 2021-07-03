gap> START_TEST("HAP library");
gap> K:=ClosedSurface(2);;
gap> A:=CohomologyRing(K,2);
<algebra of dimension 6 over GF(2)>
gap> Basis(A)[2]*Basis(A)[5];
v.6
gap> STOP_TEST( "tst.tst", 1000 );
