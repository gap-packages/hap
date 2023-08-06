gap> START_TEST("HAP library");
gap> K:=ClosedSurface(-1);;  
gap> K:=DirectProduct(K,K);;
gap> A:=CohomologyRing(K,2);;
gap> Bockstein(A,A.2);
v.4
gap> Bockstein(A,A.3);
v.4+v.6
gap> Bockstein(A,A.4);
0*v.1
gap> Bockstein(A,A.5);
v.7+v.8
gap> STOP_TEST( "tst.tst", 1000 );
