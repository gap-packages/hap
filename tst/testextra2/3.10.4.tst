gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(3,1);;
gap> f:=KnotComplementWithBoundary(ArcPresentation(K));
Map of regular CW-complexes

gap> F:=FundamentalGroup(f);
[ f1, f2 ] -> [ f2^-1*f1*f2^2*f1*f2^-1, f1 ]
gap> phi:=ChainMap(f);
Chain Map between complexes of length 2 . 

gap> H:=Homology(phi,2);
[ g1 ] -> [ g1 ]
gap> STOP_TEST( "tst.tst", 1000 );
