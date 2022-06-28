gap> START_TEST("HAP library");
gap> ap1:=ArcPresentation(PureCubicalKnot(11,1));;
gap> A:=ThreeManifoldViaDehnSurgery(ap1,5,1);;
gap> ap2:=ArcPresentation(PureCubicalKnot(11,2));;
gap> B:=ThreeManifoldViaDehnSurgery(ap2,5,1);;
gap> W:=ConnectedSum(A,B);
Regular CW-complex of dimension 3

gap> Homology(W,1);
[ 2, 594 ]
gap> Homology(W,2);
[  ]
gap> Homology(W,3);
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );

