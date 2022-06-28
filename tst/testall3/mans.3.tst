gap> START_TEST("HAP library");
gap> ap:=ArcPresentation(PureCubicalKnot(3,1));;
gap> ap:=ArcPresentation(PureCubicalKnot(3,1));;
gap> ap:=ArcPresentation(PureCubicalKnot(3,1));;
gap> p:=7;;q:=6;;
gap> W:=ThreeManifoldViaDehnSurgery(ap,p,q);
Regular CW-complex of dimension 3

gap> Homology(W,0);Homology(W,1);Homology(W,2);Homology(W,3);
[ 0 ]
[ 6 ]
[  ]
[ 0 ]
gap> F:=FundamentalGroup(W);;
#I  there are 2 generators and 2 relators of total length 46
gap> L:=LowIndexSubgroupsFpGroup(F,5);;
gap> List(L,AbelianInvariants);
[ [ 2, 3 ], [ 3, 3 ], [ 2, 2, 2 ], [ 3, 4 ], [ 3, 3, 16 ], [ 2, 2, 27 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

