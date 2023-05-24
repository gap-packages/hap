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
gap> L:=LowIndexSubgroupsFpGroup(F,5);;
gap> SortedList(List(L,AbelianInvariants));
[ [ 2, 2, 2 ], [ 2, 2, 27 ], [ 2, 3 ], [ 3, 3 ], [ 3, 3, 16 ], [ 3, 4 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

