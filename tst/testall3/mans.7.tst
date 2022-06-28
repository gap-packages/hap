gap> START_TEST("HAP library");
gap> ap:=ArcPresentation(PureCubicalKnot(3,1));;
gap> W:=ThreeManifoldViaDehnSurgery(ap,1,5);;
gap> Homology(W,1);
[ 5 ]
gap> Homology(W,2);
[  ]
gap> Homology(W,3);;
gap> F:=FundamentalGroup(W);;
#I  there are 1 generator and 1 relator of total length 5
gap> StructureDescription(F);
"C5"
gap> STOP_TEST( "tst.tst", 1000 );

