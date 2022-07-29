gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(MathieuGroup(12),2);;N:=Center(G);;
gap> RelativeGroupHomology(G,N,4);
[ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
