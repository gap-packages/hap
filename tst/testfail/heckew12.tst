gap> START_TEST("HAP library");
gap> H:=HAP_CongruenceSubgroupGamma0(1);;
gap> h:=HeckeOperator(H,2,12);
[ [ 2049, -7560, 0 ], [ 0, -24, 0 ], [ 0, 0, -24 ] ]
gap> STOP_TEST( "tst.tst", 1000 );
