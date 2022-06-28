gap> START_TEST("HAP library");
gap> gamma:=HAP_CongruenceSubgroupGamma0(11);;
gap> T1:=HomomorphismAsMatrix(HeckeOperatorWeight2(gamma,1,1));
[ [ 1, 0, 0 ], [ 0, 1, 0 ], [ 0, 0, 1 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

