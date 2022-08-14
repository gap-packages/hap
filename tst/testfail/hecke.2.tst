gap> START_TEST("HAP library");
gap> gamma:=HAP_CongruenceSubgroupGamma0(4);;
gap> T1:=HomomorphismAsMatrix(HeckeOperatorWeight2(gamma,2,1));
[ [ 1, -1 ], [ 0, 0 ] ]
gap> STOP_TEST( "tst.tst", 1000 );

