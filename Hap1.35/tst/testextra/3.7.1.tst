gap> START_TEST("HAP library");
gap> R:=ResolutionSL2Z(1,3);;
gap> IntegralCohomologyGenerators(R,0);
[ [ 1 ] ]
gap> IntegralCohomologyGenerators(R,1);
[  ]
gap> IntegralCohomologyGenerators(R,2);
[ [ 1 ] ]
gap> STOP_TEST( "tst.tst", 1000 );
