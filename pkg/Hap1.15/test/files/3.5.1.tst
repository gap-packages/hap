gap> START_TEST("HAP library");
gap> R:=ResolutionPrimePowerGroup(DihedralGroup(128),10);;
gap> p:=PoincareSeries(R,10);
(1)/(x_1^2-2*x_1+1)
gap> STOP_TEST( "tst.tst", 1000 );


