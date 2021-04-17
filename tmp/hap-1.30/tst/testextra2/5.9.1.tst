gap> START_TEST("HAP library");
gap> Q:=QuadraticNumberField(-5);;
gap> OQ:=RingOfIntegers(Q);;
gap> I:=QuadraticIdeal(OQ,7);;
gap> R:=OQ mod I;;
gap> RightTransversal(I);;
gap> IsIntegralRing(R);
false
gap> STOP_TEST( "tst.tst", 1000 );
