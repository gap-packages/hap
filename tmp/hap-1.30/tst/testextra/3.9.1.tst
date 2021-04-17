gap> START_TEST("HAP library");
gap> G:=QuaternionGroup(8);;
gap> N:=Center(G);;
gap> Einf:=LHSSpectralSequenceLastSheet(G,N);
Graded algebra GF(2)[ x_1, x_2, x_3 ] / [ x_1^2+x_1*x_2+x_2^2, x_2^3 
 ] with indeterminate degrees [ 1, 1, 4 ]
gap> Mod2CohomologyRingPresentation(G,5);
Graded algebra GF(2)[ x_1, x_2, x_3 ] / [ x_1^2+x_1*x_2+x_2^2, x_2^3 
 ] with indeterminate degrees [ 1, 1, 4 ]
gap> STOP_TEST( "tst.tst", 1000 );
