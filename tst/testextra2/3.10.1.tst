gap> START_TEST("HAP library");
gap> G:=SmallGroup(81,4);;
gap> A:=ModPSteenrodAlgebra(G,8);;
gap> List(ModPRingGenerators(A),x->Bockstein(A,x));
[ 0*v.1, 0*v.1, v.5, 0*v.1, (Z(3))*v.7+v.8+(Z(3))*v.9 ]
gap> STOP_TEST( "tst.tst", 1000 );
