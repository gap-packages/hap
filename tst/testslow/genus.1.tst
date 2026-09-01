gap> START_TEST("HAP library");
gap> K:=ContractibleGcomplex("BorelSerreSL2Z");;
gap> K:=BarycentricSubdivision(K);;
gap> dK:=ContractibleGcomplex("BorelSerreBoundarySL2Z");;
gap> dK:=BarycentricSubdivision(dK);;
gap> gamma:=HAP_CongruenceSubgroupGamma0(13);;
gap> Y:=GComplexToRegularCWComplex(K,gamma);;
gap> dY:=GComplexToRegularCWComplex(dK,gamma);;
gap> Homology(Y,1);
[ 0 ]
gap> EulerCharacteristic(Y) + Length(Homology(dY,1));
2
gap> STOP_TEST( "tst.tst", 1000 );
