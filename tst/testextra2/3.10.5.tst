gap> START_TEST("HAP library");
gap> G:=SmallGroup(32,23);;
gap> K:=CrossedModuleByNormalSubgroup(G,G);;
gap> T2:=NonabelianTensorProduct(K,K);;
gap> T3:=NonabelianTensorProduct(T2,K);;
gap> T:=Source(T3!.map);;
gap> Exponent(T);
4
gap> Order(T);
1073741824
gap> STOP_TEST( "tst.tst", 1000 );
