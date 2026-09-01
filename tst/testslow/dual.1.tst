gap> START_TEST("HAP library");
gap> K:=ContractibleGcomplex("SL(3,Z)a");;
gap> V:=DualComplex(K);;
gap> R:=FreeGResolution(V,1);; 
gap> G:=HAPCongruenceSubgroupGamma0(3,11);;
gap> S:=ResolutionFiniteSubgroup(R,G);; 
gap> C:=TensorWithIntegersSparse(S);;
gap> C:=ContractedComplex(C);;
gap> Homology(C,0);
[ 0, 0 ]
gap> R:=FreeGResolution(K,4);;
gap> S:=ResolutionFiniteSubgroup(R,G);;
gap> C:=TensorWithIntegersSparse(S);;
gap> C:=ContractedComplex(C);;
gap> Homology(C,3);
[ 2, 2, 12, 12, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
