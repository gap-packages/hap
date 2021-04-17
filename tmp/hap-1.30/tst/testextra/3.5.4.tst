gap> START_TEST("HAP library");
gap> G:=DihedralGroup(64);;
gap> FG:=GroupAlgebraAsFpGModule(G);;
gap> rad_1:=Radical(FG);;
gap> rad_2:=Radical(rad_1);;
gap> rad_3:=Radical(rad_2);;
gap> R:=ResolutionFpGModule(rad_3,11);;
gap> Homology(TensorWithIntegersModP(R,2),10);
12
gap> STOP_TEST( "tst.tst", 1000 );


