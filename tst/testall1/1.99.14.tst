gap> START_TEST("HAP library");
gap> D:=HomologicalGroupDecomposition(SymmetricGroup(5),2);;
gap> List(Flat(D),StructureDescription);
[ "D12", "S4", "C2 x C2" ]
gap> STOP_TEST( "tst.tst", 1000 );
