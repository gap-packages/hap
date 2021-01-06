gap> START_TEST("HAP library");
gap> S5:=SymmetricGroup(5);;SetName(S5,"S5");;
gap> S4:=SymmetricGroup(4);;SetName(S4,"S4");;
gap> A:=SymmetricGroup(3);;SetName(A,"S3");;
gap> AS5:=GroupHomomorphismByFunction(A,S5,x->x);;
gap> AS4:=GroupHomomorphismByFunction(A,S4,x->x);;
gap> D:=[S5,S4,[AS5,AS4]];;
gap> R:=ResolutionGraphOfGroups(D,6);;
gap> Homology(TensorWithIntegers(R),5);
[ 2, 2, 2, 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
