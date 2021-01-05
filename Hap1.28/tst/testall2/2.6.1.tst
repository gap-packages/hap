gap> START_TEST("HAP library");
gap> A:=List([1..10],i->1);;
gap> A:=List([1..10],i->StructuralCopy(A));;
gap> A:=List([1..10],i->StructuralCopy(A));;
gap> A[5][5][5]:=0;;
gap> M:=PurePermutahedralComplex(A);;
gap> Size(M);
999
gap> Y:=RegularCWComplex(M);;
gap> Size(Flat(Y!.boundaries));
125360
gap> STOP_TEST( "tst.tst", 1000 );


