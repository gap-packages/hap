#D
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"data1V2X.pdb");;
gap> K:=ReadPDBfileAsPurePermutahedralComplex(dir);
Reading chain containing 191 atoms.
Pure permutahedral complex of dimension 3.

gap> #Display(ZigZagContractedComplex(K));
gap> #See Figure 1.24 (right)
gap> C:=PureComplexComplement(K);
Pure permutahedral complex of dimension 3.

gap> Size(C);
418922
gap> D:=ZigZagContractedComplex(C);
Pure permutahedral complex of dimension 3.

gap> Size(D);
3438
gap> Y:=RegularCWComplex(D);
Regular CW-complex of dimension 3

gap> #W:=SimplifiedComplex(Y);
gap> #W!.nrCells(3);
gap> W:=ContractedComplex(Y);
Regular CW-complex of dimension 2

gap> Size(W);
69731
gap> CriticalCells(W);
[ [ 2, 1 ], [ 2, 179 ], [ 1, 13092 ], [ 1, 14050 ], [ 0, 16453 ] ]
gap> STOP_TEST( "tst.tst", 1000 );


