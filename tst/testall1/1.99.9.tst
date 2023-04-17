gap> START_TEST("HAP library");
gap> arc:=ArcPresentation(PureCubicalKnot(3,1));
[ [ 2, 5 ], [ 1, 3 ], [ 2, 4 ], [ 3, 5 ], [ 1, 4 ] ]
gap> S:=SphericalKnotComplement(arc);
Regular CW-complex of dimension 3

gap> S!.nrCells(3);
4
gap> Y:=ContractedComplex(S);
Regular CW-complex of dimension 2

gap> CriticalCells(Y);
[ [ 2, 1 ], [ 1, 9 ], [ 1, 11 ], [ 0, 22 ] ]
gap> G:=FundamentalGroup(Y);;
gap> AlexanderPolynomial(G);
x_1^2-x_1+1
gap> STOP_TEST( "tst.tst", 1000 );
