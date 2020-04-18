gap> START_TEST("HAP library");
gap> R:=ResolutionPSL2QuadraticIntegers(-6,3);
Resolution of length 3 in characteristic 0 for PSL(2,O-6) . 
No contracting homotopy available. 

gap> G:=R!.group;;
gap> M:=HomogeneousPolynomials(G,2);;
gap> C:=HomToIntegralModule(R,M);;
gap> Cohomology(C,1);
[ 2, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );
