gap> START_TEST("HAP library");
gap> 2simplices:= [[1,2,5], [2,5,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9], [4,5,8], [4,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7], [1,4,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
gap> K:=SimplicialComplex(2simplices);;
gap> C:=CochainComplex(K);;
gap> Cohomology(C,0);
[ 0 ]
gap> Cohomology(C,1);
[ 0 ]
gap> Cohomology(C,2);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


