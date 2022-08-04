gap> START_TEST("HAP library");
gap> gamma:=HAP_CongruenceSubgroupGamma0(39);;
gap> k:=2;; deg:=1;; 
gap> c:=CuspidalCohomologyHomomorphism(gamma,deg,k);;
gap> AbelianInvariants(Kernel(c));
[ 0, 0, 0, 0, 0, 0 ]
gap> STOP_TEST( "tst.tst", 1000 );

