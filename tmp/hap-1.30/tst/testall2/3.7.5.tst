gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(SymmetricGroup(12),3);;
gap> x:=ModPCohomologyGenerators(G,5);
[ [ v.1, v.2, v.3, v.4, v.7, v.8, v.9, v.10, v.17, v.18, v.19, v.34, v.35, 
      v.59, v.61 ], function( x ) ... end ]
gap> List(x[1],x[2]);
[ 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5 ]
gap> STOP_TEST( "tst.tst", 1000 );
