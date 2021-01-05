gap> START_TEST("HAP library");
gap> C:=AutomorphismGroupAsCatOneGroup(DihedralGroup(8));
Cat-1-group with underlying group Group( [ f1, f2, f3, f4, f5, f6 ] ) . 

gap> Size(C);
64
gap> Q:=QuasiIsomorph(C);
Cat-1-group with underlying group Group( [ f6, f1*f2 ] ) . 

gap> Size(Q);
4
gap> N:=NerveOfCatOneGroup(Q,6);
Simplicial group of length 6

gap> K:=ChainComplexOfSimplicialGroup(N);
Chain complex of length 6 in characteristic 0 . 

gap> Homology(K,5);
[ 2, 2, 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );
