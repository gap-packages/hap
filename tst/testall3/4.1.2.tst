gap> START_TEST("HAP library");
gap> SCTL:=EmptySCTable(13,0,"antisymmetric");;
gap> SetEntrySCTable( SCTL, 1, 6, [ 1, 7 ] );;
gap> SetEntrySCTable( SCTL, 1, 8, [ 1, 9 ] );;
gap> SetEntrySCTable( SCTL, 1, 10, [ 1, 11 ] );;
gap> SetEntrySCTable( SCTL, 1, 12, [ 1, 13 ] );;
gap> SetEntrySCTable( SCTL, 1, 7, [ -1, 6 ] );;
gap> SetEntrySCTable( SCTL, 1, 9, [ -1, 8 ] );;
gap> SetEntrySCTable( SCTL, 1, 11, [ -1, 10 ] );;
gap> SetEntrySCTable( SCTL, 1, 13, [ -1, 12 ] );;
gap> SetEntrySCTable( SCTL, 6, 7, [ 1, 2 ] );;
gap> SetEntrySCTable( SCTL, 8, 9, [ 1, 3 ] );;
gap> SetEntrySCTable( SCTL, 6, 9, [ -1, 5 ] );;
gap> SetEntrySCTable( SCTL, 7, 8, [ 1, 5 ] );;
gap> SetEntrySCTable( SCTL, 2, 8, [ 1, 12 ] );;
gap> SetEntrySCTable( SCTL, 2, 9, [ 1, 13 ] );;
gap> SetEntrySCTable( SCTL, 3, 6, [ 1, 10 ] );;
gap> SetEntrySCTable( SCTL, 3, 7, [ 1, 11 ] );;
gap> SetEntrySCTable( SCTL, 2, 3, [ 1, 4 ] );;
gap> SetEntrySCTable( SCTL, 5, 6, [ -1, 12 ] );;
gap> SetEntrySCTable( SCTL, 5, 7, [ -1, 13 ] );;
gap> SetEntrySCTable( SCTL, 5, 8, [ -1, 10 ] );;
gap> SetEntrySCTable( SCTL, 5, 9, [ -1, 11 ] );;
gap> SetEntrySCTable( SCTL, 6, 11, [ -1/2, 4 ] );;
gap> SetEntrySCTable( SCTL, 7, 10, [ 1/2, 4 ] );;
gap> SetEntrySCTable( SCTL, 8, 13, [ 1/2, 4 ] );;
gap> SetEntrySCTable( SCTL, 9, 12, [ -1/2, 4 ] );;
gap> L:=LieAlgebraByStructureConstants(Rationals,SCTL);;
gap> LieAlgebraHomology(L,2);
2
gap> LeibnizAlgebraHomology(L,2);
6
gap> 
gap> C:=Source(LieCoveringHomomorphism(L));;
gap> LieAlgebraHomology(C,2);
0
gap> LeibnizAlgebraHomology(C,2);
6
gap> Dimension(LieEpiCentre(L));
1
gap> D:=Source(LeibnizQuasiCoveringHomomorphism(L));;
gap> LeibnizAlgebraHomology(D,2);
6
gap> Source(LieExteriorSquare(L).homomorphism);
<algebra of dimension 14 over Rationals>
gap> STOP_TEST( "tst.tst", 1000 );
