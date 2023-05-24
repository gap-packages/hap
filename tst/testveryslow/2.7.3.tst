gap> START_TEST("HAP library");
gap> S:=SimplicialComplex([[1,2],[2,3],[3,1]]);;
gap> S:=RegularCWComplex(S);;
gap> T:=DirectProduct(S,S,S,S);;
gap> Cohomology(T,1);
[ 0, 0, 0, 0 ]
gap> Cohomology(T,3);
[ 0, 0, 0, 0 ]
gap> Cohomology(T,4);
[ 0 ]
gap> cup:=CupProduct(T);
function( p, q, vv, ww ) ... end
gap> cup(3,1,[0,0,1,0],[0,0,1,0]);
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


