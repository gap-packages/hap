gap> START_TEST("HAP library");
gap> Y:=PoincareDodecahedronCWComplex(
> [[1,2,3,4,5],[6,7,8,9,10]],
> [[1,11,16,12,2],[19,9,8,18,14]],
> [[2,12,17,13,3],[20,10,9,19,15]],
> [[3,13,18,14,4],[16,6,10,20,11]],
> [[4,14,19,15,5],[17,7,6,16,12]],
> [[5,15,20,11,1],[18,8,7,17,13]]);;
gap> IsClosedManifold(Y);
true
gap> List([0..3],n->Homology(Y,n));
[ [ 0 ], [  ], [  ], [ 0 ] ]
gap> StructureDescription(FundamentalGroup(Y));
"SL(2,5)"
gap> STOP_TEST( "tst.tst", 1000 );
