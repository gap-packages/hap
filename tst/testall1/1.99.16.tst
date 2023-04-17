gap> START_TEST("HAP library");
gap> L:=[ [ 1, 4, 5 ], [ 2, 6, 3 ] ];;
gap> M:=[ [ 3, 4, 5 ], [ 6, 1, 2 ] ];;
gap> N:=[ [ 2, 3, 5 ], [ 6, 4, 1 ] ];;
gap> P:=[ [ 1, 2, 5 ], [ 6, 3, 4 ] ];;
gap> Y:=PoincareOctahedronCWComplex(L,M,N,P);;
gap> IsClosedManifold(Y);
true
gap> G:=FundamentalGroup(Y);;
gap> StructureDescription(G);
"SL(2,3)"
gap> L:=[[1,2,3,4,5,6],[8,9,10,11,12,7]];;
gap> M:=[[1,7,8,2],[11,10,4,5]];;
gap> N:=[[2,8,9,3],[12,11,5,6]];;
gap> P:=[[3,9,10,4],[7,12,6,1]];;
gap> Y:=PoincarePrismCWComplex(L,M,N,P);;
gap> IsClosedManifold(Y);
true
gap> StructureDescription(FundamentalGroup(Y));
"C3 : C4"
gap> STOP_TEST( "tst.tst", 1000 );
