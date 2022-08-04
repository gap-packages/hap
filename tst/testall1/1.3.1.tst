#D
gap> START_TEST("HAP library");
gap> A:=[[0,11,10,14,22],
>        [11,0,3,13,21],
>        [10,3,0,12,20],
>        [14,13,12,0,16],
>        [22,21,20,16,0]];;
gap> G:=SymmetricMatrixToGraph(A,10);;
gap> p:=PiZero(G);;
gap> Classify([1..5],p[2]);
[ [ 1, 2, 3 ], [ 4 ], [ 5 ] ]
gap> #DisplayDendrogramMat(A,10,2);
gap> STOP_TEST( "tst.tst", 1000 );

