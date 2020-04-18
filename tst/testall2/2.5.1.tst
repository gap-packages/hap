gap> START_TEST("HAP library");
gap> SetAssertionLevel(0);;
gap> gens:= [ [ [-1,0,0,0], [0,-1,0,0], [0,0,1,1/2], [0,0,0,1] ], [ [0,-1,0,0], [1,-1,0,0], [0,0,1,1/3], [0,0,0,1] ], [ [1,0,0,1], [0,1,0,0], [0,0,1,0], [0,0,0,1] ] ];;
gap> G:=AffineCrystGroupOnLeft(gens);;
gap> G:=StandardAffineCrystGroup(G);;
gap> Y:=EquivariantEuclideanSpace(G,[0,0,0]);
Equivariant CW-complex of dimension 3

gap> Y!.dimension(0);
2
gap> Y!.dimension(1);
5
gap> Y!.dimension(2);
4
gap> Y!.dimension(3);
1
gap> C:=ChainComplexOfQuotient(Y);;
gap> Homology(C,0);
[ 0 ]
gap> Homology(C,1);
[ 0 ]
gap> Homology(C,2);
[ 0 ]
gap> Homology(C,3);
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


