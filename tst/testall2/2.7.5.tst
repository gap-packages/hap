gap> START_TEST("HAP library");
gap> SetAssertionLevel( 0);;
gap> gens:= [ [ [-1,0,0,0], [0,-1,0,0], [0,0,1,1/2], [0,0,0,1] ], [ [0,-1,0,0], [1,-1,0,0], [0,0,1,1/3], [0,0,0,1] ], [ [1,0,0,1], [0,1,0,0], [0,0,1,0], [0,0,0,1] ] ];;
gap> G:=AffineCrystGroupOnLeft(gens);;
gap> G:=StandardAffineCrystGroup(G);;
gap> G:=Image(IsomorphismPcpGroup(G));;
gap> E_G:=Resolution(G,4);;
gap> C:=HomToIntegers(E_G);;
gap> Cohomology(C,1);
[ 0 ]
gap> Cohomology(C,2);
[ 0 ]
gap> Cohomology(C,3);
[ 0 ]
gap> CupProduct(E_G,1,2,[1],[1]);
[ -1 ]
gap> STOP_TEST( "tst.tst", 1000 );


