#E
gap> START_TEST("HAP library");
gap> SetAssertionLevel( 1 );
gap> psi1:= [[1,0,1], [0,-1,0], [0,0,1]];;
gap> psi2:= [[1,0,0], [0,1,1], [0,0,1]];;
gap> SetCrystGroupDefaultAction(LeftAction);;
gap> G:=AffineCrystGroupOnLeft([psi1,psi2]);;
gap> GG:=StandardAffineCrystGroup(G);;
gap> Y:=EquivariantEuclideanSpace(GG,[1/6,1/3]);
Equivariant CW-complex of dimension 2

gap> NumbersOfCellOrbits:=List([0..2],Y!.dimension);
[ 2, 3, 1 ]
gap> VertexStabilizerOrders:= List([1,2],i->Order(Y!.stabilizer(0,i)));
[ 1, 1 ]
gap> EdgeStabilizerOrders:= List([1,2,3],i->Order(Y!.stabilizer(1,i)));
[ 1, 1, 1 ]
gap> FaceStabilizerOrders:=Order(Y!.stabilizer(2,1));
1
gap> Size(Y!.boundary(2,1));
6
gap> F:=FundamentalGroupOfQuotient(Y);
<fp group of size infinity on the generators [ x, y ]>
gap> RelatorsOfFpGroup(F);
[ x^-1*y*x^-1*y^-1 ]
gap> STOP_TEST( "tst.tst", 1000 );
