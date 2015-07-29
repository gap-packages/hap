gap> START_TEST("HAP library");
gap> psi1:= [[1,0,1], [0,-1,0], [0,0,1]];;
gap> psi2:= [[1,0,0], [0,1,1], [0,0,1]];;
gap> G:=AffineCrystGroupOnLeft([psi1,psi2]);;
gap> G:=StandardAffineCrystGroup(G);;
gap> Y:=EquivariantEuclideanSpace(G,[1/6,1/3]);
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
<fp group of size infinity on the generators [ f1, f2 ]>
gap> RelatorsOfFpGroup(F);
[ f1^-1*f2*f1^-1*f2^-1 ]
gap> STOP_TEST( "tst.tst", 1000 );


