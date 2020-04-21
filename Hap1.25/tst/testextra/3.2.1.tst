#E
gap> START_TEST("HAP library");
gap> A:=[[0,-1,0,0,],[1,0,0,0,],[0,0,0,1],[0,0,-1,0]];;
gap> B:=[[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]];;
gap> Q:=Group([A,B]);;
gap> R:=PolytopalComplex(Q,[1,0,0,0]);;
gap> R:=FreeGResolution(R,4);;
gap> List([0..4],i->R!.dimension(i));
[ 1, 3, 4, 2, 1 ]
gap> A:=GroupHomomorphismByFunction(Q,Q,x->x);;
gap> C:=HomToIntegralModule(R,A);;
gap> Cohomology(C,3);
[ 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


