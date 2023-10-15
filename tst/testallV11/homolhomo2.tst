gap> START_TEST("HAP library");
gap> G:=SymmetricGroup(4);;H:=AlternatingGroup(4);;
gap> f:=GroupHomomorphismByFunction(H,G,x->x);;
gap> p:=2;; deg:=3;;
gap> F:=ModPCohomologyRing(f,p,deg);
[ v.1, v.2, v.4+v.6, v.5 ] -> [ v.1, 0*v.1, v.4+v.5+v.6, 0*v.1 ]
gap> STOP_TEST( "tst.tst", 1000 );


