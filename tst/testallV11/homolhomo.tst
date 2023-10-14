gap> START_TEST("HAP library");
gap> G:=SymmetricGroup(5);;                    
gap> H:=AlternatingGroup(4);;
gap> f:=GroupHomomorphismByFunction(H,G,x->x);;
gap> F:=GroupHomology(f,3);;
gap> AbelianInvariants(Range(F)/Image(F));
[ 2, 2 ]
gap> Fmod2:=GroupHomology(f,3,2);;
gap> AbelianInvariants(Range(Fmod2)/Image(Fmod2));
[ 2, 2 ]
gap> STOP_TEST( "tst.tst", 1000 );


