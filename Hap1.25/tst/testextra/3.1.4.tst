gap> START_TEST("HAP library");
gap> SetAssertionLevel(0);;
gap> SetCrystGroupDefaultAction( RightAction );
gap> G:=SpaceGroup(4,122);;
gap> gens:=GeneratorsOfGroup(G);;
gap> T:=AffineCrystGroupOnRight(gens{[3,4,5,6]});;
gap> f:=GroupHomomorphismByFunction(T,G,x->x);;
gap> RG:=ResolutionCubicalCrystGroup(G,5);;
gap> RT:=ResolutionCubicalCrystGroup(T,5);;
gap> Rf:=EquivariantChainMap(RT,RG,f);;
gap> F:=TensorWithIntegers(Rf);;
gap> h:=Homology(F,2);;
gap> AbelianInvariants(Kernel(h));
[ 0, 0, 0, 0, 0 ]
gap> AbelianInvariants(Image(h));
[ 0 ]
gap> STOP_TEST( "tst.tst", 1000 );


