gap> START_TEST("HAP library");
gap> G:=Image(IsomorphismFpGroup(SmallGroup(24,3)));;
gap> CC:=SmallCatOneGroup(24,3,1);;
gap> CrossedInvariant(G,CC);
33
gap> R:=ResolutionFiniteGroup(SymmetricGroup(3),3);;
gap> IdentityAmongRelators(R,4)[2];;
gap> STOP_TEST( "tst.tst", 1000 );
