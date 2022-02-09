gap> START_TEST("HAP library");
gap> F:=FreeGroup(2);; x:=F.1;; y:=F.2;;
gap> G:=F/[ x^2, x*y*x^-1*y^-2 ];;
gap> R:=ResolutionSmallGroup(G,8);;
gap> R!.dimension(4);
1
gap> R!.dimension(8);
1
gap> R!.boundary(4,1)=R!.boundary(8,1);
true
gap> STOP_TEST( "tst.tst", 1000 );


