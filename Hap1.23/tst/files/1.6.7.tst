gap> START_TEST("HAP library");
gap> K:=PureCubicalKnot(3,1);
prime knot 1 with 3 crossings

gap> G:=WirtingerGroup(K);
<fp group on the generators [ f1, f2, f3 ]>
gap> S:=SimplifiedFpGroup(G);
<fp group of size infinity on the generators [ f2, f3 ]>
gap> RelatorsOfFpGroup(S);
[ f2*f3*f2^-1*f3^-1*f2^-1*f3 ]
gap> STOP_TEST( "tst.tst", 1000 );


