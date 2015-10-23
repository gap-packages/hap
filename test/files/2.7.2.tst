gap> START_TEST("HAP library");
gap> 2simplices:= [[1,2,4], [2,4,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9], [4,5,8], [5,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7], [1,5,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
gap> T:=RegularCWComplex(SimplicialComplex(2simplices));;
gap> D:=DiagonalApproximation(T);;
gap> p:=D.projection;;
gap> P:=ChainMap(p);;
gap> IsIsomorphismOfAbelianFpGroups(Homology(P,0));
true
gap> IsIsomorphismOfAbelianFpGroups(Homology(P,1));
true
gap> IsIsomorphismOfAbelianFpGroups(Homology(P,2));
true
gap> Size(Source(p));
180
gap> Size(Source(P));
4
gap> Size(Target(p));
54
gap> Size(Target(P));
4
gap> del:=D.inclusion;;
gap> Del:=ChainMap(del);;
gap> Size(Target(Del));
16
gap> STOP_TEST( "tst.tst", 1000 );


