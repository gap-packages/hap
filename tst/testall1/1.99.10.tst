gap> START_TEST("HAP library");
gap> K:=ComplexProjectiveSpace(2);;
gap> M:=ConnectedSum(K,K,-1);;
gap> Display(IntersectionForm(M));;
[ [   1,   0 ],
  [   0,  -1 ] ]
gap> STOP_TEST( "tst.tst", 1000 );
