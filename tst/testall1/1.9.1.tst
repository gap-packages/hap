#
gap> START_TEST("HAP library");
gap> x:=(1,2)(5,6)(7,8)(11,12);; y:=(2,3)(4,5)(8,9)(10,11);;
gap> z:=(3,4)(5,7)(6,8)(9,10);; G:=Group(x,y,z);;
gap> #CayleyGraphOfGroupDisplay(G,[x,y,z]);
gap> Y:=EquivariantTwoComplex(G);;
gap> F:=FundamentalGroupOfQuotient(Y);;
gap> Order(F);
120
gap> H:=Group(x*y,x*z,y*z);;
gap> W:=RestrictedEquivariantCWComplex(Y,H);;
gap> FH:=FundamentalGroupOfQuotient(W);;
gap> Order(FH);
60
gap> xz:=(1,2)(3,4)(5,8)(6,7)(9,10)(11,12);;
gap> yz:=(2,4,7,5,3)(6,8,10,11,9);;
gap> H:=Group(xz, yz);;
gap> W:=EquivariantTwoComplex(H);;
gap> FH:=FundamentalGroupOfQuotient(W);;
gap> Order(FH);
60
gap> STOP_TEST( "tst.tst", 1000 );


