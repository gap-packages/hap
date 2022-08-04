#
gap> START_TEST("HAP library");
gap> x:=(1,2)(5,6)(7,8)(11,12);; y:=(2,3)(4,5)(8,9)(10,11);;
gap> z:=(3,4)(5,7)(6,8)(9,10);; G:=Group(x,y,z);;
gap> #CayleyGraphOfGroupDisplay(G,[x,y,z]);
gap> Y:=EquivariantTwoComplex(G);
Equivariant CW-complex of dimension 2

gap> F:=FundamentalGroupOfQuotient(Y);
<fp group on the generators [ f1, f2, f3 ]>
gap> RelatorsOfFpGroup(F);
[ f1^2, f2^2, f3^2, f1^-1*f3^-1*f1*f3, f2^-1*f3^-1*f2^-1*f3*f2*f3, 
  (f1^-1*f2^-1)^2*f1^-1*(f2*f1)^2*f2 ]
gap> H:=Group(x*y,x*z,y*z);;
gap> W:=RestrictedEquivariantCWComplex(Y,H);
Equivariant CW-complex of dimension 2

gap> FH:=FundamentalGroupOfQuotient(W);
<fp group on the generators [ f1, f2, f3, f4, f5 ]>
gap> RelatorsOfFpGroup(FH);
[ f1, f1, f2*f3, f3*f2, f4*f5, f5*f4, f5^-1*f1*f4, f1^-1*f4^-1*f5, 
  f2^-1*f5^-1*f2^-1*f4*f3*f4, f3^-1*f4^-1*f3^-1*f5*f2*f5, f3^-2*(f2*f1)^2*f2, 
  (f1^-1*f2^-1)^2*f1^-1*f3^3 ]
gap> xz:=(1,2)(3,4)(5,8)(6,7)(9,10)(11,12);;
gap> yz:=(2,4,7,5,3)(6,8,10,11,9);;
gap> H:=Group(xz, yz);;
gap> W:=EquivariantTwoComplex(H);
Equivariant CW-complex of dimension 2

gap> FH:=FundamentalGroupOfQuotient(W);
<fp group on the generators [ f1, f2 ]>
gap> RelatorsOfFpGroup(FH);
[ f2^2, f1^5, (f1*f2)^2*f1*f2^-1 ]
gap> STOP_TEST( "tst.tst", 1000 );


