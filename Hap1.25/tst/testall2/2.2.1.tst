gap> START_TEST("HAP library");
gap> F:=FreeGroup(3);;x:=F.1;;y:=F.2;;z:=F.3;;
gap> rels:=[x*y*x*(y*x*y)^-1, y*z*y*(z*y*z)^-1, z*x*z*(x*z*x)^-1];;
gap> IsAspherical(F,rels);
Test inconclusive.

fail
gap> F:=FreeGroup(6);;
gap> x:=F.1;;y:=F.2;;z:=F.3;;a:=F.4;;b:=F.5;;c:=F.6;;
gap> rels:=[a^-1*x*y, b^-1*y*z, c^-1*z*x, a*x*(y*a)^-1, b*y*(z*b)^-1, c*z*(x*c)^-1];;
gap> IsAspherical(F,rels);
Presentation is aspherical.

true
gap> STOP_TEST( "tst.tst", 1000 );


