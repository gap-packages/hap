gap> START_TEST("HAP library");
gap> torsion:=function(n,p)
> local H, Y;
> Y:=RegularCWComplex(RandomSimplicialTwoComplex(n,p));
> H:=Homology(Y,1);
> H:=Filtered(H,x->not x=0);
> return H;
> end;
function( n, p ) ... end
gap> n:=50000;;
gap> t:=torsion(75,n/2000000);
[  ]
gap> STOP_TEST( "tst.tst", 1000 );
