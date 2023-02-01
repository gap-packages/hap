gap> START_TEST("HAP library");
gap> K:=QuillenComplex(SymmetricGroup(6),2);
Simplicial complex of dimension 2.

gap> Size(K);
1815
gap> Y:=RegularCWComplex(K);;
gap> C:=CriticalCells(Y);;
gap> n:=0;;Length(Filtered(C,c->c[1]=n));
1
gap> n:=1;;Length(Filtered(C,c->c[1]=n));
16
gap> n:=2;;Length(Filtered(C,c->c[1]=n));
0
gap> STOP_TEST( "tst.tst", 1000 );


