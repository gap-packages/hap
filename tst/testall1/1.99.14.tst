gap> START_TEST("HAP library");
gap> A:=1;;C:=2;;D:=3;;B:=4;;
gap> Ap:=5;;Cp:=6;;Dp:=7;;Bp:=8;;
gap> L:=[[A,B,D,C],[Bp,Dp,Cp,Ap]];;
gap> M:=[[A,B,Bp,Ap],[Cp,C,D,Dp]];;
gap> N:=[[A,C,Cp,Ap],[D,Dp,Bp,B]];;
gap> Ex3:=PoincareCubeCWComplex(L,M,N);
Regular CW-complex of dimension 3

gap> STOP_TEST( "tst.tst", 1000 );
