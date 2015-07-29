gap> START_TEST("HAP library");
gap> M:=ReadImageAsPureCubicalComplex(Concatenation(dir,"image1.3.3.eps"),300);;

gap> while BettiNumber(PureComplexComplement(M),0)>1 do
> M:=ThickenedPureComplex(M);
> od;

gap> F:=ConcentricFiltration(ComplementOfPureComplex(M),20);
Filtered Pure Cubical Complex of Dimension 2.

gap> P:=PersistentBettiNumbers(F,0);;
gap> BarCodeCompactDisplay(P);;
gap> STOP_TEST( "tst.tst", 1000 );


