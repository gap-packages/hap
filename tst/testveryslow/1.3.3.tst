#ED
gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/examples"),"image1.3.3.png");;
gap> t:=ReadImageAsPureCubicalComplex(dir,"matrix");;
gap> t:=Int((3/10)*Maximum(Flat(t)));;
gap> M:=ReadImageAsPureCubicalComplex(dir,t);;
gap> while BettiNumber(PureComplexComplement(M),0)>1 do
> M:=ThickenedPureComplex(M);
> od;;
gap> F:=ConcentricFiltration(ComplementOfPureComplex(M),10);
Filtered pure cubical complex of dimension 2.

gap> P:=PersistentBettiNumbers(F,0);
[ [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ], 
  [ 0, 0, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1 ], 
  [ 0, 0, 0, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 1 ], 
  [ 0, 0, 0, 0, 5, 5, 5, 4, 4, 4, 4, 4, 4, 1 ], 
  [ 0, 0, 0, 0, 0, 5, 5, 4, 4, 4, 4, 4, 4, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 5, 4, 4, 4, 4, 4, 4, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 1 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ] ]
gap> #BarCodeCompactDisplay(P);;
gap> STOP_TEST( "tst.tst", 1000 );


