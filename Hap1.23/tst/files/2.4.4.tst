gap> START_TEST("HAP library");
gap> Read(Concatenation(dir,"data244.txt"));;
gap> A:=NullMat(350,350);;
gap> for x in S do A[x[1]][x[2]]:=1; od;
gap> Y:=PureCubicalComplex(A);;
gap> F:=ThickeningFiltration(Y,100,4);
Filtered pure cubical complex of dimension 2.

gap> P:=PersistentBettiNumbers(F,1);;
gap> BarCodeCompactDisplay(P);
gap> STOP_TEST( "tst.tst", 1000 );


