gap> START_TEST("HAP library");
gap> M:=PureCubicalKnot(3,1);;
gap> N:=PureCubicalComplexToCubicalComplex(M);;
gap> Size(N);
1656
gap> ContractCubicalComplex(N);
gap> Size(N);
184
gap> N:=PureCubicalComplexToCubicalComplex(M);;
gap> D:=DVFReducedCubicalComplex(N);;
gap> NamesOfComponents(D);
[ "vectors", "rewrite", "properties", "binaryArray" ]
gap> STOP_TEST( "tst.tst", 1000 );
