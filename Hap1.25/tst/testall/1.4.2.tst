gap> START_TEST("HAP library");
gap> dir:=Filename(DirectoriesPackageLibrary("HAP","tst/testall")[1],"bing.txt");;
gap> Read(dir);
gap> BingsHouse;
Regular CW-complex of dimension 2

gap> BingsModifiedHouse;
Regular CW-complex of dimension 3

gap> Size(BingsHouse);
309
gap> Size(BingsModifiedHouse);
311
gap> Size(ContractedComplex(BingsModifiedHouse));
1
gap> Size(ContractedComplex(BingsHouse));
309
gap> STOP_TEST( "tst.tst", 1000 );


