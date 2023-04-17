#E
gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/data246.txt");;
gap> F:=ThickeningFiltration(T,20);;
gap> Y:=FiltrationTerm(F,12);;
gap> Size(Y);
1008893
gap> Y:=ZigZagContractedComplex(Y,2);;
gap> Y:=ContractedComplex(RegularCWComplex(Y));;
gap> Cohomology(Y,0);
[ 0 ]
gap> Cohomology(Y,1);
[ 0, 0 ]
gap> Cohomology(Y,2);
[ 0 ]
gap> G:=FundamentalGroup(Y);;
gap> IsAbelian(G);
true
gap> G:=NilpotentQuotient(G,1);;
gap> E_G:=Resolution(G,3);;
gap> CupProduct(E_G,1,1,[1,0],[0,1]);
[ -1 ]
gap> STOP_TEST( "tst.tst", 1000 );


