#D
gap> START_TEST("HAP library");
gap> ReadPackage("HAP","tst/testall1/data248.txt");;
gap> F:=ThickeningFiltration(Y,7);;
gap> R:=FiltrationTerm(F,6);;
gap> for i in Reversed([1..5]) do
> S:=FiltrationTerm(F,i);
> R:=ContractedComplex(R,S);
> od;
gap> K:=RegularCWComplex(Nerve(R));;
gap> M:=ContractedComplex(K);;
gap> #Display(Graph(M));
gap> STOP_TEST( "tst.tst", 1000 );


