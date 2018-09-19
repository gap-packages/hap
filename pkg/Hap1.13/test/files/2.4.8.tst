gap> START_TEST("HAP library");
gap> Read(Concatenation(dir,"data248.txt"));
gap> F:=ThickeningFiltration(Y,20);;
gap> R:=FiltrationTerm(F,19);;
gap> for i in Reversed([1..18]) do
> S:=FiltrationTerm(F,i);
> R:=ContractedComplex(R,S);
> od;
gap> K:=RegularCWComplex(Nerve(R));;
gap> M:=ContractedComplex(K);;
gap> Display(Graph(M));
gap> STOP_TEST( "tst.tst", 1000 );


