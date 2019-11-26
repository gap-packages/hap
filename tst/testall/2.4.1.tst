#D
gap> START_TEST("HAP library");
gap> n:=100;; Betti1:=[];; Betti2:=[];;
gap> for p in [0..10]*(1/200) do
> G:=RandomSimplicialGraph(n,p);;
> K:=CliqueComplex(G,3);;
> Add(Betti1,[p,BettiNumber(K,1)]);
> Add(Betti2,[p,BettiNumber(K,2)]);
> od;
gap> #ScatterPlot(Betti1);
gap> #ScatterPlot(Betti2);
gap> STOP_TEST( "tst.tst", 1000 );


