#D
gap> START_TEST("HAP library");
gap> n:=1000;; Bettis:=[];;
gap> for p in [1..20]*(1/40000) do
> K:=RandomSimplicialGraph(n,p);;
> B:=BettiNumber(K,0);;
> Add(Bettis,[p,B]);
> od;
gap> #ScatterPlot(Bettis);
gap> STOP_TEST( "tst.tst", 1000 );


