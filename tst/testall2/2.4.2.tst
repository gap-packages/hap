gap> START_TEST("HAP library");
gap> n:=15;; Betti2:=[];; Betti3:=[];;
gap> for p in [0..10]*(1/200) do
> Y:=RandomSimplicialTwoComplex(n,p);;
> K:=CliqueComplex(Y,4);;
> W:=RegularCWComplex(K);;
> W:=ContractedComplex(W);;
> b2:=BettiNumber(W,2);;
> Add(Betti2,[p,b2]);;
> b3:=BettiNumber(W,3);;
> Add(Betti3,[p,b3]);;
> od;
gap> #ScatterPlot(Betti2);
gap> #ScatterPlot(Betti3);
gap> STOP_TEST( "tst.tst", 1000 );


