#D
gap> START_TEST("HAP library");
gap> D:=[];; n:=30;;
gap> for p in [1..20]*(1/500) do
> Y:=RandomSimplicialTwoComplex(n,p);;
> G:=FundamentalGroup(Y);
> Add(D,[p,G]);
> od;
gap> P:=[];;
gap> for i in [1..Length(D)] do
> if Length(RelatorsOfFpGroup(D[i][2]))=0 then
> Add(P,[D[i][1],Length(GeneratorsOfGroup(D[i][2]))]);
> else
> Add(P,[D[i][1],Length(GeneratorsOfGroup(D[i][2])), "red"]);
> fi;
> od;
gap> #ScatterPlot(P);
gap> STOP_TEST( "tst.tst", 1000 );


