gap> START_TEST("HAP library");
gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> P:=PoincareSeriesLHS(G);
(1)/(-x_1^3+3*x_1^2-3*x_1+1)
gap> RankHomologyPGroup(G,P,1000);
251000
gap> STOP_TEST( "tst.tst", 1000 );
