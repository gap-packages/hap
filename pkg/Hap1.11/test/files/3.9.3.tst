gap> G:=SylowSubgroup(MathieuGroup(12),2);;
gap> A:=LHSSpectralSequenceLastSheet(G);;
gap> P:=HilbertPoincareSeries(A);;
gap> RankHomologyPGroup(G,P,1001);
251503
