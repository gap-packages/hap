gap> G:=DihedralGroup(64);; N:=Center(G);;
gap> R:=ResolutionNormalSeries([G,N],6);;
gap> C:=FilteredTensorWithIntegersModP(R,2);;
gap> P:=PersistentHomologyOfFilteredChainComplex(C,5,2);;
gap> BarCodeDisplay(P);;
