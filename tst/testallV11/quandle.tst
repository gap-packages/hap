gap> START_TEST("HAP library");
gap> Q:=Quandle(5,21);;
gap> IsConnected(Q);
true
gap> IsLatinQuandle(Q);
true
gap> G:=DihedralGroup(64);;
gap> Q:=ConjugationQuandle(G,1);;
gap> Size(Q);
64
gap> IsConnected(Q);
false
gap> Q:=[1..9];; e:=2;; G:=TransitiveGroup(9,15);; mu:=(1,8,7,4,9,5,3,6);;
gap> IsQuandleEnvelope(Q,G,e,mu); 
true
gap> QE:=QuandleQuandleEnvelope(Q,G,e,mu);;
gap> IsQuandle(QE); 
true
gap> IsConnected(QE);
true
gap> K:=PureCubicalKnot(3,1);;
gap> G:=GaussCodeOfPureCubicalKnot(K);;
gap> P:=PresentationKnotQuandle(G);
Quandle presentation of 3 generators and 3 relators.

gap> P!.generators;
[ 1 .. 3 ]
gap> P!.relators;
[ [ [ 3, 2 ], 1 ], [ [ 1, 3 ], 2 ], [ [ 2, 1 ], 3 ] ]
gap> PD:=PlanarDiagramKnot(3,1);
[ [ 1, 4, 2, 5 ], [ 3, 6, 4, 1 ], [ 5, 2, 6, 3 ] ]
gap> G:=PD2GC(PD);
[ [ [ -1, 3, -2, 1, -3, 2 ] ], [ -1, -1, -1 ] ]
gap> QK:=PresentationKnotQuandleKnot(12,1000);
Quandle presentation of 12 generators and 12 relators.

gap> Q:=ConnectedQuandle(30,2);;
gap> NumberOfHomomorphisms(QK,Q);
1230
gap> Q:=ConjugationQuandle(SymmetricGroup(5),1);;
gap> P:=PathComponents(Q);;
gap> P:=List(P,x->AsMagma(x));;
gap> Q:=P[5];;
gap> F:=FundamentalGroup(Q);;
gap> IdGroup(F);
[ 4, 1 ]
gap> STOP_TEST( "tst.tst", 1000 );
