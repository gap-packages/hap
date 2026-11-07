gap> START_TEST("quotient-elements.tst");
gap> hapqinfo := InfoLevel(InfoWarning);;
gap> SetInfoLevel(InfoWarning,0);
gap> hapq := QuotientGroup((1,2,3,4,5,6),[(),(1,4)(2,5)(3,6)]);;
gap> Print(IsHapQuotientElement(InverseMutable(hapq)),"\n");
true
gap> Print(InverseMutable(hapq)=hapq^-1,"\n");
true
gap> Print(hapq*InverseMutable(hapq)=One(hapq),"\n");
true
gap> Print(Size(DirectProduct(Group(hapq),Group(hapq))),"\n");
9
gap> SetInfoLevel(InfoWarning,hapqinfo);
gap> Unbind(hapqinfo);
gap> Unbind(hapq);
gap> STOP_TEST("quotient-elements.tst", 1000);
