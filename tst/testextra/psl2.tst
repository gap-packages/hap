gap> C:=ContractibleGcomplex("SL(2,O-2_a)");;
gap> PC:=QuotientOfContractibleGcomplex(C,Group(-One(C!.group)));;
gap> PR:=FreeGResolution(PC,7);;
gap> Homology(TensorWithIntegers(PR),6);
[ 2, 2, 2, 6 ]
