#
gap> START_TEST("HAP library");
gap> children:=List([1..2],i->ChildProcess());;    
gap> fn:=function(i);
> return BogomolovMultiplier(SmallGroup(3^5,i));
> end;;
gap> for s in children do
>  ChildPut(fn,"fn",s);
> od;
gap> L:=ParallelList([1..5],"fn",children); 
[ [  ], [  ], [  ], [  ], [  ] ]
gap> STOP_TEST( "tst.tst", 1000 );


