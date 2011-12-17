H1:=GroupHomology(SmallGroup(32,3),4);;

H2:=GroupHomology(AbelianGroup([2,4,6]),4);;

H3:=GroupHomology([[1,[2,3]],[2,[3,4]]],2);;

R:=ResolutionNilpotentGroup(SmallGroup(64,135),5);;
TR:=TensorWithIntegers(R);;
H4:=Homology(TR,4);;
Tensor:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3));;

Bool:=
H1=[4,4]
and
H2=[ 2, 2, 2, 2, 2, 2, 2, 2 ]
and
H3=[ 0, 0 ]
and
H4=[ 2, 2, 2, 2, 2, 2, 2 ]
and
Tensor=[ 2 ];;

if Bool then
Print("HAP seems to be working fine. \n");
else
Print("There are some problems with HAP. \n");
fi;
