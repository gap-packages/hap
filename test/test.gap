#(C) Graham Ellis, 2005-2006

#####################################################################
TestHap:=function()
local H0,H1,H2,H3,H4,H5,H6,R,TR,Tensor,Bool;

H0:=GroupHomology(AlternatingGroup(5),2,2);;

H1:=GroupHomology(SmallGroup(32,3),4,2);;

H2:=GroupHomology(AbelianGroup([2,4,6]),4);;

H3:=GroupHomology([[1,[2,3]],[2,[3,4]]],2);;

R:=ResolutionNilpotentGroup(SmallGroup(64,135),5);;
TR:=TensorWithIntegers(R);;
H4:=Homology(TR,4);;

Tensor:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3));;

P:=PresentationOfResolution(ResolutionAbelianGroup([0,2,4],3)).relators;
P:=Length(P);

R:=ResolutionAbelianGroup([0,0,0],5);;
H6:=Homology(TensorWithIntegers(R),3);

if LoadPackage("Polycyclic")=true then 
R:=ResolutionNilpotentGroup(HeisenbergPcpGroup(3),3);
H5:=Homology(TensorWithIntegers(R),2);
H5:=Length(H5);
else
H5:=14;
fi;


Bool:=
H0=[2]
and
H1=5
and
H2=[ 2, 2, 2, 2, 2, 2, 2, 2 ]
and
H3=[ 0, 0 ]
and
H4=[ 2, 2, 2, 2, 2, 2, 2 ]
and
Tensor=[ 2 ]
and
P=5
and
H5=14
and H6[1]=0
and
IsSuperperfect(PerfectGroup(120,1));

if Bool then
Print("\n\n HAP seems to be working fine. \n");
else
Print("\n\n There are some problems with HAP. \n");
fi;
end;
#####################################################################
