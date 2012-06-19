#(C) Graham Ellis, 2005-2006

#####################################################################
TestHap:=function()
local P, H0,H1,H2,H3,H4,H4a,H5,H6,H7,H8,H9,R,TR,Tensor,Tensor2,
S5,S4,A,AS5,AS4,D,Bool,Poincare;

H0:=GroupHomology(AlternatingGroup(5),2,2);;

H1:=GroupHomology(SmallGroup(32,3),4,2);;

H2:=GroupHomology(AbelianGroup([2,4,6]),4);;

H3:=GroupHomology([[1,[2,3]],[2,[3,4]]],2);;

R:=ResolutionNilpotentGroup(SmallGroup(64,135),5);;
TR:=TensorWithIntegers(R);;
H4:=Homology(TR,4);;
H4a:=IntegralRingGenerators(R,4);;

Tensor:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3));;
Tensor2:=ThirdHomotopyGroupOfSuspensionB(SymmetricGroup(3),12);;
P:=PresentationOfResolution(ResolutionAbelianGroup([0,2,4],3)).relators;
P:=Length(P);

R:=ResolutionAbelianGroup([0,0,0],5);;
H6:=Homology(TensorWithIntegers(R),3);
H7:=Cohomology(HomToIntegersModP(R,2),3);

if LoadPackage("Polycyclic")=true then 
R:=ResolutionNilpotentGroup(HeisenbergPcpGroup(3),3);
H5:=Homology(TensorWithIntegers(R),2);
H5:=Length(H5);
else
H5:=14;
fi;

S5:=SymmetricGroup(5);SetName(S5,"S5");
S4:=SymmetricGroup(4);SetName(S4,"S4");
A:=SymmetricGroup(3);SetName(A,"S3");
AS5:=GroupHomomorphismByFunction(A,S5,x->x);
AS4:=GroupHomomorphismByFunction(A,S4,x->x);
D:=[S5,S4,[AS5,AS4]];
R:=ResolutionGraphOfGroups(D,3);
H8:=Homology(TensorWithIntegers(R),2);
H9:=Homology(ChevalleyEilenbergComplex(SimpleLieAlgebra("A",2,Integers),3),2);
Poincare:=
 CoefficientsOfUnivariateRationalFunction(PoincareSeries(SmallGroup(64,134)));
Bool:=
H0=[2]
and
H1=[ 2, 2, 2, 2, 2 ]
and
H2=[ 2, 2, 2, 2, 2, 2, 2, 2 ]
and
H3=[ 0, 0 ]
and
H4=[ 2, 2, 2, 2, 2, 2, 2 ]
and
H4a=[ [ 0, 0, 1, 1, 0, -4, -4 ], [ 0, 0, 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 0, 0, 1 ]
 ]
and  
Tensor=[ 2 ]
and
Tensor2=[ 2 ]
and
P=5
and
H5=14
and H6[1]=0
and
H7=1
and
H8=[2,2]
and
H9=[ 3, 3, 3, 3, 3, 3 ]
and
IsSuperperfect(PerfectGroup(120,1))
and
Poincare=[ [ 1 ], [ 1, -3, 3, -1 ], 0 ];

if Bool then
Print("\n\n HAP seems to be working fine. \n");
else
Print("\n\n There are some problems with HAP. \n");
fi;
end;
#####################################################################
