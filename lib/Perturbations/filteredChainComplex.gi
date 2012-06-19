#G:=SmallGroup(64,134);
#L:=LowerCentralSeries(G);
#R:=ResolutionNormalSeries(L,5);
#C:=FilteredTensorWithIntegers(R);


#################################################################
#################################################################
InstallGlobalFunction(CyclesOfFilteredChainComplex,
function(C,r,p,q)
local
	SOURCE,TARGET, BasisKer, FD, FL, S,T;
					#p column
					#q row

FL:=EvaluateProperty(C,"filtration_length");
###############################
FD:=function(r,k);
if r<0 then return  0; fi;
return 
C!.filteredDimension(Minimum(r,FL),k);
end;
###############################

SOURCE:=[FD(p-1,q)+1..FD(p,q)];
TARGET:=List(SOURCE, x-> C!.boundary(p+q,x) );
S:=FD(p-r,q-1)+1;
T:=Length(TARGET[1]);
T:=[S..T];
TARGET:=List(TARGET, x-> x{T});
ConvertToMatrixRep(TARGET);
BasisKer:=LLLReducedBasis(TARGET,"linearcomb").relations;


return BasisKer;
end);
#################################################################
#################################################################


#################################################################
#################################################################
InstallGlobalFunction(BoundariesOfFilteredChainComplex,
function(C,r,p,q)
local
         BasisKer, Bds, b,x;
                                        #p column
                                        #q row
BasisKer:=CyclesOfFilteredChainComplex(C,r-1,p+r-1,q-r+2);

Bds:=[];
for b in BasisKer do
x:=C!.boundary(p+1+1);
od;
end);
#################################################################
#################################################################


