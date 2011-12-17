#####################################################################
#####################################################################
InstallOtherMethod(Dimension,
"Dimension of  simplicial complex",
[IsHapSimplicialComplex],
function(f) return EvaluateProperty(f,"dimension");
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(EulerCharacteristic,
"Euler characteristic  of a simplicial complex",
[IsHapSimplicialComplex],
function(M);
return
Sum(List([0..Dimension(M)],i->((-1)^i)*M!.nrSimplices(i)));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex],
function(M);
return
List([0..Dimension(M)],n->Homology(ChainComplex(M),n));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(M,n);
return
Homology(ChainComplex(M),n);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(Bettinumbers,
"Integral homologies  of a simplicial complex",
[IsHapSimplicialComplex],
function(M)
local H;
H:=Homology(M);
return List(H,x->Length(Filtered(x,a->a=0)));
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod(Homology,
"Homology  of a simplicial complex",
[IsHapSimplicialComplex,IsInt],
function(M,n);
return
Homology(ChainComplex(M),n);
end);
#####################################################################
#####################################################################



#####################################################################
#####################################################################
InstallOtherMethod(ChainComplex,
"Cellular chain complex of a simplicial complex",
[IsHapSimplicialComplex],
function(M)
return ChainComplexOfSimplicialComplex(M);;
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction(ChainComplexOfSimplicialComplex,
function(M)
local
        Dimsn,Boundary;

#############################
Dimsn:=function(n);
return M!.nrSimplices(n);
end;
#############################

#############################
Boundary:=function(n,i)
local V,Vhat, j, bnd;
V:=M!.simplices(n,i);

bnd:=List([1..M!.nrSimplices(n-1)],a->0);

for j in [1..Length(V)] do
Vhat:=StructuralCopy(V);
RemoveSet(Vhat,V[j]);
if IsOddInt(j) then
bnd[M!.enumeratedSimplex(Vhat)]:=1;
else
bnd[M!.enumeratedSimplex(Vhat)]:=-1;
fi;
od;

return bnd;
end;
#############################


return
Objectify(HapChainComplex,
           rec(
           dimension:=Dimsn,
           boundary:=Boundary,
           properties:=[
           ["length",EvaluateProperty(M,"dimension")],
           ["type","chainComplex"],
           ["characteristic",0]]
           ));

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(MaximalSimplicesToSimplicialComplex,
function(V)
local
	Vertices,
        NrSimplices,
	SimplicesLst,
	Simplices,
	EnumeratedSimplex,
	dim,
	n,v,x;

dim:=Maximum(List(V,v->Length(v)))-1;
Vertices:=Concatenation(V);
Vertices:=SSortedList(Vertices);
SimplicesLst:=List([1..dim+100],i->[]); #SLOPPY!!!
SimplicesLst[1]:=List(Vertices,x->[x]);

#####################
NrSimplices:=function(n);
return Length(SimplicesLst[n+1]);
end;
#####################

#####################  
Simplices:=function(n,i);
return SimplicesLst[n+1][i];
end;
#####################

#####################
EnumeratedSimplex:=function(v);
return Position(SimplicesLst[Length(v)],v);
end;
#####################

for v in V do
for n in [1..Length(v)] do
for x in Combinations(v,n) do
AddSet(SimplicesLst[n],x);
od;
od;
od;



return
Objectify(HapSimplicialComplex,
           rec(
	   vertices:=Vertices,
           nrSimplices:=NrSimplices,
           simplices:=Simplices,
	   enumeratedSimplex:=EnumeratedSimplex,
           properties:=[
           ["dimension",dim]]
           ));

end);
#################################################################
#################################################################

#################################################################
#################################################################
InstallGlobalFunction(SkeletonOfSimplicialComplex,
function(S,k)
local
	T,p,prp;
T:=Objectify(HapSimplicialComplex, rec());
T!.vertices:=StructuralCopy(S!.vertices);
T!.nrSimplices:=StructuralCopy(S!.nrSimplices);
T!.enumeratedSimplex:=StructuralCopy(S!.enumeratedSimplex);
prp:=StructuralCopy(S!.properties);
p:=PositionProperty(prp,x->"dimension" in x);
prp[p][2]:=k;
T!.properties:=prp;

return T;
end);
#################################################################
#################################################################

