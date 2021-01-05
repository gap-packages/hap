Read("homologyOperations.gd");


#####################################################################
#####################################################################
##
##  Method for applying TensorWithGModule(_,A) to one boundary homomorphism 
##  in a free ZG-resolution R. We won't document this method in the manual
##  as it is really only for internal use.
InstallMethod(TensorWithGModule,
[IsHapResolution,IsInt,IsGOuterGroup,IsGOuterGroup,IsGOuterGroup],
function(R,n,A,M,N)
local phitmp,phi,UM,UN, act, MAT,i,x,b,fn,s;

# R is a free ZG-resolution and we are considering the boundary 
# homomorphism d_n:R_n ---> R_{n-1} where n>=0.
# A is a ZG-module represented as a G-Outer group.
# M:=DirectProductGog(List([1..R!.dimension(n)],i->A));
# N:=DirectProductGog(List([1..R!.dimension(n-1)],i->A));
UM:=ActedGroup(M);  #M=R_n ----> N=R_{n-1} is represented by MAT
UN:=ActedGroup(N);  #UM is the group underlying M=AxAx...xA
act:=OuterAction(A);

MAT:=List([1..R!.dimension(n-1)],i->List([1..R!.dimension(n)],j->[]));

for i in [1..R!.dimension(n)] do
   for x in R!.boundary(n,i) do
      Add(MAT[AbsInt(x[1])][i],SignInt(x[1])*x[2]);
   od;
od;

	##########################################################
	fn:=function(j,a) #Here a is in A and i is the component of 
                          #the product
	local x,i,ans;  
	
	ans:=One(UN);
	for i in [1..R!.dimension(n-1)] do
	   for x in MAT[i][j] do
	      ans:=ans*   Image(Embedding(UN,i),
	      act(R!.elts[AbsInt(x)],a)^SignInt(x) );
	   od;
	od;
	return ans;
	end;
	##########################################################

phitmp:=GroupHomomorphismByFunction(UM,UN,x->
Product( List([1..R!.dimension(n)],j->fn(j,Image(Projection(UM,j),x))) ));

phi:=GroupHomomorphismByImagesNC(UM,UN,   
GeneratorsOfGroup(UM),
List(GeneratorsOfGroup(UM),a->Image(phitmp,a)));

return GOuterGroupHomomorphism(M,N,phi);
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
##
##  Method for applying TensorWithG(_,A) to a free ZG-resolution R.
InstallOtherMethod(TensorWithGModule,
[IsHapResolution,IsGOuterGroup],
function(R,A)
local gogs,homos,Boundary, properties, L;

L:=EvaluateProperty(R,"length");
homos:=[1..L];  #homos[1] has source in degree 1.
gogs:=[1..L+1]; #gogs[1] is the chain module C_0


########################################
Boundary:=function(n)
if IsInt(homos[n]) then
   if IsInt(gogs[n]) then gogs[n]:=
      DirectProductGog(List([1..R!.dimension(n-1)],i->A)); 
   fi;
   if IsInt(gogs[n+1]) then gogs[n+1]:=
      DirectProductGog(List([1..R!.dimension(n)],i->A)); 
   fi;
   homos[n]:=TensorWithGModule(R,n,A,gogs[n+1],gogs[n]);
fi;
return homos[n];
end;
########################################

return Objectify(HapGComplex,
                rec(
		resolution:=R,
		module:=A,
                boundary:=Boundary,
                properties:=
                [["length",L],
                ["type", "Gcomplex"],
                 ]));
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
##
##  Method for homology of a G-Complex
InstallMethod(HomologyModule,
[IsHapGComplex,IsInt],
function(C,n)
local  H, N, nat, hom; 

if n=0 then
   H:=Target(C!.boundary(n+1));
else
   H:=Kernel(C!.boundary(n));
fi;

N:=Range(C!.boundary(n+1));
nat:=NaturalHomomorphism(H,N);
hom:=Target(nat);
hom!.nat:=nat;
hom!.complex:=C;

return hom; 
end);
#####################################################################
#####################################################################
 
