
#########################################################################
#########################################################################
InstallGlobalFunction(HAP_CupProductOfSimplicialComplex,
function(K)
local Y, C, CC, D, IntCoh, fns, gns, clgy, Iterate, n, k, cup, BigToSmall, 
      SmallToBig, crit, invcrit, Equiv, E1, E2 ;


Y:=RegularCWComplex(K);

Equiv:=ChainComplexEquivalenceOfRegularCWComplex(Y);
C:=Source(Equiv[2]);  ##This has fewer generators than cells of Y
CC:=Target(Equiv[2]);  ##This has generators equal to the cells of Y
D:=HomToIntegers(C);;
IntCoh:=IntegralCohomology("CohomologyAsFpGroup",true);;


E1:=HomToIntegers(Equiv[1]);
E2:=HomToIntegers(Equiv[2]);
###############################################
###############################################
#BigToSmall:=Equiv[1]!.mapping;
BigToSmall:=E2!.mapping;
#Inputs a vector v of length Y!.dimension(k) and
#returns a vector of length C!.dimension(k).
###############################################
###############################################

###############################################
###############################################
#SmallToBig:=Equiv[2]!.mapping;
SmallToBig:=E1!.mapping;
#Inputs a vector v of length C!.dimension(k) and
#returns a vector of length Y!.dimension(k).
###############################################
###############################################



fns:=[];
gns:=[];
clgy:=[];
###############################################
###############################################
Iterate:=function(n)
local HD, F, FhomG, G, gensG, preimgensG, f, g;
if Length(Cohomology(Y,n))>0 then

HD:=IntCoh(D,n);
F:=HD!.fpgroup; #F is an fp group rep of the n-th cohomology of (the small) D
                #So F has one generator for each basis element of D_n
FhomG:=NqEpimorphismNilpotentQuotient(F,1);
G:=Image(FhomG); #G is a pcg group rep of the n-th cohomology of D
gensG:=GeneratorsOfGroup(G);
preimgensG:=List(gensG,x->PreImagesRepresentative(FhomG,x));
preimgensG:=List(preimgensG,ExtRepOfObj);

###############################
f:=function(i)  #Inputs a generator of cohomology group and returns a cocycle 
                #on (the non-contracted) K
local w,v, m, j;
w:=preimgensG[i];
v:=[1..C!.dimension(n)]*0;

for j in [1..Length(w)/2] do
m:=HD!.h2c(w[2*j-1]);
v:=v+w[2*j]*m;
od;
#Print(SmallToBig(v,n),"\n\n");
return SmallToBig(v,n);
end;
###############################


###############################
g:=function(z) #Inputs a cocycle on K and returns a cohomology group element
local w;
return Exponents(Image(FhomG,HD!.c2h(BigToSmall(z,n))));
end;
###############################

fns[n+1]:=f;
gns[n+1]:=g;
clgy[n+1]:=Cohomology(D,n);
fi;
end;
##############################################
##############################################

for n in [0..Dimension(K)] do
Iterate(n);
od;

##############################################
##############################################
cup:=function(i,j,x,y)
local k, w, a, b, c, xx, yy;

w:=[1..CC!.dimension(i+j)]*0;
xx:=[1..CC!.dimension(i)]*0;
yy:=[1..CC!.dimension(j)]*0;
for k in [1..Length(clgy[i+1])] do
xx:=xx+x[k]*fns[i+1](k);
od;
for k in [1..Length(clgy[j+1])] do
yy:=yy+y[k]*fns[j+1](k);
od;

#Print(xx,"\n",yy,"\n\n");

for k in [1..CC!.dimension(i+j)] do
c:=K!.simplicesLst[i+j+1][k];
a:=c{[1..i+1]};
b:=c{[i+1..i+j+1]};
a:=Position(K!.simplicesLst[i+1],a);
b:=Position(K!.simplicesLst[j+1],b);
#Print([xx[a],yy[b]]);
w[k]:=xx[a]*yy[b];
od;

#Print(w,"\n");
return gns[i+j+1](w);

end;
##############################################
##############################################

return cup;
end);
#########################################################################
#########################################################################


##########################################################
##########################################################
InstallOtherMethod(CupProduct,
"integral cohomology cup product for a simplicial complex",
[IsHapSimplicialComplex],
function(Y)
return HAP_CupProductOfSimplicialComplex(Y);
end);
##########################################################
##########################################################

#BELOW HERE TO BE DELETED

#K:=ConnectedSum(ComplexProjectiveSpace(2),ComplexProjectiveSpace(2),-1);
K:=SimplicialK3Surface();
#K:=ClosedSurface(2);

cup:=HAP_CupProductOfSimplicialComplex(K);
h:=Length(Cohomology(RegularCWComplex(K),2));

A:=NullMat(h,h);;
for i in [1..h] do
 for j in [1..h] do
 x:=[1..h]*0;; x[i]:=1;;
 y:=[1..h]*0;; y[j]:=1;;
 A[i][j]:=cup(2,2,x,y)[1];
 od;od;

