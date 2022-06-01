
################################################
################################################
InstallGlobalFunction(LinkingForm,
function(KK,kk)
local k,K,Y, dim, pdhom, Equiv, C, CC, D, E1, IntCoh,  
      SmallToBigCo,  
      BigToSmallHo, IntHom, H, imgensG, Dual1, Dual2,
      fnk, HD, HC, F, G, FhomG, gensG, preimgensG,n,duality,
      fc, Dual2fn, Bk, LnkForm ;

if IsHapRegularCWComplex(KK) then K:=BarycentricSubdivision(KK);
else K:=KK; fi;

k:=AbsInt(kk);
Y:=RegularCWComplex(K); 
dim:=Dimension(Y);
Equiv:=ChainComplexEquivalenceOfRegularCWComplex(Y);
C:=Source(Equiv[2]);  ##This has fewer generators than cells of Y
CC:=Target(Equiv[2]);  ##This has generators equal to the cells of Y
fc:=TransposedMat(BoundaryMatrix(CC,dim));
fc:=NullspaceIntMat(fc)[1];  #fundamental class of manifold
fc:=SignInt(kk)*fc;
D:=HomToIntegers(C);;
IntCoh:=IntegralCohomology("CohomologyAsFpGroup",true);;
IntHom:=IntegralHomology("HomologyAsFpGroup",true);;

Bk:=BoundaryMatrix(C,k);

E1:=HomToIntegers(Equiv[1]);
SmallToBigCo:=E1!.mapping;
    #Inputs a cohomology vector v of length C!.dimension(k) and
    #returns a vector of length Y!.dimension(k).
BigToSmallHo:=Equiv[1]!.mapping;
    #Inputs a homology vector v of length Y!.dimension(k) and
    #returns a vector of length C!.dimension(k).

HC:=IntHom(C,dim-k);
HD:=IntCoh(D,k);

F:=HD!.fpgroup; #F is an fp group rep of the k-th cohomology of (the small) D
                #So F has one generator for each basis element of D_k

FhomG:=NqEpimorphismNilpotentQuotient(F,1);

G:=Image(FhomG); #G is a pcg group rep of the k-th cohomology of D

gensG:=GeneratorsOfGroup(G);
preimgensG:=List(gensG,x->PreImagesRepresentative(FhomG,x));

preimgensG:=List(preimgensG,ExtRepOfObj);

H:=HC!.fpgroup; #H is a pcp-group rep of the dim-k-th homology of (the small) C



###############################
fnk:=function(i)  #Inputs a generator of k-th cohomology group and returns 
                 #a cocycle on (the non-contracted) K
local w,v, m, j;
w:=preimgensG[i];
v:=[1..C!.dimension(k)]*0;

for j in [1..Length(w)/2] do

m:=HD!.h2c(w[2*j-1]);

v:=v+w[2*j]*m;
od;

return SmallToBigCo(v,k);
end;
###############################


###############################
#function inputs a k-cohomology generator x and returns an n-k-cycle z
duality:=function(g)
local x, xx,z, i, j, s, a, b;

x:=Exponents(g);   #g is an element in a pcp-group
xx:=[1..CC!.dimension(k)]*0;
for j in [1..Length(x)] do
xx:=xx+x[j]*fnk(j);
od;

z:=[1..CC!.dimension(dim-k)]*0; 

for i in [1..K!.nrSimplices(dim)] do
s:=K!.simplicesLst[dim+1][i];
a:=s{[1..k+1]};
b:=s{[k+1..dim+1]};
a:=Position(K!.simplicesLst[k+1],a);
b:=Position(K!.simplicesLst[dim-k+1],b);

z[b]:=z[b] + fc[i]*xx[a];
od;
return HC.c2h(BigToSmallHo(z,dim-k));
end;
################################

imgensG:=List(gensG,duality);
Dual1:=GroupHomomorphismByImages(G,H,gensG,imgensG);
if Dual1=fail then Print("Poincare duality fails.\n"); return fail; fi;
Dual2:=InverseGeneralMapping(Dual1);


#############################
Dual2fn:=function(h)
local w,i, z;
z:=[1..C!.dimension(k)]*0;
w:=ExtRepOfObj(PreImagesRepresentative(FhomG,Image(Dual2,h)))  ;
for i in [1..Length(w)/2] do
z:=z+HD!.h2c(w[2*i-1])*w[2*i];
od;
return SolutionMat(Bk,z);
end;
#############################

#################################
LnkForm:=function(hh,gg)
local g,ex, v, m, a, b, i;
ex:=Exponents(hh);
v:=[1..C!.dimension(dim-k)]*0;
for i in [1..Length(ex)] do
v:=v+HC!.h2c(i)*ex[i];
od;
g:=Dual2fn(gg);
m:=v*g;
a:=NumeratorRat(m);
b:=DenominatorRat(m);
a:=a mod b;
return a/b;
end;
#################################


return rec(LinkingForm:=LnkForm, Homology:=H);
end);
################################################
################################################

################################################
################################################
InstallGlobalFunction(LinkingFormHomotopyInvariant,
function(W)
local L, Lnk, H, h, I, J;
L:=LinkingForm(W,2);
Lnk:=L!.LinkingForm;
H:=L!.Homology;

H:=L!.Homology;

I:=[];
for h in H do
Add(I,Lnk(h,h));
od;


L:=LinkingForm(W,-2);
Lnk:=L!.LinkingForm;
H:=L!.Homology;
J:=[];
for h in H do
Add(J,Lnk(h,h));
od;

return SortedList([SortedList(I), SortedList(J)]);
end);
################################################
################################################

