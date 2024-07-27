#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(HomologyVectorSpace,
function(X,n)
local  
	FasterHomology_Obj,
	Homology_Obj,
	Homology_Arr,
	HomologyAsFpGroup,
        prime,Fld,C;

#if IsHapChainMap(X) then C:=X; else C:=Source(X); fi;
if IsHapChainComplex(X) then C:=X; else C:=Source(X); fi;
prime:=EvaluateProperty(C,"characteristic");

if prime=-1/2 then Fld:=Field(1); fi;
if IsPrimeInt(prime) then Fld:=GF(prime); fi;
if not IsPrimeInt(prime) and not prime=-1/2 then
Print("The field should be the rationals or of prime characteristic.\n");
return fail;
fi;

#####################################################################
#####################################################################
FasterHomology_Obj:=function(C,n)
local
        M1, M2,
        dim,
        rankM1, rankM2,
        Dimension, Boundary,
        BasisKerd1, BasisImaged2, Rels, Rank, RankM1, RankM2,
	LengthM1,LengthM2,
        i;

Dimension:=C!.dimension;
Boundary:=C!.boundary;

if n <0 then return 0; fi;
if Dimension(n)=0 then return Fld^0; fi;

########################
if n=0 then
M1:=[];
fi;

if n>0 then
M1:=[];
 if Dimension(n-1)=0 then BasisKerd1:=Basis(Fld^Dimension(n)); 
 else

 for i in [1..Dimension(n)] do
 M1[i]:=Boundary(n,i);
 od;
 M1:=MutableCopyMat(M1);
 ConvertToMatrixRep(M1);
 fi;
fi;
#######################

if Length(M1)=0 then RankM1:=0; else
RankM1:=RankMatDestructive(M1);
fi;
LengthM1:=Dimension(n);
M1:=0;

M2:=[];
for i in [1..Dimension(n+1)] do
M2[i]:=Boundary(n+1,i);
od;
if Length(M2)=0 then RankM2:=0; else
M2:=MutableCopyMat(M2);
ConvertToMatrixRep(M2);
RankM2:=RankMatDestructive(M2);fi;
M2:=0;

Rank:= LengthM1-RankM1 -RankM2;;

return Fld^Rank;
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
Homology_Obj:=function(C,n)
local
	M1, M2, 
	dim, 
	rankM1, rankM2, 
	Dimension, Boundary,
	BasisKerd1, BasisImaged2, Rels, Rank, 
	i,k,tmp;

Dimension:=C!.dimension;
Boundary:=C!.boundary;

if n <0 then return false; fi;
if Dimension(n)=0 then return Fld^0; fi;
########################
if n=0 then
BasisKerd1:=IdentityMat(Dimension(n))*One(Fld);
fi;

if n>0 then

M1:=[];

if Dimension(n-1)=0 then
  for i in [1..Dimension(n)] do
  M1[i]:=[Zero(Fld)];
  od;
else
  for i in [1..Dimension(n)] do
  M1[i]:=Boundary(n,i);
  od;
fi;
ConvertToMatrixRep(M1);
#BasisKerd1:=NullspaceMatDestructive(M1);
#BasisKerd1:=NullspaceMat(M1);
BasisKerd1:=TriangulizedNullspaceMat(M1);
M1:=0;

fi;
#######################

M2:=[];
tmp:=[];
k:=0;

 while k<Dimension(n+1) do
 M2:=[];
 for i in [k+1..Minimum(Dimension(n+1),k+100)] do
 M2[i-k]:=Boundary(n+1,i);
 od;

 k:=Minimum(Dimension(n+1),k+100);
 Append(tmp,BaseMat(M2));
 od;
 Add(tmp,[1..Dimension(n)]*Zero(Fld) );

BasisImaged2:=BaseMat(tmp); tmp:=0;
dim:=Length(BasisImaged2);

Rels:=[];
for i in [1..dim] do
        Rels[i]:=SolutionMat(BasisKerd1,BasisImaged2[i]);
od;

Rank:=Length(BasisKerd1) - dim;

return rec(
	   basisKerd1:=BasisKerd1,
	   rank:=Rank,
	   rels:=Rels);
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
HomologyAsFpGroup:=function(C,n)
local
        F, H, FhomH, Rels, Fgens, Frels, IHC, HhomC, ChomH,
        Vector2Word, BasisKerd1, rel, i, j, FieldToInt, 
	one, Bas;

#############SAVE RECOMPUTATIONS##########
if "HomVecSpc" in NamesOfComponents(C) then
if IsBound(C!.HomVecSp[n+1]) then
return C!.HomVecSp[n+1];
fi;
else
C!.HomVecSp:=[];
fi;
##########################################

IHC:=Homology_Obj(C,n);
BasisKerd1:=IHC.basisKerd1;
Rels:=IHC.rels;

F:=Fld^Length(BasisKerd1);
Fgens:=List(CanonicalBasis(F),x->x);
Frels:=[Zero(F)];

one:=One(Fld);
#####################################################################
FieldToInt:=function(x)
local
	i;
if prime=-1/2 then return x; fi;
for i in [0..prime] do
if i*one=x then return i; fi;
od;

end;
#####################################################################

#####################################################################
Vector2Word:=function(rel)
local w,i;

return rel*one;  

end;
#####################################################################

for rel in Rels do
Append(Frels,[Vector2Word(rel)]);
od;

if Dimension(F)=0 then Frels:=[]; fi;
FhomH:=NaturalHomomorphismBySubspace(F,Subspace(F,Frels));
H:=Range(FhomH);
Bas:=List(CanonicalBasis(H),a->a);

#####################################################################
HhomC:=function(w)
local x,i,c;

x:=PreImagesRepresentative(FhomH,Bas[w]);
c:=BasisKerd1[1]*0;

for i in [1..Length(x)] do
c:=c+BasisKerd1[i]*x[i];
od;


return c;
end;
#####################################################################


#####################################################################
ChomH:=function(v)
local w;
w:=SolutionMat(BasisKerd1,v);
w:=Vector2Word(w); 
w:=Image(FhomH,w);
return w;
end;
#####################################################################

C!.HomVecSp[n+1]:=
	     rec(
            fpgroup:=H,
	    h2c:=HhomC,
	    c2h:=ChomH );

return C!.HomVecSp[n+1];
end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
Homology_Arr:=function(f,n)
local
                C,D,ChomD,
	        HC, HChomC, ChomHC, IHC,
                HD, HDhomD, DhomHD, IHD,
                HChomHD, gensHC, imageGensHC,
                x,V,W;


C:=f!.source;
D:=f!.target;
ChomD:=f!.mapping;

########################################
if C!.dimension(n)=0 and D!.dimension(n)=0 then  
V:=Fld^0;
return NaturalHomomorphismBySubspace(V,V);
fi;

if C!.dimension(n)=0 then
W:=HomologyAsFpGroup(D,n)!.fpgroup;
if Dimension(W)=0 then
V:=Fld^0;
return NaturalHomomorphismBySubspace(V,V);
fi;
V:=Subspace(W,[Zero(W)]);
return LeftModuleHomomorphismByImagesNC(V,W,GeneratorsOfVectorSpace(V),[Zero(W)]);
fi;

if D!.dimension(n)=0 then
V:=HomologyAsFpGroup(C,n)!.fpgroup;
return NaturalHomomorphismBySubspace(V,V);
fi;
########################################


IHC:=HomologyAsFpGroup(C,n);
HC:=IHC.fpgroup;
gensHC:=List(CanonicalBasis(HC),x->x);
HChomC:=IHC.h2c;
ChomHC:=IHC.c2h;

IHD:=HomologyAsFpGroup(D,n);
HD:=IHD.fpgroup;
HDhomD:=IHD.h2c;
DhomHD:=IHD.c2h;

imageGensHC:=[];

for x in [1..Length(gensHC)] do 
Append(imageGensHC,[  DhomHD(ChomD(HChomC(x),n))  ]  );
#ChomD(a,n)  take most of the time
od;


#HChomHD:=GroupHomomorphismByImagesNC(HC,HD,gensHC,imageGensHC);
if Dimension(HD)=0 then
HChomHD:=NaturalHomomorphismBySubspace(HC,HC);
else
HChomHD:=LeftModuleHomomorphismByImagesNC( HC, HD, gensHC, imageGensHC );
fi;

return HChomHD;
end;
#####################################################################
#####################################################################

if EvaluateProperty(X,"type")="chainComplex" then
 if IsPrimeInt(EvaluateProperty(X,"characteristic")) or
 EvaluateProperty(X,"characteristic")=-1/2 then
 return FasterHomology_Obj(X,n); 
 else Print("This function can only be applied over the rationals or in prime characteristic.\n");
 return fail; fi;
fi;

if EvaluateProperty(X,"type")="chainMap" then
 if IsPrimeInt(EvaluateProperty(Source(X),"characteristic")) or
 EvaluateProperty(Source(X),"characteristic") =-1/2 then
 return Homology_Arr(X,n); 
 else Print("This function can only be applied over the rationals or in prime characteristic.\n");
 return fail; fi;
fi;

end);
#####################################################################
#####################################################################
#####################################################################
