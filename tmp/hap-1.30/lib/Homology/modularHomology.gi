#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
#####################################################################
InstallGlobalFunction(ModularHomology,
function(X,n)
local  
	FasterHomology_Obj,
	Homology_Obj,
	Homology_Arr,
	HomologyAsFpGroup;

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
        prime, field,
        i;

prime:=EvaluateProperty(C,"characteristic");
if IsPrimeInt(prime) then field:=GF(prime); else field:=Rationals; fi;
Dimension:=C!.dimension;
Boundary:=C!.boundary;

if n <0 then return rec( rank:=0); fi;
if Dimension(n)=0 then return rec( rank:=0); fi;

########################
if n=0 then
#BasisKerd1:=IdentityMat(Dimension(n))*One(field);
M1:=[];
fi;

if n>0 then
M1:=[];
 if Dimension(n-1)=0 then BasisKerd1:=Basis(field^Dimension(n));
 else

 for i in [1..Dimension(n)] do
 M1[i]:=Boundary(n,i);
ConvertToVectorRep(M1[i]);
 od;
 M1:=MutableCopyMat(M1);
 ConvertToMatrixRep(M1);
# BasisKerd1:=NullspaceMatDestructive(M1);
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
ConvertToVectorRep(M2[i]);
od;
if Length(M2)=0 then RankM2:=0; else
M2:=MutableCopyMat(M2);
ConvertToMatrixRep(M2);
RankM2:=RankMatDestructive(M2);fi;
M2:=0;

Rank:= LengthM1-RankM1 -RankM2;;

return rec( rank:=Rank);
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
        prime,field,
        i,k,tmp;

prime:=EvaluateProperty(C,"characteristic");
if IsPrimeInt(prime) then field:=GF(prime); else field:=Rationals; fi;
Dimension:=C!.dimension;
Boundary:=C!.boundary;

if n <0 then return false; fi;

########################
if n=0 then
BasisKerd1:=IdentityMat(Dimension(n))*One(field);
fi;

if n>0 then

M1:=[];

if Dimension(n-1)=0 then
  for i in [1..Dimension(n)] do
  M1[i]:=[Zero(field)];
  od;
else
  for i in [1..Dimension(n)] do
  M1[i]:=Boundary(n,i);
  od;
fi;
ConvertToMatrixRep(M1);
#BasisKerd1:=NullspaceMatDestructive(M1);
BasisKerd1:=NullspaceMat(M1);
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
 Add(tmp,[1..Dimension(n)]*Zero(field) );

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
        F, H, FHhomH, Rels, Fgens, Frels, IHC, HhomC, ChomH,
        Vector2Word, BasisKerd1, rel, i, j, prime, FieldToInt, one,
        sem, z1,sol1,lngm , Hgens, FH, epim, FHgens;

IHC:=Homology_Obj(C,n);
BasisKerd1:=IHC.basisKerd1;
Rels:=IHC.rels;
prime:=EvaluateProperty(C,"characteristic");

F:=FreeGroup(Length(BasisKerd1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];

one:=Elements(GaloisField(prime))[2];
#####################################################################
FieldToInt:=function(x)
local
	i;
for i in [0..prime] do
if i*one=x then return i; fi;
od;

end;
#####################################################################

#####################################################################
Vector2Word:=function(rel)
local w,i;

w:=Identity(F);
for i in [1..Length(Fgens)] do
w:=w*Fgens[i]^FieldToInt(rel[i]);
od;

return w;

end;
#####################################################################

for rel in Rels do
Append(Frels,[Vector2Word(rel)]);
od;


for i in [1..Length(Fgens)] do
Append(Frels,[Fgens[i]^prime]);
for j in [i..Length(Fgens)] do
Append(Frels,[Fgens[i]*Fgens[j]*Fgens[i]^-1*Fgens[j]^-1]);
od;
od;

#H:=F/Frels;
#FhomH:=GroupHomomorphismByImagesNC(F,H,Fgens,GeneratorsOfGroup(H));
#Hgens:=GeneratorsOfGroup(H);

FH:=F/Frels;
FHgens:=GeneratorsOfGroup(FH);

FHhomH:=NqEpimorphismNilpotentQuotient(FH,1);

H:=Range(FHhomH);
Hgens:=GeneratorsOfGroup(H);

epim:=EpimorphismFromFreeGroup(FH);


#####################################################################
#HhomC:=function(w);
#return BasisKerd1[w];
#end;
#####################################################################
#####################################################################
HhomC:=function(w)
local ww,v,i,s;
ww:=PreImagesRepresentative(FHhomH,Hgens[w]);
ww:=PreImagesRepresentative(epim,ww);
v:=BasisKerd1[1]*0*one;
for i in [1..Length(BasisKerd1)] do
s:=ExponentSumWord(ww,PreImagesRepresentative(epim,FHgens[i]));
v:=v+s*BasisKerd1[i]*one;
od;
return v ;
end;
#####################################################################



    z1 := Zero(BasisKerd1[1][1]);
    sol1:=ListWithIdenticalEntries(Length(BasisKerd1),z1);
    lngm:=Length(BasisKerd1[1]);
    sem := SemiEchelonMatTransformation(StructuralCopy(BasisKerd1));

#####################################################################
ChomH:=function(v)
local w,i,xx, vno, z,x, row, ncols,sol,vec,solmat;

### w:=SolutionMat(BasisKerd1,v) mod expH;
ncols := Length(BasisKerd1[1]);
vec:=v;
    z := StructuralCopy(z1);
    sol := StructuralCopy(sol1);
    ConvertToVectorRepNC(sol);
    for i in [1..ncols] do
        vno := sem.heads[i];
        if vno <> 0 then
            x := vec[i];
            if x <> z then
                AddRowVector(vec, sem.vectors[vno], -x);
                AddRowVector(sol, sem.coeffs[vno], x);
            fi;
        fi;
    od;

#xx:=One(H);
#for i in [1..Length(sol)] do
#xx:=xx*Hgens[i]^FieldToInt(sol[i]);
#od;

xx:=Identity(H);
for i in [1..Length(FHgens)] do

xx:=xx*Image(FHhomH,FHgens[i])^FieldToInt(sol[i]);

od;


return xx;
end;
#####################################################################


return rec(
            fpgroup:=H,
	    h2c:=HhomC,
	    c2h:=ChomH );
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
                x;


C:=f!.source;
D:=f!.target;
ChomD:=f!.mapping;

IHC:=HomologyAsFpGroup(C,n);
HC:=IHC.fpgroup;
gensHC:=GeneratorsOfGroup(HC);
HChomC:=IHC.h2c;
ChomHC:=IHC.c2h;

IHD:=HomologyAsFpGroup(D,n);
HD:=IHD.fpgroup;
HDhomD:=IHD.h2c;
DhomHD:=IHD.c2h;

imageGensHC:=[];
for x in [1..Length(gensHC)] do
Append(imageGensHC,[  DhomHD(ChomD(HChomC(x),n))  ]  );
od;

HChomHD:=GroupHomomorphismByImagesNC(HC,HD,gensHC,imageGensHC);
return HChomHD;
end;
#####################################################################
#####################################################################


				    
		

if EvaluateProperty(X,"type")="chainComplex" then
return FasterHomology_Obj(X,n).rank; fi;

if EvaluateProperty(X,"type")="chainMap" then
return Homology_Arr(X,n); fi;

end);
#####################################################################
#####################################################################
#####################################################################
