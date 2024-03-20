#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(ModularCohomology,
function(X,n)
local
	Cohomology_Obj,
	Cohomology_Arr,
	CohomologyAsFpGroup;


#####################################################################
#####################################################################
Cohomology_Obj:=function(C,N)
local  
	M1, M2, 
	dim, 
	BasisKerd2, BasisImaged1, BasisKerd2one, one, prime,
	Rels, Smith, TorsionCoefficients,
	Dimension, Boundary,
	i,row,n,
	RankM1, RankM2,LengthM1,Rank;

#n:=N-1;

if N <0 then return rec(rank:=0); fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;
M1:=[];
M2:=[];

if N=0 then M2:=[List([1..Dimension(0)],i->0)];
else
for i in [1..Dimension(N-1)] do
M2[i]:=Boundary(N-1,i);
od;
fi;

if Length(M2)=0 then RankM2:=0; dim:=0; else
M2:=MutableCopyMat(M2);
ConvertToMatrixRep(M2);
RankM2:=RankMatDestructive(1*M2);
BasisImaged1:=BaseMat(M2);
dim:=Length(BasisImaged1);
fi;
M2:=0;


for i in [1..Dimension(N)] do
M1[i]:=Boundary(N,i);
od;

if Length(M1)=0 then RankM1:=0; else
M1:=MutableCopyMat(M1);
ConvertToMatrixRep(M1);
RankM1:=RankMatDestructive(1*M1);
fi;
LengthM1:=Length(M1);
if Length(M1)=0 then
BasisKerd2:=[];
Rels:=[];
else
BasisKerd2:=NullspaceMat(M1);
prime:=EvaluateProperty(C,"characteristic");
one:=One(GF(prime));
BasisKerd2one:=one*BasisKerd2;

M1:=0;

Rels:=[];
for i in [1..dim] do
        Rels[i]:=SolutionMat(BasisKerd2one,BasisImaged1[i]);
#Rels[i]:=SolutionMat(BasisKerd2,BasisImaged1[i]);
od;
fi;

Rank:= LengthM1-RankM1 -RankM2;;

return rec(
           basisKerd2:=BasisKerd2,
           rank:=Rank,
           rels:=Rels
           );


end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
CohomologyAsFpGroup:=function(C,n)
local  
	F, H, FhomH, Rels, Fgens, Frels, IHC, HhomC, ChomH,
	prime,FieldToInt, Vector2Word, BasisKerd2, BasisKerd2one, one,rel, i, j,
quot, Q, Qgens, iso, Qiso, Hgens;

IHC:=Cohomology_Obj(C,n);
BasisKerd2:=IHC.basisKerd2;
Rels:=IHC.rels;
prime:=EvaluateProperty(C,"characteristic");

F:=FreeGroup(Length(BasisKerd2));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];

one:=One(GaloisField(prime));
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


H:=F/Frels;
Hgens:=GeneratorsOfGroup(H);
FhomH:=GroupHomomorphismByImagesNC(F,H,Fgens,Hgens);

if Length(Hgens)=0 then
quot:=GroupHomomorphismByFunction(H,H,x->x); else
quot:=NqEpimorphismNilpotentQuotient(H,1);
fi;
Q:=Image(quot);
Qgens:=GeneratorsOfGroup(Q);
iso:=IsomorphismFpGroup(Q);
Qiso:=Image(iso);



#####################################################################
#HhomC:=function(w);
#return BasisKerd2[w];
#end;
#####################################################################
######################################################################
HhomC:=function(w)
local ww,v,i,s;

ww:=PreImagesRepresentative(quot,Qgens[w]);
ww:=PreImagesRepresentative(FhomH,ww);
v:=BasisKerd2[1]*0;
for i in [1..Length(BasisKerd2)] do
s:=ExponentSumWord(ww,Fgens[i]);
v:=v+s*BasisKerd2[i];
od;
return v ;
end;
#####################################################################



BasisKerd2one:=BasisKerd2*one;
#####################################################################
ChomH:=function(v)
local w;

w:=SolutionMat(BasisKerd2one,v);
w:=Vector2Word(w);
return Image(iso, Image(quot, Image(FhomH,w) ));
end;
#####################################################################

return rec(
	    fpgroup:=Qiso,    #Changed from H
	    h2c:=HhomC,
	    c2h:=ChomH );
end;
#####################################################################
#####################################################################



#####################################################################
#####################################################################
Cohomology_Arr:=function(f,n)
local
		C,D,ChomD,
		HC, HChomC, ChomHC, IHC,	
		HD, HDhomD, DhomHD, IHD,	
		HChomHD, gensHC, imageGensHC,
		x;

C:=f!.source;
D:=f!.target;
ChomD:=f!.mapping;

IHC:=CohomologyAsFpGroup(C,n);
HC:=IHC.fpgroup;
gensHC:=GeneratorsOfGroup(HC);
HChomC:=IHC.h2c;
ChomHC:=IHC.c2h;

IHD:=CohomologyAsFpGroup(D,n);
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

if X="CohomologyAsFpGroup" then return CohomologyAsFpGroup; fi;

if EvaluateProperty(X,"type")="cochainComplex" then
return Cohomology_Obj(X,n).rank; fi;

if EvaluateProperty(X,"type")="cochainMap" then
return Cohomology_Arr(X,n); fi;

Print("ERROR: Input should be a cochain complex or cochain map.\n");
end );
