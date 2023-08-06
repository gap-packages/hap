#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IntegralCohomology,
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
	BasisKerd1, BasisImaged2, 
	Rels, Smith, TorsionCoefficients,
	Dimension, Boundary,
	i,row,n;

n:=N-1;

if N <0 then return rec(torsionCoefficients:=[ ]); fi;
#if N=0 then return 
#rec(
#	torsionCoefficients:=[0],
#	); 
#fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;
M1:=[];
M2:=[];

if n=-1 then M2:=[List([1..Dimension(0)],i->0)];
else
for i in [1..Dimension(n)] do
M2[i]:=Boundary(n,i);
od;
#M2:=TransposedMat(M2);
fi;

for i in [1..Dimension(n+1)] do
M1[i]:=Boundary(n+1,i);
od;
#M1:=TransposedMat(M1);


BasisKerd1:=LLLReducedBasis(M1,"linearcomb").relations;
#BasisImaged2:=LLLReducedBasis(M2).basis;
BasisImaged2:=BaseIntMat(M2);
dim:=Length(BasisImaged2);

Rels:=[];
for i in [1..dim] do
	Rels[i]:=SolutionMat(BasisKerd1,BasisImaged2[i]);
od;

Smith:= SmithNormalFormIntegerMat(Rels);

TorsionCoefficients:=[];

for i in [1..Length(BasisKerd1)] do
	if i<=Length(Smith) then
	TorsionCoefficients[i]:=Smith[i][i];
	else 
	TorsionCoefficients[i]:=0;
	fi;
od;

return rec(
	   torsionCoefficients:=Filtered(TorsionCoefficients, i-> not (i=1)),
	   rels:=Rels,
	   basisKerd1:=BasisKerd1);

end;
#####################################################################
#####################################################################


#####################################################################
#####################################################################
CohomologyAsFpGroup:=function(C,n)
local  
	F, H, FhomH, Rels, Fgens, Frels, IHC, HhomC, ChomH,
	Vector2Word, BasisKerd1, rel, Hgens, Vector2Word2, i, j,
        quot, Q, iso, Qiso, Qgens;

IHC:=Cohomology_Obj(C,n);
BasisKerd1:=IHC.basisKerd1;
Rels:=IHC.rels;

F:=FreeGroup(Length(BasisKerd1));
Fgens:=GeneratorsOfGroup(F);
Frels:=[];

#####################################################################
Vector2Word:=function(rel)
local w,i;

w:=Identity(F);
for i in [1..Length(Fgens)] do
w:=w*Fgens[i]^rel[i];		
od;

return w;
end;
#####################################################################

for rel in Rels do
Append(Frels,[Vector2Word(rel)]);
od;

#I don't think we need these any more
#for i in [1..Length(Fgens)] do
#for j in [i..Length(Fgens)] do
#Append(Frels,[Fgens[i]*Fgens[j]*Fgens[i]^-1*Fgens[j]^-1]);
#od;
#od;

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
Vector2Word2:=function(rel)
local w,i;

w:=Identity(H);
for i in [1..Length(Hgens)] do
w:=w*Hgens[i]^rel[i];
od;

return w;
end;
#####################################################################


#####################################################################
#HhomC:=function(w);
#return BasisKerd1[w];
#end;
#####################################################################
######################################################################
HhomC:=function(w)
local ww,v,i,s;

ww:=PreImagesRepresentative(quot,Qgens[w]);
ww:=PreImagesRepresentative(FhomH,ww);
v:=BasisKerd1[1]*0;
for i in [1..Length(BasisKerd1)] do
s:=ExponentSumWord(ww,Fgens[i]);
v:=v+s*BasisKerd1[i];
od;
return v ;
end;
#####################################################################


#####################################################################
ChomH:=function(v)
local w;
w:=SolutionMat(BasisKerd1,v);
w:=Vector2Word2(w);
return Image(iso, Image(quot, w) );
end;
#####################################################################

return rec(
	    fpgroup:=Qiso,   #Changed from H
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
		epiNHC,epiNHD,NHC,NHD,x, gens, gensNHC;

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


if Length(GeneratorsOfGroup(HD))=0 then
return GroupHomomorphismByImagesNC(HC,HD,gensHC,List(gensHC,x->One(HD)));
fi;   ####ADDED FEB 2014



imageGensHC:=[];
for x in [1..Length(gensHC)] do
Append(imageGensHC,[  DhomHD(ChomD(HChomC(x),n))  ]  );
od;




HChomHD:=GroupHomomorphismByImagesNC(HC,HD,gensHC,imageGensHC);
#return HChomHD;



epiNHC:=NqEpimorphismNilpotentQuotient(HC,1);
epiNHD:=NqEpimorphismNilpotentQuotient(HD,1);
NHC:=Range(epiNHC);
NHD:=Range(epiNHD);



#return GroupHomomorphismByFunction(NHC,NHD,x ->
#Image(epiNHD,Image(HChomHD, PreImagesRepresentative(epiNHC,x))) );


gensNHC:=GeneratorsOfGroup(NHC);
gens:=List(gensNHC,x->Image(epiNHD,Image(HChomHD, PreImagesRepresentative(epiNHC,x)))  );
return GroupHomomorphismByImagesNC(NHC,NHD,gensNHC,gens);
end;
#####################################################################
#####################################################################

if X="CohomologyAsFpGroup" then return CohomologyAsFpGroup; fi;

if EvaluateProperty(X,"type")="cochainComplex" then
#return Cohomology_Obj(X,n).torsionCoefficients; fi;
return IntegralCohomologyOfCochainComplex(X,n); fi;

if EvaluateProperty(X,"type")="cochainMap" then
return Cohomology_Arr(X,n); fi;

Print("ERROR: Input should be a cochain complex or cochain map.\n");
end );
