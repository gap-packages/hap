#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IntegralHomology,
function(X,n)
local
	Homology_Obj,
	Homology_Arr,
	HomologyAsFpGroup;

if not EvaluateProperty(X,"characteristic")=0 then 
Print("ERROR: There is an inconsitency with characteristic of Z. \n");
return fail; fi;

#####################################################################
#####################################################################
Homology_Obj:=function(C,n)
local  
	M1, M2, 
	dim, 
	BasisKerd1, BasisImaged2, BasisKerd1cp, BasisImaged2cp, 
	Rels, Smith, TorsionCoefficients,
	Dimension, Boundary,
	i;

if n <0 then return false; fi;

Dimension:=C!.dimension;
Boundary:=C!.boundary;


########################
if n=0 then 
BasisKerd1:=IdentityMat(Dimension(n));

else
M1:=[];

for i in [1..Dimension(n)] do
M1[i]:=Boundary(n,i);
od;
ConvertToMatrixRep(M1);
BasisKerd1:=LLLReducedBasis(M1,"linearcomb").relations;
#BasisKerd1:=NullspaceIntMat(M1);
M1:=0;

fi;
#######################


M2:=[];
for i in [1..Dimension(n+1)] do
M2[i]:=Boundary(n+1,i);
od;

ConvertToMatrixRep(M2);
#BasisImaged2:=LLLReducedBasis(M2).basis;
BasisImaged2:=BaseIntMat(M2);
dim:=Length(BasisImaged2);
M2:=0;

#Rels:=[];
#for i in [1..dim] do
#	Rels[i]:=SolutionMat(BasisKerd1,BasisImaged2[i]);
#od;

#BasisKerd1:=MutableCopyMat(BasisKerd1);
#ConvertToMatrixRep(BasisKerd1);
BasisImaged2:=MutableCopyMat(List([1..dim],i->BasisImaged2[i]));
#!!!!!!!!!!!!!!!CHANGE IT
#ConvertToMatrixRep(BasisImaged2);

if Length(BasisImaged2)>0 then
Rels:=SolutionsMatDestructive(BasisKerd1, BasisImaged2);
else
Rels:=[];
fi;

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
HomologyAsFpGroup:=function(C,n)
local  
	F, H, HH, FhomH, FhomHH, HHhomH, nn,
	Rels, Fgens, Frels, IHC, HhomC, ChomH, ExpH,
	Vector2Word, BasisKerd1, rel, i, j, Htmp,FhomHtmp,HtmphomH;

if not "fpIntHom" in NamesOfComponents(C) then
C!.fpIntHom:=[1..1000];        #SLOPPPY! Some one might ask for 
			      #the 1000-dimensional homology
fi;

if n=0 then nn:=1000; else nn:=n; fi;

if IsInt(C!.fpIntHom[nn]) then

IHC:=Homology_Obj(C,n);
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

if true then #not LoadPackage("nq") then  #This does not work too well!
#################################	  #Should make "abelian groups"
for i in [1..Length(Fgens)] do		  #type in GAP
for j in [i..Length(Fgens)] do
Append(Frels,[Fgens[i]*Fgens[j]*Fgens[i]^-1*Fgens[j]^-1]);
od;
od;

H:=F/Frels;
SetIsAbelian(H,true);
ExpH:=AbelianInvariants(H);
  if Length(ExpH)>0 then ExpH:=Maximum(ExpH); else ExpH:=1; fi;
ExpH:=0;
FhomH:=GroupHomomorphismByImagesNC(F,H,Fgens,GeneratorsOfGroup(H));
################################
else
################################
HH:=F/Frels;
ExpH:=0;
ExpH:=AbelianInvariants(HH);
  if Length(ExpH)>0 then ExpH:=Maximum(ExpH); else ExpH:=1; fi;
FhomHH:=GroupHomomorphismByImagesNC(F,HH,Fgens,GeneratorsOfGroup(HH));
HHhomH:=HAP_NqEpimorphismNilpotentQuotient(HH,1);;
H:=Range(HHhomH);
FhomH:=GroupHomomorphismByImagesNC(F,H,Fgens,
List(Fgens,x->Image(HHhomH,Image(FhomHH,x)) ));
################################
fi;

#####################################################################
HhomC:=function(w);
return BasisKerd1[w];
end;
#####################################################################

if Length(BasisKerd1)=0 then
#####################################################################
ChomH:=function(v)  #Am I sure about this fix?
return Identity(H);
end;
#####################################################################
else
#####################################################################
ChomH:=function(v)
local w;

if ExpH>0 then
w:=SolutionMat(BasisKerd1,v) mod ExpH; #Am I sure about this?
else
w:=SolutionMat(BasisKerd1,v);
fi;
w:=Vector2Word(w) ; 
return Image(FhomH,w);
end;
#####################################################################
fi;

C!.fpIntHom[nn]:=rec(
	    fpgroup:=H,
	    h2c:=HhomC,
	    c2h:=ChomH );
fi;

return C!.fpIntHom[nn];
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
#return Homology_Obj(X,n).torsionCoefficients; fi; #unecessarily slow!!
return IntegralHomologyOfChainComplex(X,n); fi;

if EvaluateProperty(X,"type")="chainMap" then
return Homology_Arr(X,n); fi;

Print("ERROR: Input should be a chain complex or chain map.\n");
end );
