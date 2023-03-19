#(C) Graham Ellis, 2005-2006
#RT:=0;
#####################################################################
InstallGlobalFunction(IntegralHomology,
function(X,n)
local
	Homology_Obj,
	Homology_Arr,
	HomologyAsFpGroup, xx;


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
if n=0  or Dimension(n-1)=0 then    #CHANGED March 2023
BasisKerd1:=IdentityMat(Dimension(n));

else
M1:=[];

for i in [1..Dimension(n)] do
M1[i]:=Boundary(n,i);
od;
ConvertToMatrixRep(M1);
BasisKerd1:=LLLReducedBasis(M1,"linearcomb").relations;
M1:=0;

fi;
#######################


M2:=[];
for i in [1..Dimension(n+1)] do
M2[i]:=Boundary(n+1,i);
od;

ConvertToMatrixRep(M2);
if M2=[[]] then BasisImaged2:=[]; else
BasisImaged2:=BaseIntMat(M2); fi;   #CHANGED 2023
dim:=Length(BasisImaged2);
M2:=0;


BasisImaged2:=MutableCopyMat(List([1..dim],i->BasisImaged2[i]));
#!!!!!!!!!!!!!!!CHANGE IT

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
	F, H, HH, HHhomH, nn, expH,
	Rels, Fgens, FhomFH, FHhomH, FhomH, 
        FH, Frels, Hgens, FHgens, IHC, HhomC, ChomH, 
	Vector2Word, BasisKerd1, rel, i, j, Htmp,HtmphomH,
        sem, z1,sol1,lngm , epim;

if not "fpIntHom" in NamesOfComponents(C) then
C!.fpIntHom:=[1..1000];       #SLOPPPY! Some one might ask for 
			      #the 1000-dimensional homology
fi;

if n=0 then nn:=1000; else nn:=n; fi;

if IsInt(C!.fpIntHom[nn]) then

IHC:=Homology_Obj(C,n);
BasisKerd1:=IHC.basisKerd1;
Rels:=IHC.rels;
if Length(IHC.torsionCoefficients)>0 then  #changed August 2017
expH:=Lcm(IHC.torsionCoefficients);
else expH:=0;
fi;

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

#################################	  #Should make "abelian groups"
for i in [1..Length(Fgens)] do		  #type in GAP
for j in [i..Length(Fgens)] do
Append(Frels,[Fgens[i]*Fgens[j]*Fgens[i]^-1*Fgens[j]^-1]);
od;
od;

FH:=F/Frels;
FHgens:=GeneratorsOfGroup(FH);

FHhomH:=NqEpimorphismNilpotentQuotient(FH,1);
#FHhomH:=EpimorphismNilpotentQuotient(FH); #Commented out 20 June 2018


H:=Range(FHhomH);
Hgens:=GeneratorsOfGroup(H);
################################

epim:=EpimorphismFromFreeGroup(FH);

#####################################################################
HhomC:=function(w)
local ww,v,i,s;
ww:=PreImagesRepresentative(FHhomH,Hgens[w]);
ww:=PreImagesRepresentative(epim,ww);
v:=BasisKerd1[1]*0;
for i in [1..Length(BasisKerd1)] do
s:=ExponentSumWord(ww,PreImagesRepresentative(epim,FHgens[i]));
v:=v+s*BasisKerd1[i];
od;
return v ;
end;
#####################################################################

if Length(BasisKerd1)=0 then
#####################################################################
ChomH:=function(v)  #Am I sure about this fix?
return Identity(H);
end;
#####################################################################
else
    
    z1 := Zero(BasisKerd1[1][1]);
    sol1:=ListWithIdenticalEntries(Length(BasisKerd1),z1);
    lngm:=Length(BasisKerd1[1]);
    sem := SemiEchelonMatTransformationDestructive(1*BasisKerd1);

#####################################################################
ChomH:=function(v)
local w,i,xx, vno, z,x, row, ncols,sol,vec,solmat;

### thenw:=SolutionMat(BasisKerd1,v) mod expH; 
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
if expH>0 then sol:=sol  mod expH; fi;


xx:=Identity(H);
for i in [1..Length(FHgens)] do

xx:=xx*Image(FHhomH,FHgens[i])^sol[i];
od;


return xx;
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

#################################### ADDED MARCH 2023
if C!.dimension(n)=0 then
HC:=FreeGroup(0);
fi;
if D!.dimension(n)=0 then
HD:=FreeGroup(0);
fi;
if IsBound(HC) and not IsBound(HD) then
HD:=HomologyAsFpGroup(D);
fi;
if IsBound(HD) and not IsBound(HC) then
HC:=HomologyAsFpGroup(D);
fi;
if IsBound(HC) and IsBound(HD) then return
GroupHomomorphismByFunction(HC,HD,x->x);
fi;
####################################

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

if X="HomologyAsFpGroup" then return HomologyAsFpGroup; fi;

if not EvaluateProperty(X,"characteristic")=0 then
Print("ERROR: There is an inconsitency with characteristic of Z. \n");
return fail; fi;

if EvaluateProperty(X,"type")="chainComplex" then
xx:=IntegralHomologyOfChainComplex(X,n);
return xx; fi;

if EvaluateProperty(X,"type")="chainMap" then
xx:=Homology_Arr(X,n);
return xx; fi;

Print("ERROR: Input should be a chain complex or chain map.\n");
end );
