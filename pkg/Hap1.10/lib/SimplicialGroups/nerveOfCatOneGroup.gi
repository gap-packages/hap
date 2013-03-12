#InstallOtherMethod( GroupHomomorphismByImagesNC, "for group with no generators",
#    [IsGroup,IsGroup,IsEmpty,IsEmpty], SUM_FLAGS,
#        function(g,h,gg,gh)
#    return GroupHomomorphismByFunction(g,h,x->One(h),false,x->One(g));
#end);

########################################################################
########################################################################
InstallGlobalFunction(NerveOfCatOneGroup,
function(X,n)
local NerveOfCatOneGroup_Homopre,
	  C,NerveOfCatOneGroup_Obj,
	  NerveOfCatOneGroup_Homo,
	  NerveOfCatOneGroup_Seq;
	  
#####################################################################
#####################################################################
NerveOfCatOneGroup_Obj:=function(C,number)
local 
    ListGroups,Boundaries, Degeneracies, Nerve,
    smap,tmap,e,
	G,M,H,
	AutM,phi,
	g,m,tempprod,ConjugatorOfProd,tempBound,tempDegen,
	ComposOfGens,Gens,TempB,ImageOfGens,TempL,
	Elto,boundary,degeneracy,
	i,j,k,n,
	ElementsOfSemiDirect,BoundaryElement,CreatElement,DegenElement,
	GroupsList,BoundariesList,DegeneraciesList;

	
if not IsHapCatOneGroup(C) then
	Print("This function must be applied to a cat-1-group.\n");
	return fail;
fi;

#################
ListGroups:=[];
Boundaries:=[];
Degeneracies:=[];

smap:=C!.sourceMap;
tmap:=C!.targetMap;
Add(Boundaries,[smap,tmap]);
G:=smap!.Source;
M:=Kernel(smap);
e:=Identity(M);
Add(ListGroups,Image(tmap));
Add(ListGroups,G);
AutM:=AutomorphismGroup(M);
Gens:=GeneratorsOfGroup(G);
phi:=GroupHomomorphismByImagesNC(G,AutM,Gens,List(Gens,g->ConjugatorAutomorphismNC(M,Image(tmap,g))));
H:=SemidirectProduct(G,phi,M);
Add(ListGroups,H);

##   Create [m1,m2, ...mn] --> semidirect product of n elements
####################################
CreatElement:=function(ListM)  
local i,G,
    m,len;
	len:=Length(ListM);
if len=1 then
		m:=ListM[1];
fi;
if len>1 then
	m:=ListM[1];
	for i in [2..len] do    
		G:=ListGroups[i+1];
		m:=Image(Embedding(G,1),m)*Image(Embedding(G,2),ListM[i]);
    od;	
fi;
return m;
end;
########################################
ElementsOfSemiDirect:= function(H,n)  
local 
    gens,componentsofgens,tempL,
	emb1,emb2,pro,K,
	x,y,m,len,
	i,j;
	gens:=GeneratorsOfGroup(H);
	componentsofgens:=[];
	len:=Length(gens);
	if n=1 then
		for i in [1..len] do
			componentsofgens[i]:=[gens[i]];
		od;
	fi;
	if n>1 then 
		for i in [1..len] do
			K:=H;
			tempL:=[];
			y:=gens[i];
			for j in [1..n-1] do
				pro:=Projection(K);
				emb1:=Embedding(K,1);
				emb2:=Embedding(K,2);
				m:=PreImagesRepresentative(emb2,(Image(emb1,Image(pro,y)))^(-1)*y);	
				tempL[n-j+1]:=m;
				y:=Image(pro,y);
				K:=emb1!.Source;	
			od;
			tempL[1]:=y;
			componentsofgens[i]:=tempL;
		od;
	fi;
	return [gens,componentsofgens];
end;

####################################
BoundaryElement:=function(ListM)
local n,i,j,tempB,Bound;
n:=Length(ListM); 
if n=2 then 
	Bound:=[[Image(tmap,ListM[1])*ListM[2]],[ListM[1]*ListM[2]],[ListM[1]]];	
fi;

if n>2 then
	Bound:=[];
######Creat 1###############
	tempB:=[ Image(tmap,ListM[1])*ListM[2] ];
	for i in [2..n-1] do
		tempB[i]:=ListM[i+1];
	od;
	Add(Bound,tempB);
##### Creat 2 --> n#############
	for i in [2..n] do
		tempB:=[];
		for j in [1..i-2] do
			tempB[j]:=ListM[j];
			od;
		tempB[i-1]:= ListM[i-1]*ListM[i];
		for j in [i..n-1] do
			tempB[j]:=ListM[j+1];
		od;
	Add(Bound,tempB);
	od;
####### Creat n+1##########
	tempB:=[];
	for i in [1..n-1] do
		tempB[i]:=ListM[i];
	od;
	Add(Bound,tempB);
fi;

return Bound;
end;
###############################################################################

DegenElement:=function(ListM)
local n,i,j,tempB,Degen,g;

n:=Length(ListM);
g:=ListM[1]; 
if n=1 then 
	Degen:=[[Image(smap,g),Image(smap,g^(-1))*g],[g,e]];	
fi;

if n>1 then
	Degen:=[];
#####Creat 1############
	tempB:=[Image(smap,g),Image(smap,g^(-1))*g];
	for i in [3..n+1] do
		tempB[i]:=ListM[i-1];
	od;
	Add(Degen,tempB);
#####Creat 2 --> n+1########
	for i in [2..n+1] do
		tempB:=[];
		for j in [1..i-1] do
			tempB[j]:=ListM[j];
		od;
		tempB[i]:=e;
		for j in [i+1..n+1] do
			tempB[j]:=ListM[j-1];
		od;
	Add(Degen,tempB);
	od;
fi;
return Degen;
end;
############################################################################

for i in [2..number] do
	Elto:=ElementsOfSemiDirect(H,i);
	Gens:=Elto[1];	
	ComposOfGens:=Elto[2];
	n:=Length(Gens);
	tempBound:=[];
	for j in [1..n]do
		tempBound[j]:=BoundaryElement(ComposOfGens[j]);
	od;
	ImageOfGens:=[];
	for j in [1..i+1] do
		TempL:=[];
		for k in [1..n] do
		    Add(TempL,CreatElement(tempBound[k][j]));
	    od;
		Add(ImageOfGens,TempL);
	od;
	boundary:=[];
	for j in [1..i+1] do
		boundary[j]:=GroupHomomorphismByImagesNC(H,G,Gens,ImageOfGens[j]);
	od;
	Add(Boundaries,boundary);
	G:=H;
	
######## Create semidirect product ##########
	
	ConjugatorOfProd:=[];	
	for j in [1..n] do
		tempprod:=ComposOfGens[j][1];
		for k in [2..i] do
			tempprod:=tempprod*ComposOfGens[j][k];
		od;
		Add(ConjugatorOfProd,ConjugatorAutomorphismNC(M,Image(tmap,tempprod)));
	od;
	phi:=GroupHomomorphismByImagesNC(G,AutM,Gens,ConjugatorOfProd);  
	H:=SemidirectProduct(G,phi,M);
	if i<number then 
		Add(ListGroups,H);
	fi;
		
od;

####### Create degeneracy maps ############

phi:=GroupHomomorphismByImagesNC(ListGroups[1],ListGroups[2],GeneratorsOfGroup(ListGroups[1]),GeneratorsOfGroup(ListGroups[1]));
Add(Degeneracies,[phi]);
for i in [2..number] do
	G:=ListGroups[i];  ##G-->H
	H:=ListGroups[i+1];
	Elto:=ElementsOfSemiDirect(G,i-1);
	Gens:=Elto[1];	
	ComposOfGens:=Elto[2];
	n:=Length(Gens);
	tempDegen:=[];
	for j in [1..n]do
		tempDegen[j]:=DegenElement(ComposOfGens[j]);
	od;
	ImageOfGens:=[];
	for j in [1..i] do
		TempL:=[];
		for k in [1..n] do
		    Add(TempL,CreatElement(tempDegen[k][j]));
	    od;
		Add(ImageOfGens,TempL);
	od;
	degeneracy:=[];
	for j in [1..i] do
		degeneracy[j]:=GroupHomomorphismByImagesNC(G,H,Gens,ImageOfGens[j]);
	od;
	Add(Degeneracies,degeneracy);	
od;

#######################################################
GroupsList:=function(i)
return ListGroups[i+1];
end;
#######################################################
BoundariesList:=function(i,j)
	return Boundaries[i][j+1];
end;
#######################################################
DegeneraciesList:=function(i,j)
return Degeneracies[i+1][j+1];
end;
###########################################################
Nerve:=Objectify(HapSimplicialGroup,
       rec(
            groupsList:=GroupsList,
            boundariesList:=BoundariesList,
            degeneraciesList:=DegeneraciesList,
            properties:=[["length",number]]
          ));
return Nerve;
end;
###########################################################
###########################################################

NerveOfCatOneGroup_Homopre:=function(NG,NH,f,n)
	local 
	groupsG,groupsH,
	PG,PH,GensPG,Map,EmbG,ProG,EmbH,
	GenGs,FacGenGs,LGs,LHs,
	len,i,j,k,temp,x,p,
	CreateOneElementH,Mapping,
	ImgH;
	
	groupsG:=NG!.groupsList;
	groupsH:=NH!.groupsList;
	PG:=groupsG(0);
	PH:=groupsH(0);
	GensPG:=GeneratorsOfGroup(PG);
	Map:=[];
	Map[1]:=GroupHomomorphismByImages(PG,PH,GensPG,List(GensPG,x->Image(f,x)));
	Map[2]:=f;
	EmbG:=[];
	ProG:=[];
	EmbH:=[];
	GenGs:=[];
	LGs:=[];
	LHs:=[];
	for i in [2..n] do
		LGs[i]:=groupsG(i);
		LHs[i]:=groupsH(i);
		EmbG[i]:=[];
		EmbG[i][1]:=Embedding(LGs[i],1);
		EmbG[i][2]:=Embedding(LGs[i],2);
		EmbH[i]:=[];
		EmbH[i][1]:=Embedding(LHs[i],1);
		EmbH[i][2]:=Embedding(LHs[i],2);
		ProG[i]:=Projection(LGs[i]);
		GenGs[i]:=GeneratorsOfGroup(LGs[i]);
	od;
	FacGenGs:=[];
	for i in [2..n] do
		FacGenGs[i]:=[];
		len:=Length(GenGs[i]);
		for j in [1..len] do
			x:=GenGs[i][j];
			temp:=[];
			for k in [0..i-2] do
				p:=Image(EmbG[i-k][1],Image(ProG[i-k],x))^(-1)*x;
				temp[i-k]:=PreImagesRepresentative(EmbG[i-k][2],p);
				x:=Image(ProG[i-k],x);
			od;
			temp[1]:=x;
			FacGenGs[i][j]:=temp;	
		od;
	od;	
	
	##################################
	##Create one element by getting semidirect product of list [g1,g2, ...gn] 
	####################################
	CreateOneElementH:=function(L,k)  
			local i,
				m,len;
				m:=L[1];
				for i in [2..k] do    
					m:=Image(EmbH[i][1],m)*Image(EmbH[i][2],L[i]);
				od;	
			return m;
	end;
	########################################
	
    for i in [2..n] do
		len:=Length(FacGenGs[i]);
		ImgH:=[];
		for j in [1..len] do
			temp:=FacGenGs[i][j];
			for k in [1..i] do
				temp[k]:=Image(f,temp[k]);
			od;
			ImgH[j]:=CreateOneElementH(temp,i);
		od;
		Map[i+1]:=GroupHomomorphismByImages(LGs[i],LHs[i],GenGs[i],ImgH);
	od;
	######################
	Mapping:=function(k)
	 return Map[k+1];
	end;
	
return Objectify(HapSimplicialGroupMap,
       rec(
            source:=NG,
            target:=NH,
            mapping:=Mapping,
            properties:=[["length",n]]
          ));
end;	


NerveOfCatOneGroup_Homo:=function(Nf,n)
local 
    NG,NH,f, RewL;
	
    NG:=NerveOfCatOneGroup_Obj(Nf!.source,n);  ##G
	NH:=NerveOfCatOneGroup_Obj(Nf!.target,n);  ##H
	f:=Nf!.mapping;
	return NerveOfCatOneGroup_Homopre(NG,NH,f,n);
end;	
###########################################################
###########################################################

NerveOfCatOneGroup_Seq:=function(Lf,n)   

local len,i,NC,RewL;

len:=Length(Lf);
NC:=[];
for i in [1..len] do
	NC[i]:=NerveOfCatOneGroup_Obj(Lf[i]!.source,n);
od;
NC[len+1]:=NerveOfCatOneGroup_Obj(Lf[len]!.target,n);
RewL:=[];
for i in [1..len] do
	RewL[i]:=NerveOfCatOneGroup_Homopre(NC[i],NC[i+1],Lf[i]!.mapping,n);
od;
return RewL;
end;
####################################################################
####################################################################

if IsComponentObjectRep( X )  then
        if "TailMap" in NamesOfComponents( X ) and "HeadMap" in NamesOfComponents( X )  then
            C := Objectify( HapCatOneGroup, rec(
                  sourceMap := X!.TailMap,
                  targetMap := X!.HeadMap ) );
            return NerveOfCatOneGroup_Obj(C,n);
        fi;
    fi;


if IsHapCatOneGroup(X) then
	return NerveOfCatOneGroup_Obj(X,n);
fi;

if IsHapCatOneGroupHomomorphism(X) then
	return NerveOfCatOneGroup_Homo(X,n);
fi;
if IsList(X) then
	return NerveOfCatOneGroup_Seq(X,n);
fi;

end);







