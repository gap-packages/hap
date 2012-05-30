#C Ellis & Le

InstallOtherMethod( GroupHomomorphismByImagesNC, "for group with no generators",
    [IsGroup,IsGroup,IsEmpty,IsEmpty], SUM_FLAGS,
        function(g,h,gg,gh)
    return GroupHomomorphismByFunction(g,h,x->One(h),false,x->One(g));
end);

###########################################################
###########################################################
InstallGlobalFunction(NerveOfCatOneGroup,
function(CC,number)
local 
    ListGroups,Boundaries, Degeneracies, Nerve,
    smap,tmap,e,
	C,G,M,H,
	AutM,phi,
	g,m,tempprod,ConjugatorOfProd,tempBound,tempDegen,
	ComposOfGens,Gens,TempB,ImageOfGens,TempL,
	Elto,boundary,degeneracy,
	i,j,k,n,
	ElementsOfSemiDirect,BoundaryElement,CreatElement,DegenElement,
	GroupsList,BoundariesList,DegeneraciesList;

C:=XmodToHAP(CC);  
	
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
	Bound:=[[Image(tmap,ListM[1])*ListM[2]],[Image(tmap,ListM[1])*ListM[2]*Image(tmap,ListM[1]^(-1))*ListM[1]],[ListM[1]]];	
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
		tempB[i-1]:=Image(tmap,ListM[i-1])*ListM[i]*Image(tmap,ListM[i-1]^(-1))*ListM[i-1];
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
end);

###########################################################

