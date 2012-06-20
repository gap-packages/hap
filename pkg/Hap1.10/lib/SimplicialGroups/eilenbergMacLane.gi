InstallGlobalFunction(EilenbergMacLaneSimplicialGroup,
function(G,NN,number)
	local 
		nn,i,j,n,k,pos,len1,len2,temp,x,y,
		CoF,CoD,MN,T,N,H1,H2,GensH2,lengensH2,Prj,
		ListGensH2,GensH1,lengensH1,ListGensH1,Emb,
		NumberFace,NumberDegener,ImageGens,
		IdTri,GensG,ListId,TriGroup,
		HomoTriToTri,HomoGToTri,HomoTriToG,HomoIdG,
		Surjection,LenOfSur,
		Faces,Degeners,ListGroups,
		GroupsList,FaciesList,DegeneraciesList,
		CreateSurjective,FindComposite,CoDegeneracies,CoFaces,SearchPos;
 
#if NN<3 and IsPcpGroup(G) then
#return 
#EilenbergMacLaneSimplicialGroup_alt(G,NN,number);
#fi;


nn:=NN-1;
##############Find all surjection from [m] to [n]#############
##############Output: List of list couples [i,j]##############
CreateSurjective:=function(m,n)
local
	CreatCouple,Rew,y,M;
	
	if m<n then
		return fail;
	fi;
	if m=0 then
		return [[[0,0]]];
	fi;
##############################################	
	CreatCouple:=function(m,n)
	local 
		A,B,M,x,i;
		if n=0 then
			M:=[];
			for i in [0..m] do
				Add(M,[i,0]);
			od;
			return [M];
		fi;
		if m=n-1 then
			M:=[];
			for i in [0..m] do
				Add(M,[i,i]);
			od;
			return [M];
		fi;
		A:=CreatCouple(m-1,n-1);
			for x in A do
				Add(x,[m,n-1]);
			od;
		
		B:=CreatCouple(m-1,n);
			for x in B do
				Add(x,[m,n]);
			od;
		
		return Concatenation(A,B);		
	end;
##################################################
	Rew:=CreatCouple(m-1,n);
	for y in Rew do
		Add(y,[m,n]);
	od;
	return Rew;
end;

###########Create CoFaces:[n]-->[n+1]########################
###########Outout: List of list od couples ##################
CoFaces:=function(n)
local i,k,M,Rew;
	Rew:=[];
	for i in [0..n+1] do
		M:=[];
		for k in [0..i-1] do
			Add(M,[k,k]);
		od;
		for k in [i..n] do
			Add(M,[k,k+1]);
		od;
		Add(Rew,M);
	od;
	return Rew;
end;

#############Create CoDegeneracies:[n]-->[n-1]################
###########Output: List of list od couples ##################
CoDegeneracies:=function(n)
local i,k,M,Rew;
	Rew:=[];
	for i in [0..n-1] do
		M:=[];
		for k in [0..i-1] do
			Add(M,[k,k]);
		od;
		Add(M,[i,i]);
		for k in [i+1..n] do
			Add(M,[k,k-1]);
		od;
		Add(Rew,M);
	od;
	return Rew;
end;

###############Find the composite of [m]->[n]->>[k] ##################
###############Output:[m]->>[k] if it exist or 0 if not exist#########
FindComposite:=function(M,N)
local Rew,k,m,Temp,i,x,y;
	k:=nn;
	m:=Length(M)-1;
	x:=M[m+1][2];
	y:=N[x+1][2];
	if y<>k then
		return 0;
	fi;
	Rew:=[];
	Temp:=[];
	for i in [0..m] do
		x:=M[i+1][2];
		y:=N[x+1][2];
		Add(Rew,[i,y]);
		Add(Temp,y);
	od;
	if Length(Set(Temp))<k+1 then
		return 0;
	fi;
	return Rew;
end;

############Find the position of M in List L##################
############If yes return pos, no return 0####################
SearchPos:=function(M,L)
	local
		i,n;
	n:=Length(L);
	for i in [1..n] do	
		if M=L[i] then
			return i;
		fi;
	od;
	return 0;
end;

######################################################
######################################################			
	Surjection:=[];
	LenOfSur:=[];
	for i in [nn..number] do
		Surjection[i+1]:=CreateSurjective(i,nn); 			###[i+1]
		LenOfSur[i+1]:=Length(Surjection[i+1]);      ##[i+1]
	od;
	TriGroup:=Group(Identity(G));
	ListGroups:=[];
	for i in [0..nn-1] do
		ListGroups[i+1]:=TriGroup;               
	od;
	ListGroups[nn+1]:=G;
	for i in [nn+1..number] do
		ListGroups[i+1]:=DirectProduct(List([1..LenOfSur[i+1]],x->G));  ##[i+1]
	od;
	HomoTriToTri:=GroupHomomorphismByImages(TriGroup,TriGroup,[],[]);
	IdTri:=Identity(TriGroup);
	GensG:=GeneratorsOfGroup(G);
	ListId:=[];
	for i in [1..Length(GensG)] do
		ListId[i]:=IdTri;
	od;
	HomoGToTri:=GroupHomomorphismByImages(G,TriGroup,GensG,ListId);
	HomoTriToG:=GroupHomomorphismByImages(TriGroup,G,[],[]);
	HomoIdG:=GroupHomomorphismByImages(G,G,GensG,GensG);
	
	##########################################
	########Create Faces #####################
	Faces:=[];             ####[i][j+1]
	for i in [1..number] do
		Faces[i]:=[];
	od;
	for i in [1..nn-1] do
		for j in [0..i] do
			Faces[i][j+1]:=HomoTriToTri;
		od;
	od;
	##########Create Faces at position nn-->####
	if nn>0 then
		for j in [0..nn] do
			Faces[nn][j+1]:=HomoGToTri;
		od;
	fi;
    ############### [n-1]->[n]->>[nn] with n,i([n-1] ->>[n]),k (pos [n]->>[nn]find pos in [n-1]->>[nn]
	############### Output: pos if exist or 0 if not exist###################################
	NumberFace:=[];
	for n in [1..number] do
		NumberFace[n]:=[];
		for i in [0..n] do
			NumberFace[n][i+1]:=[];
		od;
	od;
	
	for n in [nn+1..number] do
		CoF:=CoFaces(n-1);
		MN:=Surjection[n];
		for i in [0..n] do
			for k in [1..LenOfSur[n+1]] do
				N:=Surjection[n+1][k];
				T:=FindComposite(CoF[i+1],N);
				if T=0 then
					NumberFace[n][i+1][k]:=0;
				else
					NumberFace[n][i+1][k]:=SearchPos(T,MN);  ##i+1
				fi;
			od;
		od;
	od;
	
   ###########################################
	for n in [nn+1..number] do
		H2:=ListGroups[n+1];
		H1:=ListGroups[n];
		len2:=LenOfSur[n+1];
		len1:=LenOfSur[n];
		GensH2:=GeneratorsOfGroup(H2);
		lengensH2:=Length(GensH2);
		Prj:=[];
		for k in [1..len2] do
			Prj[k]:=Projection(H2,k);
		od;
		ListGensH2:=[];
		for k in [1..lengensH2] do
			x:=GensH2[k];
			temp:=[];
			for j in [1..len2] do
				temp[j]:=Image(Prj[j],x);
			od;
			Add(ListGensH2,temp);
		od;
		######
		Emb:=[];
		if n=nn+1 then      ## at position nn, there is only G
				Emb[1]:=HomoIdG;
		else
			for k in [1..len1] do
				Emb[k]:=Embedding(H1,k);
			od;
		fi;
	########Find Faces n i##########################
		for i in [0..n] do
			ImageGens:=[];
			for j in [1..lengensH2] do	
				x:=Identity(H1);
				for k in [1..len2] do
					pos:=NumberFace[n][i+1][k];
					if pos<>0 then
						x:=x*Image(Emb[pos],ListGensH2[j][k]);
					fi;
				od;
			    ImageGens[j]:=x;
			od; 
		    Faces[n][i+1]:=GroupHomomorphismByImages(H2,H1,GensH2,ImageGens);
		od;
	od;							
	################################################
	########Create Degeneracies#####################
	Degeners:=[];             ####[i+1][j+1]
	for i in [0..number-1] do
		Degeners[i+1]:=[];
	od;
	for i in [0..nn-2] do
		for j in [0..i] do
			Degeners[i+1][j+1]:=HomoTriToTri;
		od;
	od;
	########Create Degeneracies at position nn-1#####################
	if nn>1 then
		for j in [0..nn-1] do
			Degeners[nn][j+1]:=HomoTriToG;
		od;
	fi;
	
	############### [n+1]->[n]->>[nn] with n,i([n+1]->>[n]),k (pos [n]->>[nn]find pos in [n+1]->>[nn]
	############### Output: pos if exist or 0 if not exist###################################
	NumberDegener:=[];  ####[n+1][i+1][k]
	for n in [0..number-1] do
		NumberDegener[n+1]:=[];       
		for i in [0..n] do
			NumberDegener[n+1][i+1]:=[];
		od;
	od;
	
	for n in [nn..number-1] do
		CoD:=CoDegeneracies(n+1);
		MN:=Surjection[n+2];
		for i in [0..n] do
			for k in [1..LenOfSur[n+1]] do
				N:=Surjection[n+1][k];
				T:=FindComposite(CoD[i+1],N);
				if T=0 then
					NumberDegener[n+1][i+1][k]:=0;
				else
					NumberDegener[n+1][i+1][k]:=SearchPos(T,MN);  ##i+1
				fi;
			od;
		od;
	od;
   ###########################################
	for n in [nn..number-1] do
		H1:=ListGroups[n+1];
		H2:=ListGroups[n+2];
		len1:=LenOfSur[n+1];
		len2:=LenOfSur[n+2];
		GensH1:=GeneratorsOfGroup(H1);
		lengensH1:=Length(GensH1);
		Prj:=[];
		if n=nn then   ## at position nn, there is only G
				Prj[1]:=HomoIdG;
		else
			for k in [1..len1] do
				Prj[k]:=Projection(H1,k);
			od;
		fi;
		ListGensH1:=[];
		for k in [1..lengensH1] do
			x:=GensH1[k];
			temp:=[];
			for j in [1..len1] do
				temp[j]:=Image(Prj[j],x);
			od;
			Add(ListGensH1,temp);
		od;
		######
		Emb:=[];	
		for k in [1..len2] do
			Emb[k]:=Embedding(H2,k);
		od;
		
		########Find Faces n i##########################
		for i in [0..n] do
			ImageGens:=[];
			for j in [1..lengensH1] do	
				x:=Identity(H2);
				for k in [1..len1] do
					pos:=NumberDegener[n+1][i+1][k];
					if pos<>0 then
						x:=x*Image(Emb[pos],ListGensH1[j][k]);
					fi;
				od;
			    ImageGens[j]:=x;
			od; 
		    Degeners[n+1][i+1]:=GroupHomomorphismByImages(H1,H2,GensH1,ImageGens);
		od;
    od;						
		
#######################################################
GroupsList:=function(i)
	return ListGroups[i+1];
end;
#######################################################
FaciesList:=function(i,j)
	return Faces[i][j+1];
end;
#######################################################
DegeneraciesList:=function(i,j)
	return Degeners[i+1][j+1];
end;
###########################################################
return Objectify(HapSimplicialGroup,
       rec(
            groupsList:=GroupsList,
            boundariesList:=FaciesList,
            degeneraciesList:=DegeneraciesList,
            properties:=[["length",number]]
          ));
end);		



###############################################################
###############################################################

InstallGlobalFunction(EilenbergMacLaneSimplicialGroupMap,
function(f,n,len)
	local 
		H,G,KH,KG,
		Mapf,Mapping,HH,GG,
		GensH,Prj,Emb,ImgsH,
		i,j,k,t,lengens,
		h,g;
		
H:=f!.Source;
G:=f!.Range;
KH:=EilenbergMacLaneSimplicialGroup(H,n,len);
KG:=EilenbergMacLaneSimplicialGroup(G,n,len);
Mapf:=[];
for i in [0..n-2] do
	Mapf[i+1]:=GroupHomomorphismByImages(Group(Identity(H)),Group(Identity(G)),[],[]);
od;
Mapf[n]:=f;  ##n-1
for i in [n..len] do
	HH:=KH!.groupsList(i);
	GG:=KG!.groupsList(i);
	GensH:=GeneratorsOfGroup(HH);
	Prj:=[];
	Emb:=[];
	k:=Length(HH!.DirectProductInfo!.groups);	
	for j in [1..k] do
		Prj[j]:=Projection(HH,j);
		Emb[j]:=Embedding(GG,j);
	od;
	ImgsH:=[];
	lengens:=Length(GensH);
    for j in [1..lengens] do
		h:=GensH[j];
		g:=Identity(GG);
		for t in [1..k] do
			g:=g*Image(Emb[t],Image(f,Image(Prj[t],h)));
		od;
		ImgsH[j]:=g;
	od;
	Mapf[i+1]:=GroupHomomorphismByImages(HH,GG,GensH,ImgsH);
od;
###################
Mapping:=function(i)
	return Mapf[i+1];
end;
###################
return Objectify(HapSimplicialGroupMap,
       rec(
            source:=KH,
            target:=KG,
            mapping:=Mapping,
            properties:=[["length",len]]
          ));
end);		
