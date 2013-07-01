InstallGlobalFunction(HomotopyLowerCentralSeries,
function(C,k)
    local  HomotopyOne,HomotopyTwo,D;

##############################################################
##############################################################
HomotopyOne:=function(C)
	local  n, SemidirectPro,CreateCatOneGroup,
			sC,tC,Kers,Ims,G,M,P,NatHom,Pi,PiOne,
			i,H,LG,Cat,Gens,Map,
			A,B,Gen,ProA,EmbA1,EmbA2,EmbB1,EmbB2,len,ImGen,j,g,p,m,RewL;
			
			
######################################
SemidirectPro:=function(K,H) ## Create semidirectproduct of K and H with H is normal group
	local  
	    AutH,GensK,phi;
		AutH:=AutomorphismGroup(H);
		GensK:=GeneratorsOfGroup(K);
		phi:=GroupHomomorphismByImagesNC(K,AutH,GensK,List(GensK,g->ConjugatorAutomorphism(H,g)));
		return SemidirectProduct(K,phi,H);
	end;

#######################################	
sC:=C!.sourceMap;
tC:=C!.targetMap;
M:=Kernel(sC);
P:=Image(sC);
if Order(Image(tC,M))=1 then
	NatHom:=IdentityMapping(P);
	Pi:=P;
else
	NatHom:=NaturalHomomorphismByNormalSubgroup(P,Image(tC,M));
	Pi:=NatHom!.Range;
fi;
if Order(Pi)=1 then
	return C;
fi;
PiOne:=[];
PiOne[1]:=Pi;
n:=1;
while Order(PiOne[n])>1 do
	n:=n+1;
	PiOne[n]:=CommutatorSubgroup(PiOne[n-1],Pi); 
od;

H:=[];
for i in [1..n] do
	H[i]:=PreImage(NatHom,PiOne[i]);
od;
LG:=[];
for i in [1..n] do
	LG[i]:=SemidirectPro(H[i],M);
od;
###########################################
########Create Cat One Group from semidirectProduct
CreateCatOneGroup:=function(G)
local
	Gens,ImGens,ImGent,
	Pro,Emb1,Emb2,s,t,
	g,p,m,len,i;
	Gens:=GeneratorsOfGroup(G);
	len:=Length(Gens);
	Pro:=Projection(G);
	Emb1:=Embedding(G,1);
	Emb2:=Embedding(G,2);
	ImGent:=[];
	ImGens:=[];
	for i in [1..len] do
		g:=Gens[i];
		p:=Image(Pro,g);
		m:=PreImagesRepresentative(Emb2,(Image(Emb1,p))^(-1)*g);
		ImGens[i]:=Image(Emb1,p);
		ImGent[i]:=Image(Emb1,p*Image(tC,m));
	od;
	s:=GroupHomomorphismByImages(G,G,Gens,ImGens);
	t:=GroupHomomorphismByImages(G,G,Gens,ImGent);;
	return Objectify(HapCatOneGroup,
			rec(sourceMap:=s,
			targetMap:=t));
end;
####################################################
####################################################
Cat:=[];
for i in [1..n] do
	Cat[i]:=CreateCatOneGroup(LG[i]);
od;
Map:=[]; 
for i in [2..n] do
	A:=LG[i];
	B:=LG[i-1];
	Gen:=GeneratorsOfGroup(A);
	ProA:=Projection(A);
	EmbA1:=Embedding(A,1);
	EmbA2:=Embedding(A,2);
	EmbB1:=Embedding(B,1);
	EmbB2:=Embedding(B,2);
	len:=Length(Gen);
	ImGen:=[];
	for j in [1..len] do
		g:=Gen[j];
		p:=Image(ProA,g);
		m:=PreImagesRepresentative(EmbA2,(Image(EmbA1,p))^(-1)*g);
		ImGen[j]:=Image(EmbB1,p)*Image(EmbB2,m);
	od;
	Map[i]:=GroupHomomorphismByImages(A,B,Gen,ImGen);
od;

RewL:=[];
for i in [1..n-1] do
	RewL[i]:= Objectify(HapCatOneGroupHomomorphism,
		   rec(
				source:= Cat[n-i+1],
				target:= Cat[n-i],
				mapping:= Map[n-i+1]
			  ));

od;
return RewL;
end;


##############################################################
##############################################################
HomotopyTwo:=function(C)
	local   n,
			CreateCatOneGroup,SemidirectNatHom,
			sC,tC,Kers,Ims,G,M,P,NatHom,Pi,PiTwo,
			i,H,LG,Cat,Gens,GensP,Map,
			ImM,ImgP,ImgM,x,HomP,HomM,
			A,B,Gen,ProA,EmbA1,EmbA2,EmbB1,EmbB2,len,ImGen,j,g,p,m,RewL;
	
	sC:=C!.sourceMap;
	tC:=C!.targetMap;
	M:=Kernel(sC);
	P:=Image(sC);
	GensP:=GeneratorsOfGroup(P);
	Pi:=Intersection(M,Kernel(tC));
	if Order(Pi)=1 then
		return C;
	fi;
	PiTwo:=[];
	PiTwo[1]:=Pi;
	n:=1;
	while Order(PiTwo[n])>1 do
	    n:=n+1;
		PiTwo[n]:=CommutatorSubgroup(PiTwo[n-1],Pi);
	od;
	NatHom:=[];
	ImM:=[];
	NatHom[1]:=IdentityMapping(M);
	ImM[1]:=M;
	for i in [2..n] do
		NatHom[i]:=NaturalHomomorphismByNormalSubgroup(M,PiTwo[n-i+1]);
		ImM[i]:=NatHom[i]!.Range;
	od;

	HomP:=[];
	HomM:=[];
	for i in [1..n-1] do
		Gens:=GeneratorsOfGroup(ImM[i]);
		len:=Length(Gens);
		ImgP:=[];
		ImgM:=[];
		for j in [1..len] do
			x:=PreImagesRepresentative(NatHom[i],Gens[j]);
			ImgP[j]:=Image(tC,x);
			ImgM[j]:=Image(NatHom[i+1],x);
		od;
		HomP[i]:=GroupHomomorphismByImages(ImM[i],P,Gens,ImgP);
		HomM[i]:=GroupHomomorphismByImages(ImM[i],ImM[i+1],Gens,ImgM);
	od;
	
	Gens:=GeneratorsOfGroup(ImM[n]);
	len:=Length(Gens);
	ImgP:=[];
	for j in [1..len] do
		x:=PreImagesRepresentative(NatHom[n],Gens[j]);
		ImgP[j]:=Image(tC,x);
	od;
	HomP[n]:=GroupHomomorphismByImages(ImM[n],P,Gens,ImgP);

########################################	
SemidirectNatHom:=function(k) 
	local  
		H,AutH,GensH,len,phi,ConjuAuto;
		H:=ImM[k];
		GensH:=GeneratorsOfGroup(H);
		len:=Length(GensH);
		#####################
		ConjuAuto:=function(p)
			local i,x,Img;
			Img:=[];
			for i in [1..len] do
				x:=PreImagesRepresentative(NatHom[k],GensH[i]);
				x:=p*x*p^(-1);
				Img[i]:=Image(NatHom[k],x);
			od;
			return GroupHomomorphismByImages(H,H,GensH,Img);		
		end;
		#######################
		AutH:=AutomorphismGroup(H);	 
		phi:=GroupHomomorphismByImagesNC(P,AutH,GensP,List(GensP,g->ConjuAuto(g)));
		return SemidirectProduct(P,phi,H);
	end;
########################################
LG:=[];
for i in [1..n] do
	LG[i]:=SemidirectNatHom(i);
od;
###########################################
########Create Cat One Group from semidirectProduct
CreateCatOneGroup:=function(k)
local
	G,GensG,ImGens,ImGent,
	Pro,Emb1,Emb2,s,t,
	g,p,m,len,i;
	
	G:=LG[k];
	GensG:=GeneratorsOfGroup(G);
	len:=Length(GensG);
	Pro:=Projection(G);
	Emb1:=Embedding(G,1);
	Emb2:=Embedding(G,2);
	ImGent:=[];
	ImGens:=[];
	for i in [1..len] do
		g:=GensG[i];
		p:=Image(Pro,g);
		m:=PreImagesRepresentative(Emb2,(Image(Emb1,p))^(-1)*g);
		ImGens[i]:=Image(Emb1,p);
		ImGent[i]:=Image(Emb1,p*Image(HomP[k],m));
	od;
	s:=GroupHomomorphismByImages(G,G,GensG,ImGens);
	t:=GroupHomomorphismByImages(G,G,GensG,ImGent);
	return Objectify(HapCatOneGroup,
			rec(sourceMap:=s,
			targetMap:=t));
end;
####################################################
####################################################
Cat:=[];
for i in [1..n] do
	Cat[i]:=CreateCatOneGroup(i);
od;
Map:=[]; 
for i in [1..n-1] do
	A:=LG[i];
	B:=LG[i+1];
	Gen:=GeneratorsOfGroup(A);
	ProA:=Projection(A);
	EmbA1:=Embedding(A,1);
	EmbA2:=Embedding(A,2);
	EmbB1:=Embedding(B,1);
	EmbB2:=Embedding(B,2);
	len:=Length(Gen);
	ImGen:=[];
	for j in [1..len] do
		g:=Gen[j];
		p:=Image(ProA,g);
		m:=PreImagesRepresentative(EmbA2,(Image(EmbA1,p))^(-1)*g);
		ImGen[j]:=Image(EmbB1,p)*Image(EmbB2,Image(HomM[i],m));
	od;
	Map[i]:=GroupHomomorphismByImages(A,B,Gen,ImGen);
od;

RewL:=[];
for i in [1..n-1] do
	RewL[i]:= Objectify(HapCatOneGroupHomomorphism,
		   rec(
				source:= Cat[i],
				target:= Cat[i+1],
				mapping:= Map[i]
			  ));

od;
return RewL;
end;
######################################################################

	D:=XmodToHAP(C);
    if k=1 then
		return HomotopyOne(D);
    fi;
	if k=2 then
		return HomotopyTwo(D);
	fi;
end);







