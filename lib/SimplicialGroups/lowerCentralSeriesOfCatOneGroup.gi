InstallGlobalFunction(HomotopyLowerCentralSeries,
function(C,n)
	local SemidirectPro,CreateCatOneGroup,
			sC,tC,Kers,Ims,G,M,P,NatHom,Pi,PiOne,
			i,H,LG,Cat,Gens,Map,
			A,B,Gen,ProA,EmbA1,EmbA2,EmbB1,EmbB2,len,ImGen,j,g,p,m,RewL;
			
			
			

######################################
SemidirectPro:=function(K,H) ## Create semidirectproduct of K and H with H is normal group
	local  AutH,GensK,phi;
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
NatHom:=NaturalHomomorphismByNormalSubgroup(P,Image(tC,M));
Pi:=NatHom!.Range;
PiOne:=[];
PiOne[1]:=Pi;
for i in [2..n] do
	PiOne[i]:=CommutatorSubgroup(PiOne[i-1],Pi);
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
end);

##############################################################
##############################################################
