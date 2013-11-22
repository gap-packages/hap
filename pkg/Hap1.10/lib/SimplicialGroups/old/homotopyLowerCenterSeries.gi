InstallGlobalFunction(HomotopyLowerCentralSeries,
function(C)

	local  	sC,tC,Kers,M,P,GensM,GensP,AutM,
			nOne, nTwo,
			PiOne,PiOnes,
			PiTwo,Ns,
			SemidirectPro,SemidirectNatHom,
			NatHomOne,NatHomTwos,
			Ms,Ps,HomPs,HomMs,
			ImgPs,ImgMs,x,
			CreateCatOne,CreateCatTwo,
			GOnes,GTwos,
			CatOnes,CatTwos,Gens,MapOnes,MapTwos,
			A,B,ProA,EmbA1,EmbA2,EmbB1,EmbB2,ImGens,CatMapOnes,CatMapTwos,
			i,j,g,p,m,n;

	sC:=C!.sourceMap;
	tC:=C!.targetMap;
	M:=Kernel(sC);
	P:=Image(sC);
	GensM:=GeneratorsOfGroup(M);
	GensP:=GeneratorsOfGroup(P);

	################### Pi One #####################
	################################################
	nOne:=0;
	CatMapOnes:=[];
	if Order(Image(tC,M))=1 then
		NatHomOne:=IdentityMapping(P);
		PiOne:=P;
	else
		NatHomOne:=NaturalHomomorphismByNormalSubgroup(P,Image(tC,M));
		PiOne:=NatHomOne!.Range;
	fi;
	if Order(PiOne)>1 then
		PiOnes:=[];
		PiOnes[1]:=PiOne;
		nOne:=1;
		while Order(PiOnes[nOne])>1 do
			nOne:=nOne+1;
			PiOnes[nOne]:=CommutatorSubgroup(PiOnes[nOne-1],PiOne); 
			if PiOnes[nOne-1] = PiOnes[nOne] then
				return fail;
			fi;
		od;
		PiOnes:=Reversed(PiOnes);
		Ps:=[];
		for i in [1..nOne] do
			Ps[i]:=PreImage(NatHomOne,PiOnes[i]);
		od;
		AutM:=AutomorphismGroup(M);
		######################################
		SemidirectPro:=function(k) ## Create semidirectproduct of Pk and M 
			local Pk,Gens,phi,ConjuAuto;
			
				ConjuAuto:=function(p)
					return GroupHomomorphismByImages(M,M,GensM,List(GensM,m->p^(-1)*m*p));
				end;
				Pk:=Ps[k];
				Gens:=GeneratorsOfGroup(Pk);
				phi:=GroupHomomorphismByImages(Pk,AutM,Gens,List(Gens,g->ConjuAuto(g)));
				return SemidirectProduct(Pk,phi,M);
		end;
		GOnes:=[];
		for i in [1..nOne] do
			GOnes[i]:=SemidirectPro(i);
		od;
		###########################################
		########Create Cat One Group from semidirectProduct
		CreateCatOne:=function(G)
		local
			Gens,ImGens,ImGent,
			Pro,Emb1,Emb2,s,t,
			g,p,m,n,i;
			
			Gens:=GeneratorsOfGroup(G);
			n:=Length(Gens);
			Pro:=Projection(G);
			Emb1:=Embedding(G,1);
			Emb2:=Embedding(G,2);
			ImGent:=[];
			ImGens:=[];
			for i in [1..n] do
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
		CatOnes:=[];
		for i in [1..nOne] do
			CatOnes[i]:=CreateCatOne(GOnes[i]);
		od;
		MapOnes:=[]; 
		for i in [1..nOne-1] do
			A:=GOnes[i];
			B:=GOnes[i+1];
			Gens:=GeneratorsOfGroup(A);
			ProA:=Projection(A);
			EmbA1:=Embedding(A,1);
			EmbA2:=Embedding(A,2);
			EmbB1:=Embedding(B,1);
			EmbB2:=Embedding(B,2);
			n:=Length(Gens);
			ImGens:=[];
			for j in [1..n] do
				g:=Gens[j];
				p:=Image(ProA,g);
				m:=PreImagesRepresentative(EmbA2,(Image(EmbA1,p))^(-1)*g);
				ImGens[j]:=Image(EmbB1,p)*Image(EmbB2,m);
			od;
			MapOnes[i]:=GroupHomomorphismByImages(A,B,Gens,ImGens);
		od;

		for i in [1..nOne-1] do
			CatMapOnes[i]:= Objectify(HapCatOneGroupHomomorphism,
				   rec(
						source:= CatOnes[i],
						target:= CatOnes[i+1],
						mapping:= MapOnes[i]
					  ));

		od;
	fi;

################### Pi Two #####################
################################################

	CatMapTwos:=[];
	PiTwo:=Intersection(M,Kernel(tC));
	if Order(PiTwo)>1 then
		#######M/N -->P  
		Ns:=[];
		Ns[1]:=PiTwo;
		nTwo:=1;
		while Order(Ns[nTwo])>1 do
			nTwo:=nTwo+1;
			Ns[nTwo]:=CommutatorSubgroup(Ns[nTwo-1],P);
			if Ns[nTwo-1] = Ns[nTwo] then
				return fail;
			fi;
		od;
		Ns:=Reversed(Ns);
		NatHomTwos:=[];
		NatHomTwos[1]:=IdentityMapping(M);
		Ms:=[];
		Ms[1]:=M;
		for i in [2..nTwo] do
			NatHomTwos[i]:=NaturalHomomorphismByNormalSubgroup(M,Ns[i]);
			Ms[i]:=NatHomTwos[i]!.Range;
		od;

		HomPs:=[];  ###### M/Ni-->P
		HomMs:=[];  ###### M/Ni-1 -->M/Ni
		for i in [1..nTwo-1] do
			Gens:=GeneratorsOfGroup(Ms[i]);
			n:=Length(Gens);
			ImgPs:=[];
			ImgMs:=[];
			for j in [1..n] do
				x:=PreImagesRepresentative(NatHomTwos[i],Gens[j]);
				ImgPs[j]:=Image(tC,x);
				ImgMs[j]:=Image(NatHomTwos[i+1],x);
			od;
			HomPs[i]:=GroupHomomorphismByImages(Ms[i],P,Gens,ImgPs);
			HomMs[i]:=GroupHomomorphismByImages(Ms[i],Ms[i+1],Gens,ImgMs);
		od;
		######### for M/Nn-->P
		Gens:=GeneratorsOfGroup(Ms[nTwo]);
		n:=Length(Gens);
		ImgPs:=[];
		for j in [1..n] do
			x:=PreImagesRepresentative(NatHomTwos[nTwo],Gens[j]);
			ImgPs[j]:=Image(tC,x);
		od;
		HomPs[nTwo]:=GroupHomomorphismByImages(Ms[nTwo],P,Gens,ImgPs);
		########################################	
		####Semidirect P and M/Nk
		SemidirectNatHom:=function(k) 
			local  Mk,AutMk,GensMk,n,phi,ConjuAuto;
				
				Mk:=Ms[k];
				GensMk:=GeneratorsOfGroup(Mk);
				n:=Length(GensMk);
				#####################
				ConjuAuto:=function(p)
					local i,x,Img;
					Img:=[];
					for i in [1..n] do
						x:=PreImagesRepresentative(NatHomTwos[k],GensMk[i]);
						x:=p^(-1)*x*p;
						Img[i]:=Image(NatHomTwos[k],x);
					od;
					return GroupHomomorphismByImages(Mk,Mk,GensMk,Img);
				end;
				#######################
				AutMk:=AutomorphismGroup(Mk);	 
				phi:=GroupHomomorphismByImages(P,AutMk,GensP,List(GensP,g->ConjuAuto(g)));
				return SemidirectProduct(P,phi,Mk);
			end;
		########################################
		if nOne>1 then
			GTwos:=[GOnes[nOne]];
			for i in [2..nTwo] do
				GTwos[i]:=SemidirectNatHom(i);
			od;
		else
			GTwos:=[];
			for i in [1..nTwo] do
				GTwos[i]:=SemidirectNatHom(i);
			od;
		fi;
			
		###########################################
		########Create Cats One Group from semidirectProduct
		CreateCatTwo:=function(k)
		local
			G,GensG,ImGens,ImGent,
			Pro,Emb1,Emb2,s,t,
			g,p,m,n,i;
			
			G:=GTwos[k];
			GensG:=GeneratorsOfGroup(G);
			n:=Length(GensG);
			Pro:=Projection(G);
			Emb1:=Embedding(G,1);
			Emb2:=Embedding(G,2);
			ImGent:=[];
			ImGens:=[];
			for i in [1..n] do
				g:=GensG[i];
				p:=Image(Pro,g);
				m:=PreImagesRepresentative(Emb2,(Image(Emb1,p))^(-1)*g);
				ImGens[i]:=Image(Emb1,p);
				ImGent[i]:=Image(Emb1,p*Image(HomPs[k],m));
			od;
			s:=GroupHomomorphismByImages(G,G,GensG,ImGens);
			t:=GroupHomomorphismByImages(G,G,GensG,ImGent);
			return Objectify(HapCatOneGroup,
					rec(sourceMap:=s,
					targetMap:=t));
		end;
		####################################################
		####################################################
		if nOne>1 then
			CatTwos:=[CatOnes[nOne]];
			for i in [2..nTwo] do
				CatTwos[i]:=CreateCatTwo(i);
			od;		
		else
			CatTwos:=[];
			for i in [1..nTwo] do
				CatTwos[i]:=CreateCatTwo(i);
			od;
		fi;
		
		MapTwos:=[]; 
		for i in [1..nTwo-1] do
			A:=GTwos[i];
			B:=GTwos[i+1];
			Gens:=GeneratorsOfGroup(A);
			ProA:=Projection(A);
			EmbA1:=Embedding(A,1);
			EmbA2:=Embedding(A,2);
			EmbB1:=Embedding(B,1);
			EmbB2:=Embedding(B,2);
			n:=Length(Gens);
			ImGens:=[];
			for j in [1..n] do
				g:=Gens[j];
				p:=Image(ProA,g);
				m:=PreImagesRepresentative(EmbA2,(Image(EmbA1,p))^(-1)*g);
				ImGens[j]:=Image(EmbB1,p)*Image(EmbB2,Image(HomMs[i],m));
			od;
			MapTwos[i]:=GroupHomomorphismByImages(A,B,Gens,ImGens);
		od;
		for i in [1..nTwo-1] do
			CatMapTwos[i]:= Objectify(HapCatOneGroupHomomorphism,
				   rec(
						source:= CatTwos[i],
						target:= CatTwos[i+1],
						mapping:= MapTwos[i]
					  ));

		od;
	fi;

return Concatenation(CatMapOnes,CatMapTwos);
end);

######################################################################
######################################################################
InstallGlobalFunction(PersistentHomologyOfCatOneGroup,
function(arg)
local  C,n,p,Maps,
	   PrimeOne,PrimeTwo,PrimeOneTwo,
	   LenOne,LenTwo,LenOneTwo;
	   
	C:=arg[1];
	n:=arg[2];
	PrimeOne:=PrimeDivisors(Size(HomotopyGroup(C,1)));
	PrimeTwo:=PrimeDivisors(Size(HomotopyGroup(C,2)));
	LenOne:=Length(PrimeOne);
	LenTwo:=Length(PrimeTwo);
	if LenOne >1 or LenTwo >1 then
		Print("Homotopy groups have not order of prime power");
		return;
	fi;
	PrimeOneTwo:=Set(Concatenation(PrimeOne,PrimeTwo));
	LenOneTwo:=Length(PrimeOneTwo);
	if LenOneTwo=0 then
		Maps:=HomotopyLowerCentralSeries(C);
		Maps:=NerveOfCatOneGroup(C,n+1);
		Maps:=ChainComplexOfSimplicialGroup(Maps);
		if Length(arg)=3 then
			p:=arg[3];
			Maps:=TensorWithIntegersModP(Maps,p);
		fi;
		return [Homology(Maps,n)];
	fi;	
	
	if LenOneTwo =1 then
		p:=PrimeOneTwo[1];	
		if Length(arg)=3 then
			p:=arg[3];
		fi;
		Maps:=HomotopyLowerCentralSeries(C);
		Maps:=NerveOfCatOneGroup(Maps,n+1);
		Maps:=ChainComplexOfSimplicialGroup(Maps);
		Maps:=List(Maps,f->TensorWithIntegersModP(f,p));
		Maps:=List(Maps,f->HomologyVectorSpace(f,n));
		return LinearHomomorphismsPersistenceMat(Maps);
	fi;		
	
	if LenOneTwo=2 then
		Print("Prime divisors of homotopy groups are not the same");
		return;
	fi;
end);











