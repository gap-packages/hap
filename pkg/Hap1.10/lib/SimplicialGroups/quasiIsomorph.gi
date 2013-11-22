#############################################################################
#0
#F	SubQuasiIsomorph
##	Input:	A finite cat-1-group C
##	Output: A quasi-isomorphic sub cat-1-group of C
##
InstallGlobalFunction(SubQuasiIsomorph,function(C)
local  
	s,t,G,H,Kers,Kert,Kersnt,tKers,OrdPiOne,OrdPiTwo,OrdPi,
	LS,Lx,x,sx,Ordsx,flag,
    newGens,news,newt;
	
	s:= C!.sourceMap;
    t:= C!.targetMap;
    G:=Source(s);
    Kers:=Kernel(s);
	Kert:=Kernel(t);
	Kersnt:=Intersection(Kers,Kert);
	tKers:=Image(t,Kers);
	OrdPiOne:=Order(HomotopyGroup(C,1));
	OrdPiTwo:=Order(HomotopyGroup(C,2));
	OrdPi:=OrdPiOne*OrdPiTwo;
	LS:=ConjugacyClassesSubgroups(LatticeSubgroups(G));
	if not IsMutable(LS) then
	    LS:= ShallowCopy(LS);
	fi;
	Sort(LS,function(x,y) return Size(x[1])<Size(y[1]); end);
	flag:=0;
	for Lx in LS do
	    x:=Lx[1];
		if Order(x)>= OrdPi then
	    if IsSubgroup(x,Kersnt) then
			for x in Lx do
				if  IsSubgroup(x,Image(s,x)) then
				if  IsSubgroup(x,Image(t,x)) then
					sx:=Image(s,x);
					Ordsx:=Order(sx);
					if  Ordsx = Order(Image(t,Intersection(Kers,x)))*OrdPiOne then
					if  Ordsx = Order(Intersection(sx,tKers))*OrdPiOne then
						H:=x;
						flag:=1;
						break;
					fi;
					fi;
				fi;
				fi;
			od;
		fi;
		fi;
		if flag =1 then
			break;
		fi;
	od;
	if H=G then 
		return C;
	fi;
    newGens:=GeneratorsOfGroup(H);
    news:=GroupHomomorphismByImagesNC(H,H,newGens,List(newGens,x->Image(s,x)));
    newt:=GroupHomomorphismByImagesNC(H,H,newGens,List(newGens,x->Image(t,x)));
    return Objectify(HapCatOneGroup,rec( 
				sourceMap:=news,
				targetMap:=newt));
end);
##
#################### end of SubQuasiIsomorph ################################

#############################################################################
#0
#F	QuotientQuasiIsomorph
##	Input:	A finite cat-1-group C
##	Output:	A quasi-isomorphic quotient cat-1-group of C
##
InstallGlobalFunction(QuotientQuasiIsomorph,function (C)
local 	
	s,t,G,H,Kers,Kert,Kersnt,Ims,OrdIms,Imt,OrdPiOne,OrdPiTwo,Ord,
	LN,x,n,i,
	OrderPiOneGx,OrderPiTwoGx,
	epi,newG,newGens,news,newt;

	s:=C!.sourceMap;
    t:=C!.targetMap;
    G:=Source(s);
    Kers:=Kernel(s);
	Ims:=Image(s);
	OrdIms:=Order(Ims);
	Imt:=Image(t);
	Kert:=Kernel(t);
	Kersnt:=Intersection(Kers,Kert);
	OrdPiOne:= Order(HomotopyGroup(C,1));
	OrdPiTwo:= Order(HomotopyGroup(C,2));
	Ord:=Order(G)/(OrdPiOne*OrdPiTwo);
	
    ######################################################################
	#1
	OrderPiOneGx:=function(x)
	local tsx;
		
	    tsx:=Image(t,PreImages(s,Intersection(Ims,x)));
		return (OrdIms*Order(Intersection(tsx,x)))/
				(Order(Intersection(Ims,x))*Order(tsx));
	end;
    ##	
	######################################################################
	#1
	OrderPiTwoGx:=function(x)
	local f;
		
	    f:=NaturalHomomorphismByNormalSubgroup(G,x);
	    return Order(Intersection(Image(f,PreImages(s,Intersection(Ims,x))),
				Image(f,PreImages(t,Intersection(Imt,x)))));
	end;
	##
	######################################################################
	
	LN:=NormalSubgroups(G);
	if not IsMutable(LN) then
	    LN:= ShallowCopy(LN);
	fi;
	Sort(LN,function(x,y) return Size(x)>Size(y); end);
	n:=Length(LN);
	for i in [1..n] do
		x:=LN[i];
		if Order(x) <= Ord then
		if IsSubgroup(x,Image(s,x)) then
		if IsSubgroup(x,Image(t,x)) then
		if IsSubgroup(x,CommutatorSubgroup(PreImages(s,Intersection(Ims,x)),
				PreImages(t,Intersection(Imt,x)))) then
		if Order(Kersnt) = Order(Intersection(Kersnt,x))*OrdPiTwo then
		if OrderPiTwoGx(x) = OrdPiTwo then
		if OrderPiOneGx(x) = OrdPiOne then
			H:=x;
			break;
		fi;
		fi;
		fi;
		fi;
		fi;
		fi;
		fi;
	od;
	if Order(H)=1 then
		return C;
	fi;
	
	epi:=NaturalHomomorphismByNormalSubgroup(G,H);
    newG :=Image(epi);
    newGens:=GeneratorsOfGroup(newG);
    news:=GroupHomomorphismByImagesNC(newG,newG,newGens,
		List(newGens,x->Image(epi,Image(s,PreImagesRepresentative(epi,x)))));
    newt:=GroupHomomorphismByImagesNC(newG,newG,newGens,
		List(newGens,x->Image(epi,Image(t,PreImagesRepresentative(epi,x)))));
    return Objectify(HapCatOneGroup,rec(
          sourceMap:=news,
          targetMap:=newt));
end);
##
#################### end of QuotientQuasiIsomorph ###########################

#############################################################################
#0
#F	QuasiIsomorph
##	Input:	A finite cat-1-group or a finite crossed module X
##	Output:	A quasi-isomorphism of X
##
InstallGlobalFunction(QuasiIsomorph,function (X)
local QuasiIsomorphOfCat, QuasiIsomorphOfCross;
		
	######################################################################	  
	#1
    #F	QuasiIsomorphOfCat
	##	Input:	A finite cat-1-group C
	##	Output:	A quasi-isomorphism of C
	##
	QuasiIsomorphOfCat:=function(C)
	local D;
	
		D:=QuotientQuasiIsomorph(C);
		D:=SubQuasiIsomorph(D);
		while Size(D) < Size(C)  do
			C:=D;
			D:=QuotientQuasiIsomorph(C);
			if Size(D) < Size(C) then
				D:=SubQuasiIsomorph(D);
			fi;
		od;
		return D;
	end;
	##
	############### end of QuasiIsomorphOfCat ############################

	######################################################################	
	#1
    #F	QuasiIsomorphOfCross
	##	Input:	A finite crossed module XC
	##	Output:	A quasi-isomorphism of XC
	##
	QuasiIsomorphOfCross:=function(XC)
		local C,D;
		
		C:=CatOneGroupByCrossedModule(XC);
		D:=QuasiIsomorphOfCat(C);
		return CrossedModuleByCatOneGroup(D);
	end;
	##
	############### end of QuasiIsomorphOfCross ##########################
	
	if IsHapCatOneGroup(X) then
		return QuasiIsomorphOfCat(X);
	fi;
	if  IsHapCrossedModule(X) then
		return QuasiIsomorphOfCross(X);
	fi;
end);
##	
#################### end of QuasiIsomorph ###################################

#############################################################################
#0
#F	IsomorphismCatOneGroups	
##	Input: Two cat-1-group C and D
##	Output: return a morphism between C and D if C and D are isomorphic.
##			return false for other else
##
InstallGlobalFunction(IsomorphismCatOneGroups, function(C,D)
local 
	sC,tC,G,sD,tD,GD,
	xC,xD,A,attr,actAttr,M,p,f,h,
	ActToMap,ActToPair,ActToSubgroup,FindOrbit,processOrbit;
	
	######################################################################
	#1
	ActToMap:=function(s,f)
	  return InverseGeneralMapping(f)*s*f;
	end;
	##
	######################################################################
	#1
	ActToPair:=function(p,f)
		local h;
		h:=InverseGeneralMapping(f);
		return [h*p[1]*f,h*p[2]*f];
	end;
	##
	######################################################################
	#1
	ActToSubgroup:=function(K,f) 
		return Image(f,K); 
	end;
	##
	######################################################################
	#1
	#F	FindOrbit
	##
	FindOrbit:=function(A,attr,actAttr)
	local 
		Gens,Orb,T,Dict,S,op,qs,g,p,img,h;

		Gens:=SmallGeneratingSet(A);
        Orb:=[attr(xC)];
		T:=[One(A)];
		Dict:=NewDictionary(Orb[1],true);
		AddDictionary(Dict,Orb[1],1);
		S:=TrivialSubgroup(A);
		op:=1;
		qs:=Size(A);
		while op<=Length(Orb) and Size(S)<qs do
			for g in Gens do
				img:=actAttr(Orb[op],g);
				p:=LookupDictionary(Dict,img);
				if p=fail then
					Add(Orb,img);
					AddDictionary(Dict,img,Length(Orb));
					Add(T,T[op]*g);
					qs:=Size(A)/Length(Orb);
				elif Size(S)<=qs/2 then # otherwise stabilizer cant grow
					h:=T[op]*g/T[p];
					S:=ClosureSubgroup(S,h);
				fi;
			od;
			op:=op+1;
		od;
		return [Dict,S,T];
	end;
	##
	######################################################################
	##
	processOrbit:=function()
	
		M:=FindOrbit(A,attr,actAttr);
		p:=LookupDictionary(M[1],attr(xD));
		if p=fail then
			return fail;
		fi;
		A:=M[2];
		if p<>1 then
			h:=M[3][p]^-1;
			f:=f*h;
			xD:=ActToPair(xD,h);
		fi;
	end;
	##
	######################################################################
	
	######################################################################
	sC:=C!.sourceMap;
	tC:=C!.targetMap;
	G:=Source(sC);
	sD:=D!.sourceMap;
	tD:=D!.targetMap;
	GD:=Source(sD);
	
	f:=IsomorphismGroups(GD,G);
	if f = fail then
		return fail;
	fi;
	if IsomorphismGroups(HomotopyGroup(C,1),HomotopyGroup(D,1))=fail then
		return fail;
	fi;
	if IsomorphismGroups(HomotopyGroup(C,2),HomotopyGroup(D,2))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Kernel(tC),Kernel(tD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Kernel(sC),Kernel(sD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Image(tC),Image(tD))=fail then
		return fail;
	fi;
	if IsomorphismGroups(Image(sC),Image(sD))=fail then
		return fail;
	fi;

	xC:=[sC,tC]; 
	xD:=ActToPair([sD,tD],f);
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;
	
	############## Image of first component #####################
	A:=AutomorphismGroup(G);
	attr:=s->Image(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	
	############## Kernel of first component #####################
	attr:=s->Kernel(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;
	
	############## Image of second component #####################
	attr:=s->Image(s[2]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;	
	############## Kernel of second component #####################
	attr:=s->Kernel(s[1]);
	actAttr:=ActToSubgroup;
	processOrbit();
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;	
	############## First component #####################
	attr:=s->s[1];
	actAttr:=ActToMap;
	processOrbit();
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;	
	############## Second component #####################
	attr:=s->s[2];
	actAttr:=ActToMap;
	processOrbit();
	if xC=xD then
		return InverseGeneralMapping(f);
	fi;
	return fail;
end);
##
##################### IsomorphismCatOneGroups ###############################

#############################################################################
#0
#F	XmodToHAP
##	Input: 	A cat-1-group from the Xmod package 
##	Output:	The cat-1-group with data type from the HAP package 
##
InstallGlobalFunction(XmodToHAP,function (X)
local G,Gens,s,t;

    if IsHapCatOneGroup(X)  then
        return X;
    fi;
    if IsComponentObjectRep(X)  then
        if "TailMap" in NamesOfComponents(X) and 
					"HeadMap" in NamesOfComponents(X) then
		    G:=X!.Source;
			Gens:=GeneratorsOfGroup(G);
			s:=X!.TailMap;
			t:=X!.HeadMap;
            return Objectify(HapCatOneGroup,rec(
					sourceMap:=GroupHomomorphismByImagesNC(G,G,Gens,List(Gens,
							g->Image(s,g))),
					targetMap:=GroupHomomorphismByImagesNC(G,G,Gens,List(Gens,
							g->Image(t,g)))
					));      
        fi;
    fi;
    Print("Argument is not a cat-1-group from the Xmod package.\n");
    return fail;
end);
##
#################### end of XmodToHAP #######################################

	
	
