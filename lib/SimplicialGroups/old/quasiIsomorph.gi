###################################################
##	SubQuasiIsomorph
##	QuotientQuasiIsomorph
##	QuasiIsomorph
###################################################



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


	
	
