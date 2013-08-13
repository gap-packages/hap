 InstallGlobalFunction(SubQuasiIsomorph_alt,
    function (C)
    local  s,t,G,H,Kers,Kert,Kersnt,tKers,OrdPiOne,OrdPiTwo,OrdPi,
	LS,Lx,x,sx,Ordsx,flag,
    newGens,news,newt,D;
	
	s:= C!.sourceMap;
    t:= C!.targetMap;
    G:=Source( s );
    Kers:=Kernel( s );
	Kert:=Kernel( t );
	Kersnt:=Intersection( Kers , Kert );
	tKers:=Image( t , Kers );
	OrdPiOne:=Order( HomotopyGroup( C , 1 ) );
	OrdPiTwo:=Order( HomotopyGroup( C , 2 ) );
	OrdPi:=OrdPiOne*OrdPiTwo;
	LS:=ConjugacyClassesSubgroups(LatticeSubgroups( G ));
	if not IsMutable(LS) then
	    LS:= ShallowCopy(LS);
	fi;
	Sort(LS, function(x,y) return Size(x[1])<Size(y[1]); 
	    end);
	flag:=0;
	for Lx in LS do
	    x:=Lx[1];
		if Order(x)>= OrdPi then
	    if IsSubgroup( x , Kersnt ) then
			for x in Lx do
				if  IsSubgroup( x , Image( s , x ) ) then
				if  IsSubgroup( x , Image( t , x ) ) then
					sx:=Image( s , x );
					Ordsx:=Order( sx );
					if  Ordsx = Order( Image( t , Intersection( Kers , x ) ) ) * OrdPiOne then
					if  Ordsx = Order( Intersection( sx , tKers ) ) * OrdPiOne then
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
		
    newGens:=GeneratorsOfGroup( H );
    news:=GroupHomomorphismByImagesNC( H , H , newGens , List( newGens , function ( x )
              return Image( s , x );
          end ) );
    newt:=GroupHomomorphismByImagesNC( H , H , newGens , List( newGens , function ( x )
              return Image( t , x );
          end ) );
    D:= Objectify( HapCatOneGroup , rec(  
          sourceMap:=news , 
          targetMap:=newt ) );
    return D;
end);
###############################################
InstallGlobalFunction(QuotientQuasiIsomorph_alt,
    function (C)
    local  s,t,G,H,Kers,Kert,Kersnt,Ims,OrdIms,Imt,OrdPiOne,OrdPiTwo,Ord,
	LN,x,n,i,
	OrderPiOneGx,OrderPiTwoGx,
	epi,newG,newGens,news,newt,D;
	
	
	s:=C!.sourceMap;
    t:=C!.targetMap;
    G:=Source( s );
    Kers:=Kernel( s );
	Ims:=Image( s );
	OrdIms:=Order( Ims );
	Imt:=Image( t );
	Kert:=Kernel( t );
	Kersnt:=Intersection( Kers , Kert );
	OrdPiOne:= Order( HomotopyGroup( C , 1 ) );
	OrdPiTwo:= Order( HomotopyGroup( C , 2 ) );
	Ord:=Order(G)/(OrdPiOne*OrdPiTwo);
	
	#################
	OrderPiOneGx:=function( x )
	    local tsx;
	    tsx:=Image( t , PreImages( s , Intersection( Ims , x ) ) );
		return ( OrdIms*Order( Intersection( tsx , x ) ) )/( Order( Intersection( Ims , x ) )*Order( tsx ) );
	end;	
	#################	
	OrderPiTwoGx:=function( x )
	    local f;
	    f:=NaturalHomomorphismByNormalSubgroup( G , x );
	    return Order( Intersection( Image( f , PreImages( s , Intersection( Ims , x ) ) ) , Image( f , PreImages( t , Intersection( Imt , x ) ) ) ) );
	end;;
	##################
	LN:=NormalSubgroups( G );
	if not IsMutable(LN) then
	    LN:= ShallowCopy(LN);
	fi;
	Sort(LN, function(x,y) return Size(x)>Size(y); 
	    end);
	n:=Length( LN );
	for i in [1..n] do
		x:=LN[i];
		if Order(x) <= Ord then
		if IsSubgroup( x , Image( s , x ) ) then
		if IsSubgroup( x , Image( t , x ) ) then
		if IsSubgroup( x , CommutatorSubgroup( PreImages( s , Intersection( Ims , x ) ) , PreImages( t , Intersection( Imt , x ) ) ) ) then
		if Order( Kersnt ) = Order( Intersection( Kersnt , x ) ) * OrdPiTwo then
		if OrderPiTwoGx( x ) = OrdPiTwo then
		if OrderPiOneGx( x ) = OrdPiOne then
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
	if Order( H )=1 then
		return C;
	fi;
	epi:=NaturalHomomorphismByNormalSubgroup( G ,  H );
    newG :=Image( epi );
    newGens:=GeneratorsOfGroup( newG );
    news:=GroupHomomorphismByImagesNC( newG ,  newG ,  newGens , 
       List( newGens ,  function ( x )
              return Image( epi ,  Image( s ,  PreImagesRepresentative( epi ,  x ) ) );
          end ) );
    newt:=GroupHomomorphismByImagesNC( newG ,  newG ,  newGens , 
       List( newGens ,  function ( x )
              return Image( epi ,  Image( t ,  PreImagesRepresentative( epi ,  x ) ) );
          end ) );
    D:= Objectify( HapCatOneGroup ,  rec( 
          sourceMap:=news , 
          targetMap:=newt ) );
    return D;
end);
###########################################################
InstallGlobalFunction(QuasiIsomorph_alt,
    function (C)
	local D;
	D:=QuotientQuasiIsomorph_alt(C);
	D:=SubQuasiIsomorph_alt(D);
	while Size( D ) < Size( C )  do
        C:=D;
        D:=QuotientQuasiIsomorph_alt( C );
		if Size(D) < Size(C) then
			D:=SubQuasiIsomorph_alt( D );
		fi;
    od;
return D;
end);
##################################################################################
################################################################################
InstallGlobalFunction(XmodToHAP_alt,
   function ( X )
    local  G,Gens,s,t;
    if IsHapCatOneGroup( X )  then
        return X;
    fi;
    if IsComponentObjectRep( X )  then
        if "TailMap" in NamesOfComponents( X ) and "HeadMap" in NamesOfComponents( X )  then
		    G:=X!.Source;
			Gens:=GeneratorsOfGroup(G);
			s:=X!.TailMap;
			t:=X!.HeadMap;
            return Objectify( HapCatOneGroup, rec(
                  sourceMap:=GroupHomomorphismByImagesNC(G,G,Gens,List(Gens,g->Image(s,g))),
                  targetMap:=GroupHomomorphismByImagesNC(G,G,Gens,List(Gens,g->Image(t,g))) 
				  ) );
            
        fi;
    fi;
    Print( "Argument is not a cat-1-group from the Xmod package.\n" );
    return fail;
end);
##################################################################################
##################################################################################
InstallGlobalFunction(IsomorphismCatOneGroups,
	function(C,D)
	local sC,tC,GC,AutGC,g,Gens,
		  sD,tD,GD,
		  m,x,f,flag;
		
	sC:=C!.sourceMap;
	GC:=Source(sC);
	tC:=C!.targetMap;
	sD:=D!.sourceMap;
	tD:=D!.targetMap;
	GD:=Source(sD);
	m:=IsomorphismGroups(GC,GD);
	if m = fail then
		return fail;
	else
		Gens:=GeneratorsOfGroup(GC);
		AutGC:= AutomorphismGroup(GC);
		for x in AutGC do
			f:=x*m;
			flag:=1;
			for g in Gens do
				if ( Image(f,(Image(sC,g))) <> Image(sD,(Image(f,g))) ) or ( Image(f,(Image(tC,g))) <> Image(tD,(Image(f,g))) )then
					flag:=0;
					break;
				fi;
			od;
			if flag =1 then
				return Objectify(HapCatOneGroupHomomorphism,
				   rec(
						source:= C,
						target:= D,
						mapping:=f
					  ));
			fi;
		od;
		return fail;
	fi;
end);
##################################################################################
##################################################################################
InstallGlobalFunction(IsQuasiMorphismCatOneGroups,
function(C,D)
		local 
		PiOneC,PiTwoC,PiOneD,PiTwoD,
		sC,GC,tC,KertC,KersC,ImsC,KersntC,
		sD,GD,tD,tDKerD,OrdPiOneD,OrdPiTwoD,
		Homs,Gens,AllHoms,f,g,flag;
		
		PiOneC:=HomotopyGroup(C,1);
		PiTwoC:=HomotopyGroup(C,2);
		PiOneD:=HomotopyGroup(D,1);
		PiTwoD:=HomotopyGroup(D,2);
		if IdGroup(PiOneC)<>IdGroup(PiOneD) then
			return fail;
		fi;
		if IdGroup(PiTwoC)<>IdGroup(PiTwoD) then
			return fail;
		fi;
		
		sC:=C!.sourceMap;
		GC:=Source(sC);
		tC:=C!.targetMap;
		KertC:=Kernel(tC);
		KersC:=Kernel(sC);
		ImsC:=Image(sC);
		KersntC:=Intersection(KertC,KersC);

		sD:=D!.sourceMap;
		GD:=Source(sD);
		tD:=D!.targetMap; 
		OrdPiOneD:=Order(PiOneD);
		OrdPiTwoD:=Order(PiTwoD);
		tDKerD:=Image(tD,Kernel(sD));
				
		Homs:=[];
		Gens:=SmallGeneratingSet(GC);
		AllHoms:=AllHomomorphisms(GC,GD);
		for f in AllHoms do
			flag:=1;
			for g in Gens do
				if ( Image(f,(Image(sC,g))) <> Image(sD,(Image(f,g))) ) or ( Image(f,(Image(tC,g))) <> Image(tD,(Image(f,g))) ) then
					flag:=0;
					break;
				fi;
			od;
			if flag=1 then
				Add(Homs,f);
			fi;
		od;
		
		for f in Homs do
			if Order(Image(f,KersntC))= OrdPiTwoD then
				if Order(tDKerD)=Order(Intersection(tDKerD,Image(f,ImsC)))*OrdPiOneD then
					return f;
				fi;
			fi;
		od;
	return fail;
end);

	
