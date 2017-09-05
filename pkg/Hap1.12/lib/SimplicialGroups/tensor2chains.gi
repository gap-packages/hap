InstallGlobalFunction(TensorProductOfChainComplexes,
function(arg)

	local TensorProductOfTwoChainComplexes,
		  i,LR,R;

TensorProductOfTwoChainComplexes := function(R,S)
local 
	rdim,rbound,sdim,sbound,N,
	Dimen,Bound,count,n,k,M,i,j,ii,jj,x,Rew,beg1,beg2,LH,LV,
	dim,SearchPosition,dhr,dvs,
	Dimension,Boundary;
	
        ########################## Modified by Graham
	rdim := function(n);
        if n<0 or n> Length(R) then return 0; fi;
        return R!.dimension(n);
        end;
        ##########################
	rbound := R!.boundary;
        ########################## Modified by Graham
        sdim := function(n);
        if n<0 or n> Length(S) then return 0; fi;
        return S!.dimension(n);
        end;
        ##########################

	sbound := S!.boundary;
	#N := Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length")); #CHANGED BY GRAHAM
N := 1+EvaluateProperty(R,"length")+EvaluateProperty(S,"length");

##############################################
dim := function(i,j)
	if i<0 or j <0 then
		return 0;
	fi;
	return rdim(i)*sdim(j);
end;
###############################################
	Dimen := [];
	for n in [0..N] do
		count := 0;
		for j in [0..n] do;
			count := count+dim(n-j,j);
		od;
		Dimen[n+1] := count;
	od;
######################################################
Dimension := function(n)
	return Dimen[n+1];
end;
###################################################
###################################################
SearchPosition := function(n,k) 
	local count,ord,j,b,x,ii,jj; 
	
	if k>Dimension(n) then
		return fail;
	fi;
	count := 0;
	for j in [0..n] do
		ord := k-count;
		count := count+dim(n-j,j);
		if k <=count then
			b := sdim(j);
			x := Int((ord-1)/b);
			ii := x+1;
			jj := ord-x*b;
			return [n-j,j,ii,jj];
			break;
		fi;
	od;
end;	
########################################################
########################################################
dhr := function(i,j,ii,jj)
local a,b,L,x,LR;

    if i=0 then
	   return [];
    fi;
    a := rdim(i-1);
    b := sdim(j);
    L :=  List([1..a*b],x->0);
	LR := rbound(i,ii);
	for x in [1..a] do
       L[(x-1)*b+jj] := LR[x];
	od;
	return L;
end;
########################################################
dvs := function(i,j,ii,jj)
local a,b,L,x,LS,sign;
    if j =0 then
	   return [];
    fi;
	if (i mod 2)=0 then
		sign := 1;
	else
		sign := -1;
	fi;
    a := rdim(i);
    b := sdim(j-1);
    L := List([1..a*b],x->0);
    LS := sbound(j,jj);
    for x in [1..b] do
       L[(ii-1)*b+x] := sign*LS[x];
	od;
	return L;
end;	
##########################################################
	Bound := List([0..N],x->[]); 
	Bound[1][1] := [0];
	for n in [1..N] do 
		for k in [1..Dimension(n)] do
			M :=	SearchPosition(n,k);
			i := M[1];
			j := M[2];
			ii := M[3];
			jj := M[4];
			Rew := List([1..Dimension(n-1)],x->0);
			beg1 := 0;
			for x in [0..(j-2)] do
				beg1 := beg1+ dim((n-1)-x,x);
			od;
			beg2 := beg1+	dim(i,j-1);			
			LV := dvs(i,j,ii,jj);
			LH := dhr(i,j,ii,jj);
			for x in [1..dim(i,j-1)]do
				Rew[beg1+x] := LV[x];
			od;
			for x in [1..dim(i-1,j)] do
				Rew[beg2+x] := LH[x];
			od;
		    Bound[n+1][k] := Rew;
		od;
	od;
################################################################
Boundary := function(n,k)
	return Bound[n+1][k];
end;
################################################################

return Objectify( HapChainComplex, rec(
				boundary := Boundary,
				dimension := Dimension,
				properties :=  [ [ "length",N-1],
                    [ "type", "chainComplex" ], 
                    [ "characteristic",0 ] ] ) );
end;
################################################################
	
	if Length(arg)=1 then
		LR := arg[1];
		if IsList(LR) then
			R := LR[1];
			for i in [2..Length(LR)] do
				R := TensorProductOfTwoChainComplexes(R,LR[i]);
			od;
			return R;
		else
			return LR;
		fi;
    else
		R := arg[1];
		for i in [2..Length(arg)] do
			R := TensorProductOfTwoChainComplexes(R,arg[i]);
		od;
		return R;		
	fi;
end);					
		
	
			
			
		
	
   
   
   
    




