#code due to Sebastian Schoennenbeck (with minor modifications)

####################################################################
####################################################################
ResolutionFromWellRoundedComplex:=function(length)
	local stabilizer,IsEq, k,j, evenstabilizer, action, dimension, i, PseudoBoundary, boundary,g,group, actionREC;
#stabilizers,elts,DIMS,evenstabilizers,BoundaryComponent,Gens;
	#Read(file);
	stabilizer:=function(k,j) 
		if Length(stabilizers[k+1][AbsInt(j)]) = 0 then
		  return TrivialSubgroup(Group(elts[1]));
		else
		  return Group(stabilizers[k+1][AbsInt(j)]);
		fi;
	end;

	evenstabilizer:=function(k,j)
		if Length(evenstabilizers[k+1][AbsInt(j)]) = 0 then
		  return TrivialSubgroup(Group(elts[1]));
		else
		  return Group(evenstabilizers[k+1][AbsInt(j)]);
		fi;
	end;

	
actionREC:=List([1..Length(DIMS)],jj->List([1..DIMS[jj]],x->[])  );;
         
   
	action:=function(k,j,g)
		local absj,id,r,u,H;
		absj:=AbsInt(j);
##########################
if IsBound(actionREC[k+1][absj][g])  then
return actionREC[k+1][absj][g]; fi;
##########################
		H:=stabilizer(k,absj);
		id:=CanonicalRightCosetElement(H,Identity(H));
		r:=CanonicalRightCosetElement(H,elts[g]^-1);
		r:=id^-1*r;
		u:=r*elts[g];
		if u in evenstabilizer(k,absj) then 
                actionREC[k+1][absj][g]:=1;
			return 1;
		else
                actionREC[k+1][absj][g]:=-1;
			return -1;
		fi;
	end;
	
	dimension:=function(i)
		if 1<= i+1 and i+1 <=Length(DIMS) then
			return DIMS[i+1];
		fi;
		return 0;
	end;
	
	PseudoBoundary:=List([1..Length(DIMS)],i->[1..dimension(i)]);
	
	#Read("myprocs.gi");
	
	#This procedure will compare elements in permutationmodules with nontrivial stabilizers
	IsEq:=function(k,v,w) 
	#k being the dimension, v and w two elements of the form [i,g] in the usual HAP-structure 
	local i;
	if not AbsInt(v[1])=AbsInt(w[1]) 
		then return false; fi;  
			#makes sure both elements are multiples of the same "free" generator
	i:=AbsInt(v[1]);
	if not elts[v[2]]^(-1)*elts[w[2]] in stabilizer(k,i) 
		then return false; fi;   
			#if ge_i neq he_i returns false
	if v[1]=w[1] and elts[v[2]]^(-1)*elts[w[2]] in evenstabilizer(k,i) 
		then return true; fi;       
			#if both elements act the same way on the orientation and v and w are of the same orientation returns true
	if v[1]=-w[1] and not elts[v[2]]^(-1)*elts[w[2]] in evenstabilizer(k,i)
		then return true; fi; 
			#if the group elements operate not in the same way on the orientation and v[1]=-w[1] this returns true
	return false;
	end;

	#elts should if possible be duplicate free (I dont know if this is necessary but better safe than sorry)

	#The next procedure will compute the boundary homomorphism. We will previously need a list BoundaryComponent which lists all k-1 dimensional cells in the boundary of a k dimensional cell (duplicate free) as well as the list PseudoBoundary:=List([1..lngth],i->[1..Dimension(i)]); in which we will later save the images under the boundary

	boundary:=function(k,mm) 
		#Print(PseudoBoundary);
		#k the dimension, mm corresponds to the abs(mm)-th "free" generator
	local b,bndbnd,x,y,z,n,bnd,signedbnd,tmp,m,bool,pos;
	if k=0 then return []; fi;
	m:=AbsInt(mm); 
	if not IsInt(PseudoBoundary[k][m]) 
		then if mm>0
			then return PseudoBoundary[k][m];
			else return NegateWord(PseudoBoundary[k][m]); 
		fi;
	fi;
	#This is the case if k or m is out of bounds or we have already computed the correct boundary and saved it in PB. This way we will later use a lookup-table rather than run through the entire recursion each time. It also asures consistency in case we have duplicates in elts
	bnd:=StructuralCopy(BoundaryComponent[k][m]); 
	#Now we will insert the signs via recursion and the condition that boundary^2=0
	if k=1
		then bnd[1][1]:=-bnd[1][1]; 
	fi;

	#Initialising the recursion by choosing some orientation on the 1-dim cells

	if k>1 
		then bndbnd:=[]; 
			#in bndbnd we will save all we know about boundary^2 so far
		bnd:=SSortedList(bnd); 
		signedbnd:=[bnd[1]];
		RemoveSet(bnd,bnd[1]);
		x:=ShallowCopy(signedbnd[1]);
		#The following works even if not all products are already in elts.
		b:=StructuralCopy(boundary(k-1,x[1]));
		
		for j in [1..Length(b)] do
		  y:=b[j];
		  g:=elts[x[2]]*elts[y[2]];
		  pos:=Position(elts,g);
		  if pos = fail then
		    Append(elts,[g]);
		    pos:=Length(elts);
		  fi;
		  b[j][2]:=pos;
		od;
		    
		#b:=List(b,y->[y[1],Position(elts,elts[x[2]]*elts[y[2]])]); 
			#This finds the correct elements in the boundary of x
		Append(bndbnd,b);
		while Length(bnd)>0 
			do x:=Random(bnd); bool:=true;    
				#We are taking random elements because we can only work with elements that have non-trivial intersection with those we already looked at.
			b:=StructuralCopy(boundary(k-1,x[1]));
			
			for j in [1..Length(b)] do
			  y:=b[j];
			  g:=elts[x[2]]*elts[y[2]];
			  pos:=Position(elts,g);
			  if pos = fail then
			    Append(elts,[g]);
			    pos:=Length(elts);
			  fi;
			  b[j][2]:=pos;
			od;
			#b:=List(b,y->[y[1],Position(elts,elts[x[2]]*elts[y[2]])]);
			for y in bndbnd do 
				for z in b do	
					if bool and IsEq(k-2,y,z) 
						then Append(signedbnd,[[-x[1],x[2]]]);
						Append(bndbnd,NegateWord(b)); 
						RemoveSet(bnd,x);
						bool:=false; 
							#This is the case where an element in the boundary of x already exists in bndbd and we therefore need to take x negatively to assure that boundary^2 will be zero.
						else if bool and IsEq(k-2,[-y[1],y[2]],z)
                                        	then Append(signedbnd,[[x[1],x[2]]]);
                                        	Append(bndbnd,b);
                                        	RemoveSet(bnd,x);
						bool:=false;
							#The dual case. The negative of an element in the boundary of x already appears in bndbnd and we have to take x positively to balance it out.
                                		fi;
					fi;
				od;
			od;
		od;
		bnd:=signedbnd;
	fi;
	PseudoBoundary[k][m]:=bnd;    
		#Write to our lookup table 
	return boundary(k,mm); 
	end;	


	return Objectify(HapNonFreeResolution,rec(dimension:=dimension,elts:=elts, group:=Group(Gens),stabilizer:=stabilizer,homotopy:=fail, action:=action, hasse:=fail,boundary:=boundary, properties:=[["type","resolution"],["reduced",true],["length",length],["characteristic",0]]));

end;
##########################################################################
##########################################################################

HAP_KOMPLEKS123:=ResolutionFromWellRoundedComplex(1000);
