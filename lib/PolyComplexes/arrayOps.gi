#(C) Graham Ellis 2009

#####################################################################
#####################################################################
InstallGlobalFunction(ArrayValue,
function(A,x);
# A horrible piece of code!! It inputs an array A and list of integers x.
# It returns the value A[x[1]][x[2]][x[3]]...[x[n]] where n is the length
# of x. 

if Length(x)=1 then return A[x[1]];fi;
if Length(x)=2 then return A[x[2]][x[1]];fi;
if Length(x)=3 then return A[x[3]][x[2]][x[1]];fi;
if Length(x)=4 then return A[x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=5 then return A[x[5]][x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=6 then return A[x[6]][x[5]][x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=7 then return A[x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=8 then return A[x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=9 then return A[x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]];fi;
if Length(x)=10 then return A[x[10]][x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]];fi;

if Length(x)>10 then Print("ArrayValue needs to be implemented for longer lists\n"); return fail;fi;
end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ArrayAssign,
function(A,x,N);
# A horrible piece of code!! It inputs an array A and list of integers x and an object N.
# It sets the value A[x[1]][x[2]][x[3]]...[x[n]] equal to N where n is the length
# of x.

if Length(x)=1 then A[x[1]]:=N;fi;
if Length(x)=2 then A[x[2]][x[1]]:=N;fi;
if Length(x)=3 then A[x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=4 then A[x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=5 then A[x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=6 then A[x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=7 then A[x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=8 then A[x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=9 then A[x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;
if Length(x)=10 then A[x[10]][x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N;fi;

if Length(x)>10 then Print("ArrayValue needs to be implemented for longer lists\n"); return fail;fi;
end);
#####################################################################
#####################################################################

#############################################
#############################################
InstallGlobalFunction(UnboundedArrayAssign,
function(A,x,v)
local

      i,B;

B:=A;
for i in Reversed(x) do
  if not IsBound(B[i]) then B[i]:=[]; fi;
  B:=B[i];
od;

ArrayAssign(A,x,v);
end);
############################################
############################################

#####################################################################
#####################################################################
InstallGlobalFunction(ArrayIterateBreak,
function(Dim);


######
if Dim=1 then return
function(Dims,Func) local i,b;
for i in Dims[1] do b:=Func([i]); 
if b then return [i]; fi; 
od;
end;
return fail;
fi;
######
######
if Dim=2 then return
function(Dims,Func)
local i,j,b;
for i in Dims[1] do
for j in Dims[2] do
b:=Func([i,j]);
if  b then return [i,j]; fi;
od;od;
return fail;
end;
fi;
######
######
if Dim=3 then return
function(Dims,Func)
local i,j,k,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
b:=Func([i,j,k]);
if b then return [i,j,k]; fi; 
od;od;od;
return fail;
end;
fi;
######
if Dim=4 then return
function(Dims,Func)
local i,j,k,l,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
b:=Func([i,j,k,l]);
if b then return [i,j,k,l]; fi; 
od;od;od;od;
return fail;
end;
fi;
######
if Dim=5 then return
function(Dims,Func)
local i,j,k,l,m,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
b:=Func([i,j,k,l,m]);
if b then return [i,j,k,l,m]; fi; 
od;od;od;od;od;
return fail;
end;
fi;
######
if Dim=6 then return
function(Dims,Func)
local i,j,k,l,m,n,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
b:=Func([i,j,k,l,m,n]);
if b then return [i,j,k,l,m,n]; fi; 
od;od;od;od;od;od;
return fail;
end;
fi;
######
if Dim=7 then return
function(Dims,Func)
local i,j,k,l,m,n,p,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
for p in Dims[7] do
b:=Func([i,j,k,l,m,n,p]);
if b then return [i,j,k,l,m,n,p]; fi;
od;od;od;od;od;od;od;
return fail;
end;
fi;
######
if Dim=8 then return
function(Dims,Func)
local i,j,k,l,m,n,p,q,b;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
for p in Dims[7] do
for q in Dims[8] do
b:=Func([i,j,k,l,m,n,p,q]);
if b then return [i,j,k,l,m,n,p,q]; fi;
od;od;od;od;od;od;od;od;
return fail;
end;
fi;
######




end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ArrayIterate,
function(Dim);

######
if Dim=1 then return
function(Dims,Func) local i;
for i in Dims[1] do Func([i]); od;
end;
fi;
######
if Dim=2 then return
function(Dims,Func)
local i,j;
for i in Dims[1] do
for j in Dims[2] do 
Func([i,j]); 
od;od;
end;
fi;
######
if Dim=3 then return
function(Dims,Func)
local i,j,k;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
Func([i,j,k]);
od;od;od;
end;
fi;
######
if Dim=4 then return
function(Dims,Func)
local i,j,k,l;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
Func([i,j,k,l]);
od;od;od;od;
end;
fi;
######
if Dim=5 then return
function(Dims,Func)
local i,j,k,l,m;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
Func([i,j,k,l,m]);
od;od;od;od;od;
end;
fi;
######
if Dim=6 then return
function(Dims,Func)
local i,j,k,l,m,n;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
Func([i,j,k,l,m,n]);
od;od;od;od;od;od;
end;
fi;
######
if Dim=7 then return
function(Dims,Func)
local i,j,k,l,m,n,p;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
for p in Dims[7] do
Func([i,j,k,l,m,n,p]);
od;od;od;od;od;od;od;
end;
fi;
######
if Dim=8 then return
function(Dims,Func)
local i,j,k,l,m,n,p,q;
for i in Dims[1] do
for j in Dims[2] do
for k in Dims[3] do
for l in Dims[4] do
for m in Dims[5] do
for n in Dims[6] do
for p in Dims[7] do
for q in Dims[8] do
Func([i,j,k,l,m,n,p,q]);
od;od;od;od;od;od;od;od;
end;
fi;



end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ArrayAssignFunctions,
function(x);

if x=1 then return function(A,x,N); A[x[1]]:=N; end; fi;
if x=2 then return function(A,x,N); A[x[2]][x[1]]:=N; end; fi;
if x=3 then return function(A,x,N); A[x[3]][x[2]][x[1]]:=N; end; fi;
if x=4 then return function(A,x,N); A[x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=5 then return function(A,x,N); A[x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=6 then return function(A,x,N); A[x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=7 then return function(A,x,N); A[x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=8 then return function(A,x,N); A[x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=9 then return function(A,x,N); A[x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;
if x=10 then return function(A,x,N); A[x[10]][x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]:=N; end; fi;

if Length(x)>10 then Print("ArrayValueFunctions needs to be implemented for longer lists\n"); return fail;fi;

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(ArrayValueFunctions,
function(x);

if x=1 then return function(A,x); return A[x[1]]; end; fi;
if x=2 then return function(A,x); return A[x[2]][x[1]]; end; fi;
if x=3 then return function(A,x); return A[x[3]][x[2]][x[1]]; end; fi;
if x=4 then return function(A,x); return A[x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=5 then return function(A,x); return A[x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=6 then return function(A,x); return A[x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=7 then return function(A,x); return A[x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=8 then return function(A,x); return A[x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=9 then return function(A,x); return A[x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;
if x=10 then return function(A,x); return A[x[10]][x[9]][x[8]][x[7]][x[6]][x[5]][x[4]][x[3]][x[2]][x[1]]; end; fi;

if Length(x)>10 then Print("ArrayValueFunctions needs to be implemented for longer lists\n"); return fail;fi;

end);
#####################################################################
#####################################################################

##############################################################
##############################################################
InstallGlobalFunction(ContractArray,
function(AA);

if ArrayDimension(AA)=2 then return
HomotopyEquivalentSmallerSubMatrix(AA,AA*0);
fi;

if ArrayDimension(AA)=3 then return
HomotopyEquivalentSmallerSubArray3D(AA,AA*0);
fi;

return HomotopyEquivalentSmallerSubArray(AA,AA*0);

end);
##############################################################
##############################################################

##############################################################
##############################################################
InstallGlobalFunction(ContractPermArray,
function(AA);

if ArrayDimension(AA)=2 then return
HomotopyEquivalentSmallerSubPermMatrix(AA,AA*0);
fi;

if ArrayDimension(AA)=3 then return
HomotopyEquivalentSmallerSubPermArray3D(AA,AA*0);
fi;

#return HomotopyEquivalentSmallerSubPermArray(AA,AA*0);

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubArray,
function(AA,SS)
local
	Dim, dim, dim1,
	Dims,dims, dimsSet, revdimsSet,
        IsRemovableCube,
	ArrayValueDim,ArrayIt,ArrayAssignDim,ArrayValueDim1,
	cart, sizecart, 
	Ball,
	bool,
	move,Fun,Fun2,
	correction,
	Elts,
	x,w;

if ArrayDimension(AA)=2 then  return 
HomotopyEquivalentSmallerSubMatrix(AA,SS); fi;
if ArrayDimension(AA)=3 then  return
HomotopyEquivalentSmallerSubArray3D(AA,SS); fi;

AA:=FrameArray(AA);
SS:=FrameArray(SS);
dim:=ArrayDimension(AA);
dim1:=dim-1;
dims:=ArrayDimensions(AA);
dimsSet:=List(dims,d->[2..d-1]);
revdimsSet:=List(dimsSet,d->Reversed(d));
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
cart:=Cartesian(List([1..dim],a->[-1,0,1]));
RemoveSet(cart,List([1..dim],i->0));
sizecart:=Size(cart);

Ball:=[-1,0,1]; #!!!!!!!!!!!
for x in [2..dim] do
Ball:=List(Ball,x->StructuralCopy(Ball));
od;
Ball:=ArrayToPureCubicalComplex(Ball,1);
correction:=List([1..dim],i->2);
##############################
IsRemovableCube:=function(A,x);
if ArrayValueDim(SS,x)=1 then return false; fi;
return IsContractibleCube_higherdims(A,A,dims,x,dim,dim1,ArrayValueDim,ArrayAssignDim,Ball,correction,cart);
end;
##############################

######################
Fun:=function(x);
if IsRemovableCube(AA,x) then
ArrayAssignDim(AA,x,0);
bool:=true;
fi;
end;
######################

bool:=true;
##########################
while bool and ArraySum(AA)>10000 do  
				#10^5 could be replaced by a better value?
bool:=false;

ArrayIt(dimsSet,Fun);
if bool then
ArrayIt(revdimsSet,Fun);
fi;

od;
##########################

Elts:=[];
Fun2:=function(x);
if ArrayValueDim(AA,x)=1 then Add(Elts,x); fi;
end;

##########################
while bool do
bool:=false;

Elts:=[];
ArrayIt(dimsSet,Fun2);
for x in Elts do
Fun(x);
od;
if bool then
for x in Reversed(Elts) do
Fun(x);
od;
fi;

od;
##########################

####################################

SS:=UnframeArray(SS);
AA:=UnframeArray(AA);
return AA;
end);
##############################################################
##############################################################

##############################################################
##############################################################
InstallGlobalFunction(HomotopyEquivalentSmallerSubPermArray,
function(AA,SS);

if ArrayDimension(AA)=2 then  return
HomotopyEquivalentSmallerSubPermMatrix(AA,SS); fi;
if ArrayDimension(AA)=3 then  return
HomotopyEquivalentSmallerSubPermArray3D(AA,SS); fi;
end);
##############################################################
##############################################################



##############################################################
##############################################################
InstallGlobalFunction(ContractibleSubArray,
function(AAA);
return HomotopyEquivalentLargerSubArray(AAA,AAA*0);
end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(HomotopyEquivalentLargerSubArray,
function(AAA,SSS)
local
        AA,dim, dim1,
        dims,dimsSet,revdimsSet,
        IsAddableCube,
	ArrayValueDim,ArrayValueDim1,
	ArrayAssignDim,ArrayIt,
        cart, sizecart, 
        bool,
        move,
	Fun1,Fun2,
	S,
	Ball,correction,tst,
	x,w,start;

if ArrayDimension(AAA)=2 then 
return HomotopyEquivalentLargerSubMatrix(AAA,SSS);fi;

if ArrayDimension(AAA)=3 then
return HomotopyEquivalentLargerSubArray3D(AAA,SSS);fi;

AA:=FrameArray(AAA);
S:=FrameArray(SSS);
dim:=ArrayDimension(AA);
dim1:=dim-1;
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayValueDim1:=ArrayValueFunctions(dim1);
ArrayAssignDim:=ArrayAssignFunctions(dim);
ArrayIt:=ArrayIterate(dim);
dims:=ArrayDimensions(AA);
dimsSet:=List(dims,d->[1..d]);
revdimsSet:=List(dimsSet,d->Reversed(d));
cart:=Cartesian(List([1..dim],a->[-1,0,1]));
sizecart:=Size(cart);

Ball:=[-1,0,1]; #!!!!!!!!!!!!!!
for x in [2..dim] do
Ball:=List(Ball,x->StructuralCopy(Ball))*0;
od;
Ball:=PureCubicalComplex(Ball);
correction:=List([1..dim],i->2);
##############################
IsAddableCube:=function(A,S,x);
if ArrayValueDim(S,x)=1 then return false; fi;
return IsContractibleCube_higherdims(A,S,dims,x,dim,dim1,ArrayValueDim,ArrayAssignDim,Ball,correction,cart);
end;
##############################

#################
Fun1:=function(x);
if ArrayValueDim(AA,x)=1 then
start:=x;  fi;
end;
#################

#################
Fun2:=function(x);
if IsAddableCube(AA,S,x) then
#w:=ArrayValueDim1(S,x{[2..dim]});
#w[x[1]]:=1;
ArrayAssignDim(S,x,1);
bool:=true;
fi;
end;
#################

bool:=true;

##########If S is empty then ...#########
if ArraySum(S)=0 then 
#S:=AA*0;

ArrayIt(revdimsSet,Fun1);

if IsBound(start) then
ArrayAssignDim(S,start,1);
else
bool:=false;
fi;

fi;
##########Now is is probably not empty##

while bool do
bool:=false;
ArrayIt(dimsSet,Fun2);
if bool then
ArrayIt(revdimsSet,Fun2);
fi;
od;

return UnframeArray(S);
end);
##############################################################
##############################################################

##############################################################
##############################################################
InstallGlobalFunction(HomotopyEquivalentLargerSubPermArray,
function(AAA,SSS);

if ArrayDimension(AAA)=2 then
return HomotopyEquivalentLargerSubPermMatrix(AAA,SSS);fi;

if ArrayDimension(AAA)=3 then
return HomotopyEquivalentLargerSubPermArray3D(AAA,SSS);fi;
end);
##############################################################
##############################################################



##############################################################
##############################################################
InstallGlobalFunction(FrameArray,
function(A);

#if IsInt(A[1]) then
if not IsList(A[1]) then 
return Concatenation([0],A,[0]);
else
return
Concatenation([FrameArray(A[1]*0)],
List(A,a->FrameArray(a)),
[FrameArray(A[1]*0)]);
fi;

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(UnframeArray,
function(A);

#if IsInt(A[1]) then
if not IsList(A[1]) then 
return A{[2..Length(A)-1]};
else
return
List(A{[2..Length(A)-1]},a->UnframeArray(a));
fi;

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(ArraySum,
function(A) 
local 
	sz;

sz:=Sum(A);
if IsInt(sz) then return sz;
else return ArraySum(sz);
fi;

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(ArrayDimension,
function(A);

#if IsInt(A) then return 0;
if not IsList(A) then return 0;
else return 1 +ArrayDimension(A[1]);
fi;

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(ArrayDimensions,
function(A) ;

#if IsInt(A) then return []; fi;
if not IsList(A) then return []; fi;
return Concatenation(ArrayDimensions(A[1]),[Length(A)]);

end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(IsContractibleCube_higherdims,
function(B,A,dims,x,dim,dim1,ArrayValueDim,ArrayAssignDim,Ball,correction,cart)  #A is a sub array of B
local yy,y,z,cnt,flt,CC,w;

if ArrayValueDim(B,x)=0 then return false; fi;
Ball:=PureCubicalComplex(Ball!.binaryArray*0);

cnt:=0;
for y in cart do
z:=x+y;yy:=y+correction;
flt:=Flat(y);
if ArrayValueDim(A,z) = 1 then
ArrayAssignDim(Ball!.binaryArray,yy,1);
  if  Length(Filtered(flt,a->a=0))=dim-1 then cnt:=cnt+1;fi;
fi;
od;

################
################Dimensions 2 and 3
if dim<=3 then
if not PathComponentOfPureCubicalComplex(Ball,0)=1 then return false; fi;
if not EulerCharacteristic(Ball)=1 then return false; fi;
return true;
fi;
################
################

if IsZero(Flat(Ball!.binaryArray)) then return false;fi;
if cnt=2*dim then return false; fi;
if cnt=2*dim-1 then return true; fi;


if PathComponentOfPureCubicalComplex(Ball,0)>1 then
Unbind(Ball!.pathCompBinList);return false;fi;

if Sum(Flat(Ball!.binaryArray))<4 then return true; fi;
if not EulerCharacteristic(Ball)=1 then  return false; fi;

CC:=TensorWithIntegersModP(ChainComplex(Ball),2);

for z in [1..dim1] do
if not Homology(CC,z)=0 then return false; fi;
#AM I SURE THAT WORKING MOD 2 IS OK!!!
od;

return true;
end);
##############################################################
##############################################################


##############################################################
##############################################################
InstallGlobalFunction(Array,
function(A,f);

if IsList(A) and IsInt(A[1]) then return List(A,f);fi;
#if IsList(A) and not IsList(A[1]) then return List(A,f);fi;
return List([1..Length(A)],i->Array(A[i],f));
end);

##############################################################
##############################################################

##################################################
##################################################
InstallGlobalFunction(PermuteArray,
function(A,pi)
local B,x,dim,dims,dimsSet,Fun,
      ArrayValueDim,ArrayIt,ArrayAssignDim,
      d, NewDimsSet;


dim:=ArrayDimension(A);
dims:=ArrayDimensions(A);
dimsSet:=List(dims,d->[1..d]);
ArrayValueDim:=ArrayValueFunctions(dim);
ArrayIt:=ArrayIterate(dim);
ArrayAssignDim:=ArrayAssignFunctions(dim);

NewDimsSet:=List([1..dim],n->dimsSet[n^pi]);
B:=0;
for d in [1..dim] do
B:=List(NewDimsSet[d],i->StructuralCopy(B));
od;


################
Fun:=function(x) local y;
y:=List([1..dim],a->x[a^pi]);
ArrayAssignDim(B,x,ArrayValueDim(A,y));
end;
################
ArrayIt(NewDimsSet,Fun);

return B;
end);
##################################################
##################################################


