#(C) Graham Ellis, 2005-2006


#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionPrimePowerGroup,
function(arg)
local
	G,n,
	eltsG,
	gensG,
	Dimension,
	Boundary,
	PartialHomotopy,
	Homotopy,
	PseudoBoundary,
	PseudoBoundaryAsVec,
	WordToVectorList,
	VectorListToWord,
	prime, pp,
	BoundaryMatrix,
	MT,
	GactMat,
	ZGbasisOfKernel,
	one,
	InverseFlat,
	ComplementaryBasis,
	zero,
	pcgens,
	BoundaryMatrices,
	Echelonize,
	EchelonMatrices,
	SolutionMatBoundaryMatrices,
	Toggle,
	SMBM,
	SpaceSave,
	g,h,i,x,tmp;


G:=arg[1];
n:=arg[2];
if Length(arg)>2 then SpaceSave:=true; else SpaceSave:=false; fi;
tmp:=SSortedList(Factors(Order(G)));
if Length(tmp)>1 then 
Print("This function can only be applied to small prime-power groups. \n");
return fail;
fi;
prime:=tmp[1];
pp:=Order(G);
one:=Identity(GaloisField(prime));
zero:=0*one;

gensG:=ReduceGenerators(GeneratorsOfGroup(G),G);
eltsG:=Elements(G);
MT:=MultiplicationTable(eltsG);

pcgens:=Pcgs(G);
pcgens:=ReduceGenerators(pcgens,G);
pcgens:=List(pcgens,x->Position(eltsG,x));

PseudoBoundary:=[];
PseudoBoundaryAsVec:=[];
for i in [1..n] do
PseudoBoundary[i]:=[];
PseudoBoundaryAsVec[i]:=[];
od;

#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return Length(PseudoBoundaryAsVec[i]);
end;
#####################################################################

#####################################################################
Boundary:=function(i,j);
if i<0 then return []; fi;
if j>0 then
return PseudoBoundary[i][j]; 
else return
 NegateWord(PseudoBoundary[i][-j]);
 fi;
end;
#####################################################################

#####################################################################
WordToVectorList:=function(w,k)	#w is a word in R_k. 
local v,x,a;			#v is a list of vectors mod p.
#v:=List([1..Dimension(k)],x->List([1..pp],y->0) );
v:=ListWithIdenticalEntries(Dimension(k),ListWithIdenticalEntries(pp,0) );


for x in w do
a:=AbsoluteValue(x[1]);
v[a][x[2]]:=v[a][x[2]] + SignInt(x[1]);
od;

return v mod prime;
end;
#####################################################################

for x in gensG do
Add(PseudoBoundary[1], [[-1,1],[1,Position(eltsG,x)]]  );
Add(PseudoBoundaryAsVec[1], Flat(WordToVectorList
( [[-1,1],[1,Position(eltsG,x)]] ,0))*one);
od;

#####################################################################
VectorListToWord:=function(v)
local w, i, x;

w:=[];
for x in [1..Length(v)] do
for i in [1..Length(v[x])] do
if not v[x][i]=0 then 
Append(w, MultiplyWord(v[x][i],[ [x,i]   ]));
fi;
od;
od;

return w;
end;
#####################################################################


#####################################################################
GactMat:=function(g,tB)
local k,q,h,C;

C:=[];

k:=0;
for q in [0..(-1+Length(tB)/pp)] do

for h in [1..pp] do
C[k+MT[g][h]]:=tB[k+h];
od;
k:=k+pp;
od;

ConvertToMatrixRepNC(C);

return C;
end;
#####################################################################



#####################################################################
BoundaryMatrix:=function(k)	#Returns the matrix of d_k:R_k->R_k-1
local B,M,r, b, i, g,j,gB;		
				#M is actually the transpose of the matrix!
B:=TransposedMat(PseudoBoundaryAsVec[k]);

M:=[];

for g in [1..pp] do
gB:=TransposedMat(GactMat(g,B));

for i in [0..Dimension(k)-1] do
M[i*pp+g]:=gB[i+1];
od;
od;


return M;
end;
#####################################################################


#####################################################################
InverseFlat:=function(v)
local w,x,cnt,i;

w:=[];
cnt:=0;
while cnt< Length(v) do
x:=[];
for i in [cnt+1..cnt+pp] do
Add(x,v[i]);
cnt:=cnt+1;
od;
Add(w,x);
od;
return w;
end;
#####################################################################

#####################################################################
ComplementaryBasis:=function(arg)
local B, NS, BC, heads,ln, i, v,zeroheads;


B:=arg[1];
if Length(arg)>1 then
NS:=arg[2];
fi;

ConvertToMatrixRepNC(B,prime); 
heads:=SemiEchelonMat(B).heads;
ln:=Length(B[1]);
BC:=[];

zeroheads:=Filtered([1..ln],i->heads[i]=0);

if Length(arg)=1 then
	for i in zeroheads do
	v:=ListWithIdenticalEntries(ln,zero);
	v[i]:=one;
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
else
	for i in zeroheads do
	v:=NS.vectors[NS.heads[i]];
	ConvertToVectorRep(v,prime);
	Add(BC,v);
	od;
fi;



return BC;
end;
#####################################################################

#####################################################################
ZGbasisOfKernel:=function(k)		#The workhorse!
local  	tB,i, v, g, h, b,bb, ln, B, B1, B2,NS, 
	Bcomp, BndMat;

BndMat:=BoundaryMatrix(k); 


ConvertToMatrixRepNC(BndMat,prime); 
NS:=SemiEchelonMat(NullspaceMat(BndMat));
B:=NS.vectors;

tB:=TransposedMat(B);
Bcomp:=ComplementaryBasis(B);


					
for g in pcgens do     	 
b:=TransposedMat(tB-GactMat(g,tB));
for i in [1..Length(b)] do
Add(Bcomp,b[i]);
od;
if SpaceSave then
#Include the next two lines to save space.
ConvertToMatrixRep(Bcomp);
Bcomp:=MutableCopyMat(SemiEchelonMat(Bcomp).vectors);
fi;
od;							


B1:=ComplementaryBasis(Bcomp,NS);

B1:=List(B1,v->List(v,x->IntFFE(x)));
return B1;
end;
#####################################################################

for i in [2..n] do


for x in ZGbasisOfKernel(i-1) do
Add(PseudoBoundary[i], VectorListToWord(InverseFlat(x))   );
Add(PseudoBoundaryAsVec[i], x*one   );
od;
od;


#####################################################################

#####################################################################

#####################################################################


#####################################################################
Echelonize:=function()
local i,  Mt, T;

BoundaryMatrices:=[];
for i in [1..n] do
BoundaryMatrices[i]:=TransposedMat(BoundaryMatrix(i)*one);
od;

EchelonMatrices:=[];

for i in [1..n] do
#We want to solve XM=W, so work on (Mt)(Xt)=(Wt)
# where
ConvertToMatrixRepNC(BoundaryMatrices[i],prime);
T:=SemiEchelonMatTransformation(BoundaryMatrices[i]);
EchelonMatrices[i]:=[T.coeffs,Reversed(T.heads),T.vectors];
od;

end;
#####################################################################


#####################################################################
SolutionMatBoundaryMatrices:=function(m,w)
local h,u,v,i,cnt,pos,col,row,diff;

ConvertToVectorRep(w);
v:=EchelonMatrices[m][1]*w;
h:=StructuralCopy(EchelonMatrices[m][2]);
u:=ListWithIdenticalEntries(Length(h),0*one);

while not Sum(h)=0 do

col:=PositionProperty(h,x->not x=0);
row:=h[col];
h[col]:=0;
diff:=EchelonMatrices[m][3][row]*u;
pos:=Length(u)+1-col;
u[pos]:=v[row]-diff;
ConvertToVectorRep(u);
od;

#return
#SolutionMat((BoundaryMatrices[m]),w*one);

return u;
end;
#####################################################################


Toggle:=true;
#####################################################################
SMBM:=function(m,w);

if Toggle then Echelonize(); Toggle:=false; fi;

return SolutionMatBoundaryMatrices(m,w);
end;
#####################################################################

#####################################################################
Homotopy:=function(k,w)         #assume w is in kernel d_n
local u,v,s;

if Toggle then Echelonize(); Toggle:=false; fi;

v:=Flat(WordToVectorList([w],k))*one;
ConvertToVectorRep(v,prime);
if k=0 then 

u:=StructuralCopy(v);
Apply(u,i->IntFFE(i));
s:=Sum(u);
u:=ListWithIdenticalEntries(Length(v),zero);
u[1]:=one*s;
ConvertToVectorRep(u);
v:=v-u;
v:=SolutionMatBoundaryMatrices(k+1,v);

else

v:=v-SolutionMatBoundaryMatrices(k,BoundaryMatrices[k]*v);
v:=SolutionMatBoundaryMatrices(k+1,v);

fi;
Apply(v,i->IntFFE(i));

if not v=fail then v:=VectorListToWord(InverseFlat(v)); fi;
return v;
end;
#####################################################################



return Objectify(HapResolution,
	        rec(
		dimension:=Dimension,
		boundary:=Boundary,
		homotopy:=Homotopy,
		elts:=eltsG,
		group:=G,
		properties:=
			[["length",n],
			 ["reduced",true],
			 ["type","resolution"],
			 ["characteristic",prime],
			 ["isMinimal",true]],
		solutionMatBoundaryMatrices:=
		  SMBM));
end);
#####################################################################
#####################################################################
