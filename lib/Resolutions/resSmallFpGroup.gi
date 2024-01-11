#(C) Irina Kholodna, 2001 (converted to HAP format by Graham Ellis)

#####################################################################
InstallGlobalFunction(ResolutionSmallGroup,
function(arg)
local
	G,n,
	table,
	Table,
	NullList,
	LastFactor,
	WordsToIntegers,
	Der,
	IntegersToWords,
	VectorPermutations,
	GenerateAsZGModule,
	FoxMatrix,
	Is_ZLinearComb,
	MinimalModuleGenerators,
	IdentityModuleZBasis,
	IdentityModuleGenerators,
	Resolution,
	IrinasRes,
	ExtendIrinasRes,
	Dimension,
	Boundary,
	EltsG, OrderG, 
	Int2Pair,
	FoxM,
	Charact,
	FpToPerm,
	pG; 

G:=arg[1];
n:=arg[2];
if Length(arg)>2 then Charact:=arg[3]; else Charact:=0; fi;
OrderG:=Order(G);

if IsFpGroup(G) then
FpToPerm:=IsomorphismPermGroup(G);
else
G:=Group(MinimalGeneratingSet(G));
FpToPerm:=IsomorphismPermGroup(G);
fi;
pG:=Image(FpToPerm,G);

EltsG:=List(Elements(pG),x->PreImagesRepresentative(FpToPerm,x));

########################START OF IRINA'S CODE #######################
#####################################################################
#  GAP V.4 FUNCTIONS FOR COMPUTING THE MODULE OF IDENTITIES AND THE LOW-
#  DIMENSIONAL HOMOLOGY OF A FINITELY PRESENTED GROUP

#                      (Written by Irina Kholodna)



#  This file contains GAP V.4 functions for computing the module of
#  identities of a finite presentation of a finite group G, and for
#  computing a free ZG-resolution. 



#  For all the functions:
#  by  **the list of records representing an element zg in 
#          ZG+...+ZG (k copies of ZG)** 
#  we mean the constraction that becomes clear if we consider 
#  the following example.
#  Let G has a presentation  <x, y| x^2, y^2, (x*y)^2>.
#  The command Elements(G) gives the list [1, x, y, xy].
#  Let zg be an element in the free module ZG+ZG, e.g.
#  zg=(x+y+2xy)e1+(x+3y)e2.
#  **The list of records representing this element**  is:
#       [ [ rec(coeff:=1, word:=x),
#           rec(coeff:=1, word:=y),
#           rec(coeff:=2, word:=xy) ],
#         [ rec(coeff:=1, word:=x,
#           rec(coeff:=3, word:=y)  ]   ]





#_____________________GLOBAL__VARIABLE__________________________

table:=[];
     
#---------------------------------------------------------------    


#***************************************************************
Table:=function(G)
#***************************************************************

#   Given a finite group G, this function constracts
#   a square s x s matrix  (where s = Size(G)) in the following way.
#   Let g_i = Elements(G)[i], and [n_1, n_2, ..., n_s] is the
#   list of integers representing an element z in ZG in the usual
#   way. Then the ith row  [m_1, m_2, ..., m_s] of the matrix 
#   is a list of integers such that [n_(m_1), n_(m_2), ..., n_(m_s)]
#   represents g_i*z.


local #pG, 
      product, permutation, index,
      i, z, pg;

table:=[];
#pG:=Image(IsomorphismPermGroup(G));

for pg in Elements(pG) do
  product:=[];
  for z in Elements(pG) do
    Add(product,Position(Elements(pG), pg*z));
  od;
  permutation:=MappingPermListList(product,[1..Size(pG)]);
  index:=[];
  for i in [1..Size(pG)] do
    Add(index,i^permutation);
  od;
  Add(table,index);
od;

end;



#***************************************************************
NullList:=function(k)
#***************************************************************

local i, list;

  list:=[];
for i in [1..k] do
  list[i]:=0;
od;

return list;

end;


 
#***************************************************************
LastFactor:=function(G,w)
#***************************************************************

#   Given a finitely presented group G and a word w
#   in terms of generators of G, this function returns
#   the extreme right letter in this word.
 
local z, d;

for z in FreeGeneratorsOfFpGroup(G) do
  d:=w*z^-1;
  if  d<w 
  then  return z;
  else  d:=w*z;
        if  d<w 
        then  return z^-1;
        fi;
  fi;
od;

end;


#***************************************************************
WordsToIntegers:=function(G,zg)
#***************************************************************

#   Given a finite group G and the list of records zg (IN FREE GENERATORS!)
#   representing an element in ZG, this
#   function returns the integer vector of length |G|
#   corresponding to this element.

local #FpToPerm, #pG,
      position, l,
      d, i;


l:=NullList(Size(pG));
for d in zg do
  d.word:=MappedWord(d.word,FreeGeneratorsOfFpGroup(G),GeneratorsOfGroup(G));
  position:=Position(Elements(pG),d.word^FpToPerm);
  l[position]:=l[position]+d.coeff;
  od;
 

return l;

end;


#***************************************************************
Der:=function(G,w,i)
#***************************************************************

#   Given a finitely presented group G, a word w in terms of
#   freee generators of F:=FreeGroupOfFpGroup(G) and integer number i 
#   (0 < i < |FreeGeneratorsOfFpGroup(G)|+1), this function returns an
#   integer list representing (as an element in ZG) the Fox derivative 
#   of w with respect to the generator of F whose position in the list
#   of generators is i.
  
local F,
      d, rem, apply;

F:=FreeGroupOfFpGroup(G);

apply:=function(F,w,i,rem)

local l, next, item;

if Length(w)=1 
then if  w<>GeneratorsOfGroup(F)[i] and w<>GeneratorsOfGroup(F)[i]^-1  
     then  return rec(coeff:=0, word:=[]);
     else  if w=GeneratorsOfGroup(F)[i]
           then return rec(coeff:=1, word:=One(F));
           else return rec(coeff:=-1, word:=w);
           fi; 
     fi;
fi;

l:=LastFactor(F,w);
w:=w*l^-1;
item:=rec(coeff:=apply(F,l,i,rem).coeff,
            word:=w*apply(F,l,i,rem).word);
if item.word<>[] 
then  Add(rem, item); 
fi;

next:=apply(F,w,i,rem);

return next;

end;


rem:=[];
d:=apply(F,w,i,rem);

if  d.word<>[] 
then  Add(rem,d); 
fi;

return WordsToIntegers(G,rem);

end;


#***************************************************************
IntegersToWords:=function(G,v)
#***************************************************************

#   Given a finite group G and an integer vector v whose
#   length |v| is an integer multiple of the order of G,
#   this function constract the list of records representing
#   the corresponding (in the usual way) element in
#   ZG+...+ZG (|v|/|G| copies).
 
local zg, 
      k, l,
      i, j;



k:=Length(v)/Size(G);
zg:=[];
for i in [1..k] do
  zg[i]:=[];
  l:=0;
  for j in [Size(G)*(i-1)+1..Size(G)*i] do
    l:=l+1;
    if v[j]<>0 
    then Add(zg[i], rec(coeff:=v[j],
                          word:=EltsG[l]));
    fi;
  od;
od;

return zg;

end;




#***************************************************************
VectorPermutations:=function(G, v)
#***************************************************************

#   Given a finite group G and an integer
#   vector v whose length |v| is an integer multiple of the
#   order of G, this function returns a list of integer
#   vectors U={u_i | i=1,...,|G|} by permuting the coordinates
#   of v as follows. The vector v is identified in the usual
#   way with an element in the free module ZG+...+ZG 
#   (|v|/|G| copies of ZG). The vector u_i represents the 
#   element g_i*v in this module (where {g_i | i=1,...,|G|}
#   is the set of elements of G.

local #pG,
      k, U, 
      Index, i, j;

#pG:=Image(IsomorphismPermGroup(G));
k:=Size(v)/Size(pG);
U:=[];

for j in [1..Size(pG)] do
  
  Index:=StructuralCopy(table[j]);
  for i in [1..k-1] do 
    Append(Index, table[j]+i*Size(pG));
  od;
  Add(U,v{Index});

od; 

return U;

end;

#***************************************************************
GenerateAsZGModule:=function(G,V)
#*************************************************************** 

local products, v;

products:=[];
for v in V do
  Append(products, VectorPermutations(G, v));
od;

return products;

end;

         
#***************************************************************
FoxMatrix:=function(G)
#***************************************************************

#  Given a finitely presented finite group G=<X|R>, the command 
#                   A:=FoxMat(G);
#  constructs the integer matrix A representing the homomorphism
#                   delta: C2 --> C1,
#  where C2 is the direct sum  ZG+ZG+...+ZG  (|R| copies of ZG),
#        C1 is the direct sum  ZG+ZG+...+ZG  (|X| copies of ZG).


local #pG,
      der, row, r, i;

Table(G);

der:=[];
for r in RelatorsOfFpGroup(G) do
  row:=[];
  for i in [1..Length(FreeGeneratorsOfFpGroup(G))] do
    Append(row, Der(G, r, i));
  od;  
  Add(der,row);
od;

return GenerateAsZGModule(G,der);

end;
  



#***************************************************************
Is_ZLinearComb:=function(list, vector)
#***************************************************************

#   Given a list of integer vectors of equal length 
#   and an integer vector of the same length, this
#   function returns true if the vector is Z-linear
#   combination of vectors in the list and false
#   otherwise.


local B,
      list_and_vector, reduced_list,
      space, basis, coeff,
      normal_list_and_vector, normal_list;

B:=false;
if list<>[]
  then
    list_and_vector:=StructuralCopy(list);
    Add(list_and_vector,vector);
    if  RankMat(list_and_vector) = RankMat(list)
    then normal_list_and_vector:=NormalFormIntMat(list_and_vector,2);
         normal_list:=NormalFormIntMat(list,2);
         if 
      normal_list_and_vector.normal{[1..normal_list_and_vector.rank]}=
      normal_list.normal{[1..normal_list.rank]}
         then B:=true;
         fi;
    fi;

fi;

return B;

end;


#***************************************************************
MinimalModuleGenerators:=function(G, S)
#***************************************************************

#   Given a finitely presented finite group G=<X|R> and 
#   a set S = {v_1,...,v_n} of integer vectors of equal length 
#   k=|v_i| with k an integer multiply of |G|, the command
#              T:=MinimalModuleGenerators(G,S)
#   consracts a subset T of S as follows. 
#   The vector v_i are identified with elements in the free 
#   module ZG+...+ZG (k/|G| copies of ZG). Let A denote the 
#   abelian group consisting of all finite integer linear 
#   combinations of the elements in S. Let B denote the 
#   ZG-module generated by the elements of T. Then T is 
#   constructed to have the property A < B; moreover, no subset 
#   of T has this property.


local T, 
      temp, G_temp,
      sublist_G_temp, except, index,
      i, j;

temp:=[];
G_temp:=[];
for i in [1..Length(S)] do
  if not (S[i] in G_temp) 
  then  if not Is_ZLinearComb(G_temp,S[i]) 
        then Add(temp, S[i]);
             Append(G_temp, VectorPermutations(G,S[i]));
        fi;
  fi;
od;



T:=[];
except:=[];
for i in [1..Length(temp)] do

  index:=[];
  for j in [1..Length(temp)] do
    if  not(j in except) and i<>j 
    then  Append(index, [Size(G)*(j-1)+1..Size(G)*j]);
    fi;
  od;

  sublist_G_temp:=G_temp{index};
  if  Is_ZLinearComb(sublist_G_temp, temp[i])
  then  Add(except,i);
  else  Add(T, temp[i]);
  fi;

od;

return T;


end;


#***************************************************************
IdentityModuleZBasis:=function(G)
#***************************************************************

#   Given a finitely presented finite group G=<X|R>, the command 
#               S:=IdentityModuleZBasis(G);
#   constract a set S of vectors that form a basis for the free 
#   abelian group underlying the module of identities \pi.
#   Thus the vectors in S have length equal to |G|*|R|;
#   the implicit ordering on the element of G is that given by
#   the standart GAP command
#                      Elements(G);
#   for the listing the elements of G.

local A, V, v,nf;



A:=FoxMatrix(G);

nf:=NormalFormIntMat(A,6);
V:=nf.rowtrans{[nf.rank+1..Length(A)]};

return LLLReducedBasis(V).basis;
#return BaseIntMat(V);

end;


#***************************************************************
IdentityModuleGenerators:=function(G)
#***************************************************************

#   Given a finitely presented finite group G=<X|R>, the command 
#               T:=IdentityModuleGenerators(G)
#   constracts a set T of integer vectors v_i of length
#   |v_i|=|G|*|R|. The set T represents a minimal set of 
#   generators  for the module of identities \pi.


local m;



m:= MinimalModuleGenerators(G, IdentityModuleZBasis(G));


return m;

end;

#***************************************************************
Dimension:=function(n);
if n=0 then return 1; fi;
if n=1 then return Length(GeneratorsOfGroup(G)); fi;
if n= 2 then return Length(RelatorsOfFpGroup(G)); fi;
if n>2 then return Length(IrinasRes[n]); fi;
end;
#*

 
#***************************************************************
Resolution:=function(G,n)
#***************************************************************

#   Given a finitely presented finite group G=<X|R>, 
#   and an integer number n>2, the command
#          RES:=Resolution(G,n);
#   construct a list RES such that for 2<i<n
#   RES[i] is a minimal set of generators
#   for the i-th term of a free ZG-resolution of Z. 


local res, matrix,
      nf, V,
      mgen, i, j;

Table(G);
res:=[];
res[1]:=[];
for i in [1..Dimension(1)] do
V:=[];
        for j in [1..Dimension(1)*OrderG] do
        if j=(1+(i-1)*OrderG) then V[j]:=-1;else

        if j=Position(EltsG,GeneratorsOfGroup(G)[i])+(i-1)*Order(G) then
        V[j]:=1;
        else V[j]:=0;
        fi; fi;
        od;
V:=List([1..Order(G)],x->V[x+(i-1)*Order(G)]);
Append(res[1],[V]);
od;

if not IsFpGroup(G) then
      matrix:=[];
      for mgen in res[1] do
        Add(matrix, mgen);
      od;
      matrix:=GenerateAsZGModule(G,matrix);
      nf:=NormalFormIntMat(matrix,6);
      V:=nf.rowtrans{[nf.rank+1..Length(matrix)]};
      V:=LLLReducedBasis(V).basis;
      res[2]:=MinimalModuleGenerators(G,V);
fi;

for i in [3..n] do

if IsFpGroup(G) and i=3 
then  res[i]:=IdentityModuleGenerators(G);
else  
      matrix:=[];
      for mgen in res[i-1] do
        Add(matrix, mgen);
      od;
      matrix:=GenerateAsZGModule(G,matrix);
       
      nf:=NormalFormIntMat(matrix,6);
      V:=nf.rowtrans{[nf.rank+1..Length(matrix)]};
      V:=LLLReducedBasis(V).basis;
      #V:=BaseIntMat(V);
      res[i]:=MinimalModuleGenerators(G,V);
fi;

od;

return res;

end;
   
#####################################################################   
########################END OF IRINA'S CODE##########################




IrinasRes:=Resolution(G,n);

#####################################################################
Dimension:=function(n);
if n=0 then return 1; fi;
if n>2 or (not IsFpGroup(G)) then return Length(IrinasRes[n]); fi;
if n=1 then return Length(GeneratorsOfGroup(G)); fi;
if n=2 then return Length(RelatorsOfFpGroup(G)); fi;
end;
#####################################################################


#####################################################################
ExtendIrinasRes:=function()	#This function will add the first
local i, j, V;			#two dimensions to Irinas Resolution.

IrinasRes[1]:=[];
IrinasRes[2]:=[];

if IsFpGroup(G) then
for i in [1..Dimension(2)] do
V:=[];
	for j in [1..Dimension(2)*OrderG] do
	if j=(1+(i-1)*OrderG) then V[j]:=-1; else V[j]:=0; fi;
	od;
Append(IrinasRes[2], [TransposedMat(FoxM)*V]);
od;
fi;

for i in [1..Dimension(1)] do
V:=[];
	for j in [1..Dimension(1)*OrderG] do
	if j=(1+(i-1)*OrderG) then V[j]:=-1;else
	
	if j=Position(EltsG,GeneratorsOfGroup(G)[i])+(i-1)*Order(G) then 
	V[j]:=1;
	else V[j]:=0;
	fi; fi;
	od;
V:=List([1..Order(G)],x->V[x+(i-1)*Order(G)]);
Append(IrinasRes[1],[V]);
od;

end;
#####################################################################

if IsFpGroup(G) then  FoxM:=FoxMatrix(G); ExtendIrinasRes(); fi;


#####################################################################
Int2Pair:=function(i,m)    	#m=Order(G), i=Position in a vector
local j, x;

j:=i mod m;
x:=(i-j)/m;

if j=0 then return [x,m];
else
return [x+1,j];
fi;

end;
#####################################################################

#####################################################################
Boundary:=function(n,kk)
local B, i, j, p, w, k, StandardForm;

if n<1 then return []; fi;

k:=AbsoluteValue(kk);
B:=IrinasRes[n][k];
StandardForm:=[];

for i in [1..Length(B)] do
if not B[i]=0 then
   for j in [1..AbsoluteValue(B[i])] do
   p:=Int2Pair(i,OrderG);
   Append(StandardForm, [   [SignInt(B[i])*p[1],AbsoluteValue(p[2])]   ]);
   od;
fi;
od;

if kk>0 then return StandardForm;
else
return NegateWord(StandardForm);
fi;
end;
#####################################################################



return Objectify(HapResolution,
	    rec(
	    dimension:=Dimension,
	    boundary:=Boundary,
	    homotopy:=fail,
	    elts:=EltsG,
	    group:=G,
	    properties:=
	     [["type","resolution"],
	     ["length",n],
	     ["characteristic", Charact],
	     ["reduced",true] ]));
end);
#####################################################################
