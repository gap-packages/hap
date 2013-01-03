#(C) Graham Ellis, 2005-2006


#####################################################################
InstallGlobalFunction(ResolutionDirectProduct,
function(arg)
local
	R,S,T,
	G,H,E,K,GhomE,HhomE,EhomG,EhomH,EltsE,
	pcpEmod,Emod,Ehom, homE,
 	DimensionR,BoundaryR,HomotopyR,
        DimensionS,BoundaryS,HomotopyS,
	Lngth,Dimension,Boundary,Homotopy,
	PseudoBoundary,
	DimPQ,
	Int2Pair, Pair2Int,
	Charact,
	AddWrds,
	Int2Vector,
	Vector2Int,
	Elts2Int,
	HomotopyRec,
	HomotopyGradedGen,
	HomotopyOfWord,
	FinalHomotopy,
	HorizontalBoundaryGen,
	HorizontalBoundaryWord,
	F,FhomE,
	gensE, gensE1, gensE2, 
	AppendToElts, RappendToElts, SappendToElts,
	Boole,
	i,j,k,g,h,fn;

R:=arg[1];
S:=arg[2];

if "appendToElts" in NamesOfComponents(R) then RappendToElts:=R!.appendToElts;
else
RappendToElts:=function(x); Append(R!.elts,[x]);  end;
fi;

if "appendToElts" in NamesOfComponents(S) then SappendToElts:=S!.appendToElts;
else
SappendToElts:=function(x); Append(S!.elts,[x]);  end;
fi;


G:=R!.group; 
H:=S!.group; 

####################### DIRECT PRODUCT OF GROUPS ###########
if Length(arg)=2 then

E:=DirectProduct(G,H);
if Size(G)=infinity or Size(H)=infinity then SetSize(E,infinity);fi;
GhomE:=Embedding(E,1);
HhomE:=Embedding(E,2);
EhomG:=Projection(E,1);
EhomH:=Projection(E,2);

else  	#if G and H both lie in a group K, and if they commute and have
	#have trivial intersection then we create their direct product as
	#a subgroup of K. We treat pcp groups as a seperate case.


#####PCP CASE #######################
if IsPcpGroup(G)   
then
K:=PcpGroupByCollector(Collector(Identity(G)));
#G:=Subgroup(K,GeneratorsOfGroup(G));
#H:=Subgroup(K,GeneratorsOfGroup(H));

gensE:=(Concatenation(GeneratorsOfGroup(G),GeneratorsOfGroup(H)));
E:=Subgroup(K,gensE);
       
GhomE:=GroupHomomorphismByFunction(G,E,x->x);
HhomE:=GroupHomomorphismByFunction(H,E,x->x);


EhomG:=GroupHomomorphismByImagesNC(E,G,
GeneratorsOfGroup(E),
Concatenation(GeneratorsOfGroup(G),
List(GeneratorsOfGroup(H),x->One(G))));

EhomH:=GroupHomomorphismByImagesNC(E,H,
GeneratorsOfGroup(E),
Concatenation(
List(GeneratorsOfGroup(G),x->One(G)),
GeneratorsOfGroup(H)));


fi;
############PCP CASE DONE ###########

############OTHER CASE##############
if not IsPcpGroup(G) then
gensE:=Concatenation(GeneratorsOfGroup(G),GeneratorsOfGroup(H));
E:=Group(gensE);

gensE1:=Concatenation(GeneratorsOfGroup(G),
 List([1..Length(GeneratorsOfGroup(H))],x->Identity(G)));
 gensE2:=Concatenation(List([1..Length(GeneratorsOfGroup(G))],x->Identity(H)),
              GeneratorsOfGroup(H));

GhomE:=GroupHomomorphismByFunction(G,E,x->x);
HhomE:=GroupHomomorphismByFunction(H,E,x->x);
EhomG:=GroupHomomorphismByImagesNC(E,G,gensE,gensE1);
EhomH:=GroupHomomorphismByImagesNC(E,H,gensE,gensE2);
fi;
###########OTHER CASE DONE#########

fi;
################ DIRECT PRODUCT OF GROUPS CONSTRUCTED #########



EltsE:=[Identity(E)];

#####################################################################
	AppendToElts:=function(x);
	EltsE[Length(EltsE)+1]:=x;
	end;
#####################################################################

PseudoBoundary:=[];
DimensionR:=R!.dimension; 
DimensionS:=S!.dimension; 
BoundaryS:= S!.boundary;
	   
BoundaryR:=R!.boundary;  
HomotopyR:=R!.homotopy;
HomotopyS:=S!.homotopy;  

#################DETERMINE VARIOUS PROPERTIES########################
Lngth:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")>0
then Charact:=EvaluateProperty(S,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")>0
then Charact:=Product(Intersection([
DivisorsInt(EvaluateProperty(R,"characteristic")),
DivisorsInt(EvaluateProperty(S,"characteristic"))
]));
fi;

if Charact=0 then AddWrds:=AddFreeWords; else
        AddWrds:=function(v,w);
        return AddFreeWordsModP(v,w,Charact);
        end;
fi;
####################PROPERTIES DETERMINED############################


#####################################################################
Dimension:=function(i)
local D,j;
if i<0 then return 0; fi;
if i=0 then return 1; else
D:=0;

for j in [0..i] do
D:=D+DimensionR(j)*DimensionS(i-j);
od;

return D; fi;
end;
#####################################################################

for i in [1..Lngth] do
PseudoBoundary[i]:=[1..Dimension(i)];
od;

#####################################################################
DimPQ:=function(p,q)
local D,j;

if (p<0) or (q<0) then return 0; fi;

D:=0;
for j in [0..q] do
D:=D+DimensionR(p+q-j)*DimensionS(j);
od;

return D;
end;
#####################################################################

#####################################################################
Int2Pair:=function(i,p,q)       #Assume that x<=DimR(p)*DimS(q).
local s,r,x;
                                #The idea is that the generator f_i in F
				#corresponds to a tensor (e_r x e_s)
x:=AbsInt(i)-DimPQ(p+1,q-1);     #with e_r in R_p, e_s in S_q. If we
s:= x mod DimensionS(q);                #input i we get output [r,s].
r:=(x-s)/DimensionS(q);

if s=0 then return [SignInt(i)*r,DimensionS(q)];
else return [SignInt(i)*(r+1),s]; fi;

end;
#####################################################################

#####################################################################
Pair2Int:=function(x,p,q)
local y;                        #Pair2Int is the inverse of Int2Pair.

y:=[AbsInt(x[1]),AbsInt(x[2])];
return SignInt(x[1])*SignInt(x[2])*((y[1]-1)*DimensionS(q)+y[2]+DimPQ(p+1,q-1));end;
#####################################################################

#####################################################################
Int2Vector:=function(k,j)
local tmp,p,q;

p:=k;q:=0;
while j>=DimPQ(p,q)+1 do
p:=p-1;q:=q+1;
od;				#p,q are now computed from k,j

tmp:=Int2Pair(j,p,q);

return [p,q,tmp[1],tmp[2]];
end;
#####################################################################

#####################################################################
Vector2Int:=function(p,q,r,s);
return Pair2Int([r,s],p,q);
end;
#####################################################################

T:=0;
#####################################################################
Elts2Int:=function(x)
local pos;

pos:=Position(EltsE,x);
if IsPosInt(pos) then return pos;
else
	Append(EltsE,[x]);
	if not IsInt(T) then
	Append(T!.elts,[x]); fi;
	return Length(EltsE); 
fi;
end;
#####################################################################

#####################################################################
Boundary:=function(k,jj)
local j, p,q,r,s,tmp, horizontal, vertical;

if k<1 then return []; fi;
j:=AbsInt(jj);
	#################IF BOUNDARY NOT ALREADY COMPUTED############
if IsInt(PseudoBoundary[k][j]) then
tmp:=Int2Vector(k,j);
p:=tmp[1]; q:=tmp[2]; r:=tmp[3]; s:=tmp[4];

horizontal:=ShallowCopy(BoundaryR(p,r));
Apply(horizontal,x->[x[1],Elts2Int(   Image(GhomE,R!.elts[x[2]])   )  ]);
Apply(horizontal,x->[Vector2Int(p-1,q,x[1],s),x[2]]);

vertical:=ShallowCopy(BoundaryS(q,s));
Apply(vertical,x->[x[1],Elts2Int(   Image(HhomE,S!.elts[x[2]])  )    ]);
Apply(vertical,x->[Vector2Int(p,q-1,r,x[1]),x[2]]);
if IsOddInt(p) then
vertical:=NegateWord(vertical);
fi;


PseudoBoundary[k][j]:= Concatenation(horizontal, vertical);
fi;
	################IF ENDS######################################

if SignInt(jj)=1 then
return PseudoBoundary[k][j];
else
return NegateWord(PseudoBoundary[k][j]);
fi;
end;
#####################################################################

#####################################################################
HorizontalBoundaryGen:=function(n,y)
local a,i, p,q,r,s, tmp,horizontal;

a:=AbsInt(y[1]);
tmp:=Int2Vector(n,a);

p:=tmp[1]; q:=tmp[2]; r:=tmp[3]; s:=tmp[4];

horizontal:=StructuralCopy(BoundaryR(p,r));

Apply(horizontal,x->[x[1],Elts2Int( EltsE[y[2]]*Image(GhomE,R!.elts[x[2]]) )]);


Apply(horizontal,x->[Vector2Int(p-1,q,x[1],s),x[2]]);

return horizontal;

end;
#####################################################################

#####################################################################
HorizontalBoundaryWord:=function(n,w)
local x, bnd;

bnd:=[];
for x in w do
Append(bnd,HorizontalBoundaryGen(n,x));
od;
return bnd;

end;
#####################################################################

#####################################################################
HomotopyGradedGen:=function(g,p,q,r,s,bool)    #Assume EltsE[g] exists!
local aa,hty, hty1, Eg, Eg1, Eg2, g1, g2;	#bool=true for vertical homotopy

#This function seems to work! But I should really check the maths again!!

Eg:=EltsE[g];
Eg1:=Image(EhomG,Eg);
	if not Eg1 in R!.elts then RappendToElts(Eg1); fi; 
				    
Eg2:=Image(EhomH,Eg);
	if not Eg2 in S!.elts then SappendToElts(Eg2); fi;
				    
g2:=Position(S!.elts,Eg2); 
g1:=Position(R!.elts,Eg1); 
Eg1:=Image(GhomE,Eg1);
Eg2:=Image(HhomE,Eg2);


hty:=HomotopyS(q,[s,g2]);
Apply(hty,x->[ Vector2Int(p,q+1,r,x[1]), Image(HhomE,S!.elts[x[2]])]); 
Apply(hty,x->[ x[1], Elts2Int(Eg1*x[2])]);
if IsOddInt(p) then
hty:=NegateWord(hty); fi;

if (p=0 and q>0) or bool then return hty; fi;

if p>0 then
hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty)),false);
Append(hty, NegateWord(hty1));
fi;

if q>0 then return hty; fi;


hty1:=HomotopyR(p,[r,g1]);
Apply(hty1,x->[ Vector2Int(p+1,q,x[1],s), Image(GhomE,R!.elts[x[2]])]);
Apply(hty1,x->[ x[1], Elts2Int(x[2])]); #Here

Append(hty,hty1);

hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty1)),true);

Append(hty,NegateWord(hty1));

hty1:=HomotopyOfWord(p+q,StructuralCopy(HorizontalBoundaryWord(p+q+1,hty1)),false);


Append(hty,hty1); 	#I think this perturbation term is always zero and
			#thus not necessary.

return hty;
end;
#####################################################################

#####################################################################
Homotopy:=function(n,x,bool)
local vec,a;


a:=AbsInt(x[1]);
vec:=Int2Vector(n,a);
if SignInt(x[1])=1 then 
return HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool);
else
return NegateWord(HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool));
fi;
end;
#####################################################################

#####################################################################
HomotopyOfWord:=function(n,w,bool)
local x, hty;

hty:=[];
for x in w do
Append(hty,Homotopy(n,x,bool));
od;
return hty;

end;
#####################################################################

HomotopyRec:=[];
for i in [1..Lngth] do
HomotopyRec[i]:=[];
for j in [1..Dimension(i-1)] do
HomotopyRec[i][j]:=[];
od;od;

#####################################################################
FinalHomotopy:=function(n,x)
local a;

a:=AbsInt(x[1]);
if not IsBound(HomotopyRec[n+1][a][x[2]]) then
HomotopyRec[n+1][a][x[2]]:=Homotopy(n,[a,x[2]],false);
fi;

if SignInt(x[1])=1 then
return StructuralCopy(HomotopyRec[n+1][a][x[2]]);
else
return NegateWord(StructuralCopy(HomotopyRec[n+1][a][x[2]]));
fi;

end;
#####################################################################

if HomotopyR=fail or HomotopyS=fail then
FinalHomotopy:=fail;
fi;

for i in [1..Lngth] do
for j in [1..Dimension(i)] do
g:=Boundary(i,j);
od;
od;

#####################################################################
AppendToElts:=function(x)
Append(T!.elts,[x]);
end;
#####################################################################

Boole:=false;
if EvaluateProperty(R,"reduced")=true
and  EvaluateProperty(S,"reduced")=true
then Boole:=true;
fi;


T:= rec(
            dimension:=Dimension,
	    boundary:=Boundary,
	    homotopy:=FinalHomotopy,
	    elts:=EltsE,
	    group:=E,
	    firstProjection:=EhomG,
	    secondProjection:=EhomH,
	    properties:=
	    [["type","resolution"],
	    ["length",Lngth],
	    ["characteristic",Charact],
	    ["reduced",true]],
	    appendToElts:=AppendToElts);


return Objectify(HapResolution,T);
end);

#####################################################################

