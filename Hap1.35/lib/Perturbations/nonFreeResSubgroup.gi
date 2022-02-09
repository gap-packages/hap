#(C) Graham Ellis, 2005-2006

######################################################
#######################################################
InstallGlobalFunction(RelativeRightTransversal,
function(G,H,U)
local T, fn, R, RR, poscan;
T:=RightTransversal(G,H);;

####################
fn:=function(g)
local u,L;
L:=[];;
for u in U do
Add(L,PositionCanonical(T,g*u));
od;
return Minimum(L);
end;;
####################

R:=Classify(T,fn);;
RR:=List(R,x->x[1]);


####################
poscan:=function(g)
local k, x;
x:=T[PositionCanonical(T,g)];
for k in [1..Length(R)] do
if x in R[k] then return k; fi;
od;
end;
####################

return Objectify( NewType( FamilyObj( G ),
                    IsHapRightTransversalSL2ZSubgroup and IsList and
                    IsDuplicateFreeList and IsAttributeStoringRep ),
          rec( group := G,
               subgroup := H,
               cosets:=RR,
               poscan:=poscan ));
end);
#######################################################
#######################################################

#####################################################################
InstallGlobalFunction(NonFreeResolutionFiniteSubgroup,
function(arg)
local 
		R,gensG,gensK,
		DimensionR, BoundaryR, HomotopyR, EltsG,
		Dimension, Boundary, BoundaryRec, Homotopy, EltsK,
		G, K, TransK, sK, 
		Gword2Kword, G2K, G2KRec, Pair2Int, Int2Pair,
		Mult, MultRec, FIN, i, Stabilizer, Action,
                StabRec, Conj, ConjRec, n, k, len;
		
if Length(arg)=3 then
R:=arg[1]; gensG:=arg[2]; gensK:=arg[3];
fi;
if Length(arg)=2 then
R:=arg[1]; gensG:=R!.group; gensK:=arg[2];
fi;
				#gensG and gensK originally had to be
				#generating sets. Later I allowed them to be
				#the groups G, K themselves. Sloppy!
DimensionR:=R!.dimension;
BoundaryR:=R!.boundary;
HomotopyR:=R!.homotopy;
EltsG:=R!.elts;

if IsList(gensG) then G:=Group(gensG); else G:=gensG; fi;
if IsList(gensK) then K:=Group(gensK); else K:=gensK; fi;
if Size(G)=Length(EltsG) then FIN:=true; else FIN:= false; fi;
if FIN then
   if IsPseudoList(EltsG) then EltsK:=EltsG;
   else
   EltsK:=Elements(K);
   fi;
else
EltsK:=[];
fi;


len:=0;
while DimensionR(len+1)>0 do len:=len+1; od;
TransK:=[];
sK:=[];
for n in [0..len] do
TransK[n+1]:=[];
for k in [1..DimensionR(n)] do
TransK[n+1][k]:=RelativeRightTransversal(G,K,R!.stabilizer(n,k));
od;
sK[n+1]:=Sum(List(TransK[n+1],x->Length(x)));
od;




if FIN then
#####################################################################
Mult:=function(n,k,i,j);
return Position(EltsG,TransK[n+1][k][i]*EltsG[j]);
end;
#####################################################################
else
#####################################################################
Mult:=function(n,k,i,j) local x,r;

x:=TransK[n+1][k][i]*EltsG[j];

r:=Position(EltsG,x);

if r=fail then Add(EltsG,x); r:=Length(EltsG); fi;
return r;
end;
#####################################################################
fi;

if FIN then
#####################################################################
Conj:=function(n,k,i,j);
return Position(EltsG,TransK[n+1][AbsInt(k)][i]^-1*EltsG[j]*TransK[n+1][AbsInt(k)][i]);
end;
#####################################################################
else
#####################################################################
Conj:=function(n,k,i,j) local x,r;

x:=TransK[n+1][AbsInt(k)][i]^-1*EltsG[j]*TransK[n+1][AbsInt(k)][i];

r:=Position(EltsG,x);

if r=fail then Add(EltsG,x); r:=Length(EltsG); fi;
return r;
end;
#####################################################################
fi;



#####################################################################
Dimension:=function(n);
if n>len then return 0; fi;
return sK[n+1];
end;
#####################################################################



if FIN then
#####################################################################
G2K:=function(n,k,g)
local t,kk;
t:=PositionCanonical(TransK[n+1][k],EltsG[g]);
kk:=Position(EltsK,EltsG[g]*TransK[n+1][k][t]^-1);
return [kk,t];
end;
#####################################################################
else
#####################################################################
G2K:=function(n,k,g)
local t,kk,x;

t:=PositionCanonical(TransK[n+1][k],EltsG[g]);
x:=EltsG[g]*TransK[n+1][k][t]^-1;

kk:=Position(EltsK,x);

if kk=fail then Add(EltsK,x); kk:=Length(EltsK); fi;

return [kk,t];

end;
#####################################################################
fi;

#####################################################################
Pair2Int:=function(n,x)
local i,t;
i:=x[1]; t:=x[2];
return SignInt(i)*((AbsoluteValue(i)-1)*sK[n+1] + t);
end;
#####################################################################

#####################################################################
Int2Pair:=function(n,i)
local j,k, x;
j:=AbsoluteValue(i);
x:=j mod sK[n+1];
k:=(j-x)/sK[n+1];
if not x=0 then return [SignInt(i)*(k+1),x]; else
return [SignInt(i)*k,sK[n+1]]; fi;
end;
#####################################################################

#####################################################################
Gword2Kword:=function(n,w)
local x, y, v;

v:=[];
for x in w do
y:=G2K(n,AbsInt(x[1]),x[2]);
y:=[Pair2Int(n,[x[1],y[2]]),y[1]];
Add(v,y);
od;
return v;
end;
#####################################################################

BoundaryRec:=[];
for i in [1..EvaluateProperty(R,"length")] do
BoundaryRec[i]:=[];
od;

#####################################################################
Boundary:=function(n,ii)
local x, w, i;

if n<=0 then return []; fi;

i:=AbsInt(ii);
if not IsBound(BoundaryRec[n][i]) then
x:=Int2Pair(n,i);
w:=StructuralCopy(BoundaryR(n,x[1]));
#Apply(w, y->[y[1],Mult(n,x[1],x[2],y[2])]);
Apply(w, y->[R!.action(n-1,AbsInt(y[1]),x[2])*y[1],Mult(n,x[1],x[2],y[2])]);
BoundaryRec[n][i]:= Gword2Kword(n-1,w);
fi;

if ii>0 then return BoundaryRec[n][i]; 
else return NegateWord(BoundaryRec[n][i]); fi;
end;
#####################################################################


######################
Homotopy:=function(n,e)
local x,g,pos;

x:=Int2Pair(n,e[1]);
g:=EltsK[e[2]]*TransK[x[2]]; #Need to chack the maths again here!


pos:=Position(EltsG,g);


if pos =fail then Add(EltsG,g); pos:=Length(EltsG);fi;
return Gword2Kword( R!.homotopy(n,[x[1],pos]));

end;
######################

StabRec:=[];
######################
Stabilizer:=function(n,k)
local i,t;
if not IsBound(StabRec[n+1]) then
StabRec[n+1]:=[];
for i in [1..R!.dimension(n)] do
for t in TransK[n+1][k] do
Add(StabRec[n+1],Intersection(R!.stabilizer(n,i)^(t^-1),gensK));
od;
od;
fi;

return StabRec[n+1][k];

end;
######################

######################
Action:=function(n,k,g)
local x, c, p, gg; 
x:=Int2Pair(n,k);
gg:=Position(EltsG,EltsK[g]);
if gg=fail then Add(EltsG,EltsK[g]); gg:=Length(EltsG); fi;
c:=Conj(n,x[1],x[2],gg);;
#p:=Position(R!.elts,EltsG[c]);
#if p=fail then Add(R!.elts, EltsG[c]); p:=Length(R!.elts); fi;
return R!.action(n,x[1],c);
end;
######################

return Objectify(HapNonFreeResolution,
	     rec(
	     dimension:=Dimension,
	     boundary:=Boundary,
	     homotopy:=Homotopy,
             stabilizer:=Stabilizer,
             action:=Action,
	     elts:=EltsK,
	     group:=K,
             Int2Pair:=Int2Pair,
	     properties:=
	     [["length",EvaluateProperty(R,"length")],
	      ["characteristic",EvaluateProperty(R,"characteristic")],
	      ["reduced",false],
	      ["type","resolution"] ]));
end);
#####################################################################


