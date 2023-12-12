#(C)2009 Graham Ellis

#########################################################
InstallGlobalFunction(HAPTietzeReduction_OneStep,
function(R,N,bound)
local
	CR,
        Dimension,
        Boundary,
        Homotopy,
        HomotopyN,
        HomotopyNplus1,
        HomotopyNminus1,
        PseudoBoundary,
        PseudoBoundaryN,
        PseudoBoundaryNplus2,
	FindFreeFace,
	Action, ActionInv, Elts,
 	modN, modNplus1,
	triple, phi, nwrdd, wrdd,
	tmp, trp3,pos, 
        homotopyrec, homotopyrecpos,
        newb,bool,char,i,j,b,e,g,x,w,D,D2,gg;

char:=EvaluateProperty(R,"characteristic");
Elts:=R!.elts;
modN:=[];
modNplus1:=[];

PseudoBoundary:=1*List([1..R!.dimension(N+1)],i->
                              R!.boundary(N+1,i));
PseudoBoundaryN:=1*List([1..R!.dimension(N)],i->
                              R!.boundary(N,i));
if Length(R)>N+1 then
PseudoBoundaryNplus2:=1*List([1..R!.dimension(N+2)],i->
                              R!.boundary(N+2,i));
else
PseudoBoundaryNplus2:=[];
fi;

#####################################################################


#####################################################################
Action:=function(g,l)
local pos, h;
h:=Elts[g]*Elts[l[2]];
pos:=Position(Elts,h);
if pos=fail then
Add(Elts, h);
pos:=Length(Elts);
fi;
return [1*l[1],pos];
end;
#####################################################################

#####################################################################
ActionInv:=function(g,l)
local pos, h;
h:=Elts[g]^-1*Elts[l[2]];
pos:=Position(Elts,h);
if pos=fail then
Add(Elts, h);
pos:=Length(Elts);
fi;
return [1*l[1],pos];
end;
#####################################################################

if IsPseudoList(Elts) then
if IsBound(Elts!.mult) then
#####################################################################
Action:=function(g,l)
local pos, h;
pos:=Elts!.mult(g,l[2]); 
return [1*l[1],pos];
end;
#####################################################################

#####################################################################
ActionInv:=function(g,l)
local pos, h;
pos:=Elts!.mult(Elts!.inv(g),l[2]);
return [1*l[1],pos];
end;
#####################################################################

fi;fi;
########################################################
FindFreeFace:=function()
local i,b,pos,y, wrd,ee, g,trp3;
#returns either fail, or a triple [e,wrd,i] such that we can
#replace all g multiples of the generating cell [e,1] with g multiples 
#of the word wrd in boundaries of N+1-dimensional cells, and delete the 
#i-th free cell in dimension N+1;.
for i in [1..Length(PseudoBoundary)] do
b:=List(PseudoBoundary[i],x->AbsInt(x[1]));
b:=Collected(b);
pos:=PositionProperty(b,x->x[2]=1);

if IsInt(pos) and Length(b)<=bound then
y:=b[pos][1];
b:=1*PseudoBoundary[i];
pos:=PositionProperty(b,x->AbsInt(x[1]) =y); 
wrd:=b{Concatenation([1..pos-1],[pos+1..Length(b)])};
g:=b[pos][2];
if b[pos][1]>0 then
wrd:=NegateWord(wrd);trp3:=i;
ee:=b[pos][1];
else ee:=-b[pos][1]; trp3:=-i;
fi;
wrd:=List(wrd,x->ActionInv(g,x));
return [ee,wrd,trp3,g];
fi;
od;
return fail;
end;
########################################################

bool:=false;
####################
####################
while true do
if bool then  break; fi;
triple:=FindFreeFace();
if triple=fail then  break; fi;
bool:=true;

e:=1*triple[1];
wrdd:=1*triple[2];
nwrdd:=NegateWord(wrdd);
trp3:=1*triple[3];
gg:=triple[4];
Add(modNplus1,AbsInt(trp3));
Add(modN,e);

#od;
##############################
##############################
#if triple=fail then return R; fi;



###############
for j in [1..Length(PseudoBoundary)] do
b:=1*PseudoBoundary[j];
newb:=[];
for i in [1..Length(b)] do
   if not AbsInt(b[i][1])=e then
   Add(newb,(b[i]));
   else
   w:=1*wrdd;
   if b[i][1]<0 then
   w:=NegateWord(w);
   fi;
   w:=List(w,x->Action(b[i][2],x));
   Append(newb,w);
   fi;
od;
PseudoBoundary[j]:=AlgebraicReduction(newb,char);
od;
######################

modN:=Difference([1..R!.dimension(N)],modN);
modNplus1:=Difference([1..R!.dimension(N+1)],modNplus1);
PseudoBoundary:=PseudoBoundary{modNplus1};


#####################################
for b in PseudoBoundary do
for x in b do
x[1]:=SignInt(x[1])*Position(modN,AbsInt(x[1]));
od;od;
#####################################


PseudoBoundaryN:=PseudoBoundaryN{modN};

#####################################
if Length(R)>N+1 then
for j in [1..Length(PseudoBoundaryNplus2)] do
b:=PseudoBoundaryNplus2[j];
newb:=[];
for x in b do
if  AbsInt(x[1]) in modNplus1 then
Add(newb,[SignInt(x[1])*Position(modNplus1,AbsInt(x[1])),x[2]]);
fi;
od;
PseudoBoundaryNplus2[j]:=newb;
od;
fi;
#####################################

####################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if not i in [N,N+1] then return R!.dimension(i); fi;
if i=N then return  Length(PseudoBoundaryN); fi;
return Length(PseudoBoundary); 
end;
####################################################

####################################################
Boundary:=function(n,i);
if not n in [N,N+1,N+2]  then return R!.boundary(n,i); fi;
if i>0 then
if n=N then return 1*PseudoBoundaryN[i]; fi;
if n=N+1 then return 1*PseudoBoundary[i]; fi;
if n=N+2 then return 1*PseudoBoundaryNplus2[i]; fi;
fi;
if i<0 then
if n=N then return NegateWord(1*PseudoBoundaryN[AbsInt(i)]); fi;
if n=N+1 then return NegateWord(1*PseudoBoundary[AbsInt(i)]); fi;
if n=N+2 then return NegateWord(1*PseudoBoundaryNplus2[AbsInt(i)]); fi;
fi;
end;
####################################################

if not R!.homotopy=fail then
####################################################
HomotopyNplus1:=function(x)
local ans, bnd, w, h;
ans:= 1*R!.homotopy(N+1,[ SignInt(x[1])*modNplus1[AbsInt(x[1])], x[2] ]);
bnd:=R!.boundary(N+1, SignInt(x[1])*modNplus1[AbsInt(x[1])]);

for w in bnd do
if AbsInt(w[1])=e then 
h:=R!.homotopy(N+1, [-SignInt(w[1])*trp3, Position(R!.elts,R!.elts[x[2]]*R!.elts[w[2]]*R!.elts[gg]^-1)   ]);
Append(ans,h);
fi;
od;
return AlgebraicReduction(1*ans, char);
end;
####################################################

####################################################
HomotopyN:=function(x)
local h;
h:=1*R!.homotopy(N,[SignInt(x[1])*modN[AbsInt(x[1])], x[2]]);
h:=Filtered(h,a-> not AbsInt(a[1])=AbsInt(trp3) );
h:=List(h,a -> [ SignInt(a[1])*Position(modNplus1,AbsInt(a[1])), a[2]]);
return h;
end;
####################################################

########################
########################
phi:=function(V)
local v, w, ans;
ans:=[];
for v in V do
   if AbsInt(v[1])=e then
      if v[1]>0 then
         Append(ans, List(wrdd,a->Action(v[2],a))  );
      else
         Append(ans, List(nwrdd,a->Action(v[2],a))  );
      fi;
   else
      Add(ans,v);
   fi;
od;
return ans;
end;
########################
########################


####################################################
HomotopyNminus1:=function(x)
local h;
h:= phi(1*R!.homotopy(N-1,x));
h:=List(h,a -> [ SignInt(a[1])*Position(modN,AbsInt(a[1])), a[2]]);
return AlgebraicReduction(h,char);
end;
####################################################

homotopyrec:=[ ];
homotopyrec[1]:=List([1..Dimension(N-1)],i->[]);
homotopyrec[2]:=List([1..Dimension(N)],i->[]);
homotopyrec[3]:=List([1..Dimension(N+1)],i->[]);

####################################################
Homotopy:=function(k,x)
local htpy;

if k<N-1 or k>N+1 then return R!.homotopy(k,x); fi;
if k=N-1 then 
 if not IsBound(homotopyrec[1][AbsInt(x[1])][x[2]]) then
 homotopyrec[1][AbsInt(x[1])][x[2]]:= HomotopyNminus1([AbsInt(x[1]),x[2]]);
 fi;  
 htpy:=homotopyrec[1][AbsInt(x[1])][x[2]]; 
fi;

if k=N then 
 if not IsBound(homotopyrec[2][AbsInt(x[1])][x[2]]) then
 homotopyrec[2][AbsInt(x[1])][x[2]]:= HomotopyN([AbsInt(x[1]),x[2]]);
 fi;
 htpy:=homotopyrec[2][AbsInt(x[1])][x[2]];
fi;
if k=N+1 then
 if not IsBound(homotopyrec[3][AbsInt(x[1])][x[2]]) then
 homotopyrec[3][AbsInt(x[1])][x[2]]:= HomotopyNplus1([AbsInt(x[1]),x[2]]);
 fi;
 htpy:=homotopyrec[3][AbsInt(x[1])][x[2]];
fi;

if x[1]>0 then return htpy; else return NegateWord(htpy); fi;
end;
####################################################

else
Homotopy:=fail;
fi;

CR:=Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=Homotopy,
                elts:=R!.elts,
                group:=R!.group,
                properties:=
                   [["length",EvaluateProperty(R,"length")],
                    ["reduced",EvaluateProperty(R,"reduced")],
                    ["type","resolution"],
                    ["characteristic",EvaluateProperty(R,"characteristic")]  ]));

return CR;

od;
##############################
##############################
if triple=fail then return R; fi;

end);
#########################################################





#########################################################
InstallGlobalFunction(HAPTietzeReduction_OneLevel,
function(R,i,bound)
local T, s;


s:=R!.dimension(i);
T:=HAPTietzeReduction_OneStep(R,i,bound);

while s>T!.dimension(i) do
s:=T!.dimension(i);
T:=HAPTietzeReduction_OneStep(T,i,bound);
od;

return T;

end);
#########################################################

#########################################################
InstallGlobalFunction(HAPTietzeReduction_Inf,
function(arg)
local R,bound, T, s;
R:=arg[1];
if Length(arg)=2 then bound:=arg[2]; else bound:=infinity; fi;

T:=R;
for s in Reversed([0..Length(R)-1]) do
T:=HAPTietzeReduction_OneLevel(T,s,bound);
od;

return T;
end);
#########################################################

