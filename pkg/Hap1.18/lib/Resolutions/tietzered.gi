#(C)2009 Graham Ellis

#########################################################
InstallGlobalFunction(HAPTietzeReduction_OneStep,
function(R,N,bound)
local
	CR,
        Dimension,
        Boundary,
        Homotopy,
        PseudoBoundary,
        PseudoBoundaryN,
        PseudoBoundaryNplus2,
	FindFreeFace,
	Action, ActionInv, Elts,
 	modN, modNplus1,
	triple,
	tmp,
        newb,bool,i,j,b,e,g,x,w,D,D2;;

Elts:=R!.elts;
modN:=[];
modNplus1:=[];

PseudoBoundary:=List([1..R!.dimension(N+1)],i->
                              StructuralCopy(R!.boundary(N+1,i)));
PseudoBoundaryN:=List([1..R!.dimension(N)],i->
                              StructuralCopy(R!.boundary(N,i)));
if Length(R)>N+1 then
PseudoBoundaryNplus2:=List([1..R!.dimension(N+2)],i->
                              StructuralCopy(R!.boundary(N+2,i)));
fi;

###################################################################################


#####################################################################
Action:=function(g,l)
local pos, h;
h:=Elts[g]*Elts[l[2]];
pos:=Position(Elts,h);
if pos=fail then
Add(Elts, h);
pos:=Length(Elts);
fi;
return [l[1],pos];
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
return [l[1],pos];

end;
#####################################################################

########################################################
FindFreeFace:=function()
local i,b,pos,y, wrd,e, g;
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
b:=StructuralCopy(PseudoBoundary[i]);
pos:=PositionProperty(b,x->AbsInt(x[1]) =y); 
wrd:=b{Concatenation([1..pos-1],[pos+1..Length(b)])};
g:=b[pos][2];
if b[pos][1]>0 then
wrd:=NegateWord(wrd);
e:=b[pos][1];
else e:=-b[pos][1];
fi;
wrd:=List(wrd,x->ActionInv(g,x));
return [e,wrd,i,g];
fi;
od;
return fail;
end;
########################################################


####################
####################
while true do
triple:=FindFreeFace();
if triple=fail or Length(modN)>0 
then break; fi;
e:=triple[1];

Add(modNplus1,triple[3]);
Add(modN,e);

###############
for j in [1..Length(PseudoBoundary)] do
b:=PseudoBoundary[j];
newb:=[];
for i in [1..Length(b)] do
   if not AbsInt(b[i][1])=e then
   Add(newb,StructuralCopy(b[i]));
   else
   w:=StructuralCopy(triple[2]);
   if b[i][1]<0 then
   w:=NegateWord(w);
   fi;
   w:=List(w,x->Action(b[i][2],x));
   Append(newb,w);
   fi;
od;
PseudoBoundary[j]:=AlgebraicReduction(newb);
od;
######################



od;
##############################
##############################

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
if not i in [N,N+1] then return R!.dimension(i); fi;
if i=N then return  Length(PseudoBoundaryN); fi;
return Length(PseudoBoundary); 
end;
####################################################



####################################################
Boundary:=function(n,i);
if not n in [N,N+1,N+2]  then return R!.boundary(n,i); fi;
if i>0 then
if n=N then return PseudoBoundaryN[i]; fi;
if n=N+1 then return PseudoBoundary[i]; fi;
if n=N+2 then return PseudoBoundaryNplus2[i]; fi;
fi;
if i<0 then
if n=N then return NegateWord(PseudoBoundaryN[AbsInt(i)]); fi;
if n=N+1 then return NegateWord(PseudoBoundary[AbsInt(i)]); fi;
if n=N+2 then return NegateWord(PseudoBoundaryNplus2[AbsInt(i)]); fi;
fi;

end;
####################################################

CR:=Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=fail,
                elts:=R!.elts,
                group:=R!.group,
                properties:=
                   [["length",EvaluateProperty(R,"length")],
                    ["reduced",EvaluateProperty(R,"reduced")],
                    ["type","resolution"],
                    ["characteristic",EvaluateProperty(R,"characteristic")]  ]));

return CR;
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

