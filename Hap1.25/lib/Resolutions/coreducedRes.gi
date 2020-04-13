#(C)2009 Graham Ellis

#########################################################
InstallGlobalFunction(TietzeReducedResolution,
function(arg)
local
        R,N,
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
	HmtpyNminus1, HmtpyN, NewHmtpyN, HmtpyNplus1,
	hmtpyrec,tmp,
        newb,bool,i,j,b,e,g,x,w,D,D2;;


##############################################
##############################################
#####
if Length(arg)=1 then R:=StructuralCopy(arg[1]);

####No homotopy present########
#if R!.homotopy=fail then
return HAPTietzeReduction_Inf(R);
#fi;
####No homotopy case done######
#So now this is just a renaming of HAPTietzeReduction_Inf
#and the following lines are irrelevant!


D:=List([0..Length(R)],R!.dimension);

for i in [1..Length(R)-1] do
R:=TietzeReducedResolution(R,i);
od;

D2:=List([0..Length(R)],R!.dimension);

while D2<D do
  D:=D2;
  for i in [1..Length(R)-1] do
    R:=TietzeReducedResolution(R,i);
  od;
  D2:=List([0..Length(R)],R!.dimension);
od;

return R;
fi;
#####
##############################################
##############################################

R:=arg[1];
N:=arg[2];
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
HmtpyNminus1:=[];
for i in [1..R!.dimension(N-1)] do
HmtpyNminus1[i]:=[];
for g in [1..Length(Elts)] do
HmtpyNminus1[i][g]:=StructuralCopy(R!.homotopy(N-1,[i,g]));
od;
od;

hmtpyrec:=List([1..R!.dimension(N+1)],x->[]);;

if Length(R)>N+1 then
HmtpyN:=[];
NewHmtpyN:=[];
for i in [1..R!.dimension(N)] do
HmtpyN[i]:=[];
NewHmtpyN[i]:=[];
od;
fi;

if Length(R)>N+2 then
HmtpyNplus1:=[];
for i in [1..R!.dimension(N+1)] do
HmtpyNplus1[i]:=[];
od;
fi;
########################################################################################



#####################################################################
Action:=function(g,l);
return [l[1],Position(Elts,Elts[g]*Elts[l[2]])];
end;
#####################################################################

#####################################################################
ActionInv:=function(g,l);
return [l[1],Position(Elts,Elts[g]^-1*Elts[l[2]])];
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

if IsInt(pos) then
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
   Add(hmtpyrec[j],[-SignInt(b[i][1])*triple[3],b[i][2],triple[4]]);
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


######################
for j in [1..R!.dimension(N-1)] do
for g in [1..Length(Elts)] do
b:=HmtpyNminus1[j][g];
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
HmtpyNminus1[j][g]:=AlgebraicReduction(newb);

od;
od;
##################################

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

for j in [1..R!.dimension(N-1)] do
for g in [1..Length(Elts)] do
w:=HmtpyNminus1[j][g];
for  i in [1..Length(w)] do
w[i][1]:=SignInt(w[i][1])*Position(modN,AbsInt(w[i][1]));
od;
od;od;

#####################################

PseudoBoundaryN:=PseudoBoundaryN{modN};

#####################################
if Length(R)>N+2 then
for i in [1..Length(modNplus1)] do
for g in [1..Length(Elts)] do
HmtpyNplus1[i][g]:=StructuralCopy(R!.homotopy(N+1, [modNplus1[i],g]));
od;
od;

for i in [1..Length(modNplus1)] do
tmp:=[];
for g in [1..Length(Elts)] do
for x in hmtpyrec[modNplus1[i]] do
Append(HmtpyNplus1[i][g],R!.homotopy(N+1,[x[1],
               Position(Elts,Elts[g]*Elts[x[2]]*Elts[x[3]]^-1)]));
od;
od;
od;
fi;
#####################################

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


for j in modN do
for g in [1..Length(Elts)] do
b:=R!.homotopy(N,[j,g]);
newb:=[];
for x in b do
if  AbsInt(x[1]) in modNplus1 then
Add(newb,[SignInt(x[1])*Position(modNplus1,AbsInt(x[1])),x[2]]);
fi;
od;
HmtpyN[j][g]:=newb;
od;
od;

for j in [1..Length(modN)] do
for g in [1..Length(Elts)] do
NewHmtpyN[j][g]:=HmtpyN[modN[j]][g];
od;
od;
HmtpyN:=NewHmtpyN;


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



####################################################
Homotopy:=function(n,x);
if not n in [N-1,N,N+1] then return R!.homotopy(n,x); fi;
if x[1]>0 then 
if n=N-1 then return HmtpyNminus1[x[1]][x[2]]; fi;
if n=N then return HmtpyN[x[1]][x[2]]; fi;
if n=N+1 then return HmtpyNplus1[x[1]][x[2]]; fi;
fi;
if x[1]<0 then
if n=N-1 then return NegateWord(HmtpyNminus1[AbsInt(x[1])][x[2]]); fi;
if n=N then return NegateWord(HmtpyN[AbsInt(x[1])][x[2]]); fi;
if n=N+1 then return NegateWord(HmtpyNplus1[AbsInt(x[1])][x[2]]); fi;
fi;

end;
####################################################


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
end);
#########################################################
