#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(IdentityAmongRelators,
function(arg)
local  
		R,idnum,
		Dimension,
		Boundary,
		Elts, Mult, Inv,
		Frels, rels,
		Fgens,gens,
		FirstBoundaryHomomorphism,
		Boundary2Relator, 
		start,
		ActWord,
		idnt,
		gennum,
		CommonChord,
		Amalgamate,
		b, r, x,i,X;

R:=arg[1];
idnum:=arg[2];

if not (IsHapResolution(R) or IsHapNonFreeResolution(R)) then
Print("This function must be applied to a resolution. \n");
return fail;
fi;

if not EvaluateProperty(R,"reduced")=true then
if R!.dimension(0)>1 then
Print("This function must be applied to a REDUCED resolution. \n");
return fail; fi;
fi;

if not EvaluateProperty(R,"characteristic")=0 then
Print("This function only works in characteristic 0. \n");
return fail;
fi;


Dimension:=R!.dimension;
Boundary:=R!.boundary;
Elts:=R!.elts;
Frels:=[];
start:=List([1..Dimension(2)],x->List(Boundary(2,x),y->y[2]));
start:=SortedList(Intersection(start))[1];
gens:=[];


#####################################################################
Mult:=function(g,h);
return Position(Elts,Elts[g]*Elts[h]);
end;
#####################################################################

#####################################################################
Inv:=function(g);
return Position(Elts,Elts[g]^-1);
end;
#####################################################################

#####################################################################
FirstBoundaryHomomorphism:=function(x)
local r;
r:=Boundary(1,x[1]);
r:=List(r,y->Mult(x[2],y[2]));
if x[1]>0 then return r;
else return Reversed(r); fi;
end;
#####################################################################

#####################################################################
Boundary2Relator:=function(b)
local c, rel, w;

b:=SortedList(AlgebraicReduction(b));
rel:=[start];

while Length(b)>0 do
	for x in b do
	w:=FirstBoundaryHomomorphism(x);
	if w[1]= rel[Length(rel)] then
	Append(rel, [w[2]]); RemoveSet(b,x); break; 
	else
	   if w[2]= rel[Length(rel)] then
	   Append(rel, [w[1]]); RemoveSet(b,x); break;
	   fi;
	fi;
	od;
od;

return rel;
end;
#####################################################################

for r in [1..Dimension(2)] do
Append(Frels,[Boundary2Relator(Boundary(2,r))]);
od;

for r in Frels do
if (not Inv(r[2]) in gens) then AddSet(gens,r[2]);fi;
if (not Inv(r[Length(r)-1]) in gens) then AddSet(gens,r[Length(r)-1]);fi;
od;


#####################################################################
gennum:=function(r)
local g,h;

        for g in gens do
        if Mult(r[1],g)=r[2] then h:=Position(gens,g); break; fi;
        if Mult(r[1],Inv(g))=r[2] then h:=Position(gens,g); break; fi;
        od;

h:=(h-1) mod 6;
h:=h+1;
return h;
end;
#####################################################################

#####################################################################
CommonChord:=function(A,BB)
		#Returns the common contractible chord between loops
		#A and B if such a chord exists. Otherwise it returns
		#fail. This is a clumsy piece of code! 

local B,	C, i, chord,chordA,chordB,comp ;

B:=StructuralCopy(BB);
C:=Intersection(A,B);
if Length(C)<2 then return fail; fi;

chordA:=[];
chordB:=[];

for i in [1..Length(B)-1] do
if B[i] in C and B[i+1] in C then
Add(chordB, [B[i],B[i+1]]);
fi;
od;

for i in [1..Length(A)-1] do
if A[i] in C and A[i+1] in C then
Add(chordA, [A[i],A[i+1]]);
fi;
od;

chord:=Intersection(chordA,chordB);

if not Length(chord)=Length(C)-1 then
return fail;
fi;

#########################
if not B[1] in C  then
chord:=List(chordB,x->x[1]);
Add(chord,chordB[Length(chordB)][2]     );

chordA:=Reversed(B{[1..Position(B,chord[1])]});
chordB:=Reversed(B{[Position(B,chord[Length(chord)])..Length(B)-1]});
comp:=Concatenation(chordA,chordB);
fi;
########################

########################
if  B[1] in C then 
chordA:=[];
chordB:=[];
for x in B{[2..Length(B)]} do
if x in C then Add(chordA,x); else break;fi;
od;
for x in Reversed(B) do
if x in C then Add(chordB,x); else break;fi;
od;
chord:=Concatenation(Reversed(chordB),chordA);

if B[2] in C and B[Length(B)-1] in C and B[2] in C then 
comp:=Reversed(B{[Position(B,chord[Length(chord)])..Position(B,chord[1])]});
fi;

if B[2] in C and not B[Length(B)-1] in C then
comp:=Reversed(B{[Position(B,chord[Length(chord)])..Length(B)]});
fi;

if not B[2] in C and  B[Length(B)-1] in C then
comp:=Reversed(B{[1..Position(B,chord[1])]});
fi;


fi;
#######################



return [chord,comp];

end;
#####################################################################

#####################################################################
Amalgamate:=function(A,B)
local 	C,D,		#Here A is the big loop and B the small loop.
	Begin,End,amalg;

D:=CommonChord(A,B);
C:=D[2];;
if C=fail then return fail; fi;

if not A[1] in D[1] then
Begin:=A{[1..Position(A,C[1])-1]};
End:=A{[Position(A,C[Length(C)])+1..Length(A)]};
amalg:=Concatenation(Begin,C,End);
fi;

if A[1] in D[1] then
End:=A{[Position(A,C[Length(C)])+1..Position(A,C[1])-1]};

amalg:=Concatenation(C,End);
Add(amalg,C[1]);
fi;


return amalg;
end;
#####################################################################
#####################################################################

idnt:=[];

for b in Boundary(3,idnum) do
x:=Boundary(2,AbsInt(b[1]));
x:=Boundary2Relator(x);
if SignInt(b[1])=-1 then x:=Reversed(x); fi;
x:=List(x,t->Mult(b[2],t));
Add(idnt,x);
od;
####################################################################


#return CommonChord;
return [Amalgamate,idnt];
end);
#####################################################################


