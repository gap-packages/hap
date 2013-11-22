#(C) Graham Ellis, 2..5-2006

#####################################################################
#####################################################################
InstallGlobalFunction(ResolutionAsphericalPresentation,
function(arg)
local
	F,rels,N,Dimension,
	Boundary, PseudoBoundary,
	Homotopy, HomotopyRecord,
	G,
	EltsG,
	gensF, gensG,
	FhomG,
	toggle,
	i,R,B,w,ww,x;

if Length(arg)=3 then
F:=arg[1]; rels:=arg[2]; N:=arg[3]; 
         if Length(rels)>0 then G:=F/rels;
         else G:=F; fi;
fi;
if Length(arg)=2 then
G:=arg[1]; 
F:=FreeGroupOfFpGroup(G);
rels:=RelatorsOfFpGroup(G); N:=arg[2];
fi;

toggle:=false;
gensF:=GeneratorsOfGroup(F);
#if Length(rels)>0 then G:=F/rels;
#else G:=F; fi;
gensG:=GeneratorsOfGroup(G);
FhomG:=GroupHomomorphismByImagesNC(F,G,gensF,gensG);
EltsG:=[Identity(G)];

if IsFreeGroup(G) and Length(gensG)=1 then toggle:=true;
x:=HAPconstant; ##CHANGE x IF THIS CAUSES PROBLEMS!
for i in [1..x]  do
Append(EltsG,[F.1^i,F.1^-i]);
od;
fi;

#####################################################################
Dimension:=function(n);

if n=0 then return 1; fi;
if n=1 then return Length(gensF); fi;
if n=2 then return Length(rels); fi;
if n<0 or n>2 then return 0; fi;

end;
#####################################################################

PseudoBoundary:=[[],[]];

for i in [1..Dimension(1)] do
Append(PseudoBoundary[1], [   [ [-1,1], [1,i+1] ]   ]); 
Append(EltsG,[gensG[i]]);
od;

for i in [1..Dimension(2)] do
R:=LetterRepAssocWord(rels[i]);
B:=[];
w:=Identity(F);

for x in R do
Append(EltsG,[Image(FhomG,w)]);
#Append(B,[[x,Length(EltsG)]]);
w:=w*gensF[AbsoluteValue(x)]^SignInt(x);
if SignInt(x)=1 then
Append(B,[[x,Length(EltsG)]]);
else
Append(B,[[x,Length(EltsG)+1]]);
fi;

od;
Append(EltsG,[Image(FhomG,w)]);

Append(PseudoBoundary[2],[B]);


od;

#####################################################################
Boundary:=function(n,i);
if n<1 or n>2 then return [ ]; fi;

if i>0 then 
return PseudoBoundary[n][i];
else return NegateWord(PseudoBoundary[n][-i]);
fi;
end;
#####################################################################

if not toggle then Homotopy:=fail;
else

#####################################################################
HomotopyRecord:=function(i,x)
local g,k,j,hty;
if i>0 then return []; fi;
g:=EltsG[x[2]];
if g=Identity(F) then return []; fi;
for k in [-100..100] do			#THIS SLOPPINESS MIGHT CAUSE PROBLEMS
if F.1^k=g then j:=k; break; fi;
od;

if j>0 then k:=Position(EltsG,F.1^(j-1)); 
hty:=Concatenation(ShallowCopy(HomotopyRecord(i,[x[1],k])),
				[[1,k]]); 

else k:=Position(EltsG,F.1^(j+1)); 
hty:=Concatenation(ShallowCopy(HomotopyRecord(i,[x[1],k])),
                                [[-1,x[2]]]); 
fi;  #Need to think more about this!


return hty;
end;
#####################################################################

#####################################################################
Homotopy:=function(i,x);
if SignInt(x[1])=1 then return HomotopyRecord(i,x);
else
return NegateWord(HomotopyRecord(i,x));
fi;
end;
#####################################################################
fi;

return Objectify(HapResolution,
	    rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=Homotopy,
            elts:=EltsG,
            group:=G,
            properties:=
            [["length",N],
	     ["characteristic",0],
	     ["type","resolution"],
	     ["reduced",true]]  ));
end);
#####################################################################
#####################################################################
