#(C) Graham Ellis, 2005-2006


#####################################################################
InstallGlobalFunction(ResolutionFiniteGroup,
function(arg)

local
	R,
	Gens,
	K,
	tietze,
	G,
	Elts, ExtendedElts, 
	N,
	MT,
	Action,
	ChangeSign,
	Abs,
	MaxComplex,
	Dimension,
	Boundary,
	PseudoBoundary,
	ContractionMatrix,
	Contraction,
	Homotopy,
	ComputedContractions,
	Differential,
	CellValue,
	FirstZero,
	IsFinished,
	FindConsequences,
	NextResTerm,
	Spheres,
	DiffLengths,
	InitComputedContractions,
	saveSpace,
	Charact,
	AlgebraicRed,
	ExtendRes, Extendible,
	i, ii, iso,
#################################
AbsInt,				#
SignInt;			#
				#
AbsInt:=AbsInt_HAP;		#
SignInt:=SignInt_HAP;		#
#################################

if IsGroup(arg[1]) then Gens:=GeneratorsOfGroup(arg[1]);
if Length(Gens)=0 then Gens:=[Identity(arg[1])]; fi;
else Gens:=arg[1]; fi;
Gens:=SSortedList(StructuralCopy(Gens));
K:=StructuralCopy(arg[2]);
if Length(arg)>2 then tietze:=arg[3]; else tietze:=false; fi;

if Length(arg)>3 then 
		if IsInt(arg[4]) then Charact:=arg[4]; fi;
else Charact:=0; fi;

if Length(arg)>4 then 
 if arg[5]="extendible" then Extendible:=true; saveSpace:=false;
 else saveSpace:=arg[5]; Extendible:=false; fi;
else saveSpace:=false; Extendible:=false;
fi;

G:=Group(Gens);
N:=Order(G);
#if IsCyclic(G) then return ResolutionFiniteCyclicGroup(G,K); fi; #Added 15 August 2019

######################################################
if IsMatrixGroup(G) then

iso:=IsomorphismPermGroup(G);
#R:=ResolutionFiniteGroup(Image(iso),N); #CHANGED 26/11/2018
R:=ResolutionFiniteGroup(Image(iso,G),K);
R!.elts:=List(R!.elts,x->PreImagesRepresentative(iso,x));
R!.group:=G;
return R;

fi;
######################################################

Elts:=Elements(G);

#############
	##### If this piece of code is used then Elts will not be
	##### strictly sorted and things may get slow!! I should fix this.
	#####
if not Elts[1]=Identity(G) then
Elts:=[];
for i in Elements(G) do
Add(Elts,i);
od;

i:=Position(Elts,Identity(G));
Elts[i]:=Elts[1];
Elts[1]:=Identity(G);
fi;

#if IsCyclic(G) then
#if not Order(Elts[2])=Order(G) then
#Elts:=List(Elts,x->x);
#i:=PositionProperty(Elts,x->Order(x)=Order(G));
#ii:=Elts[i];
#Elts[i]:=Elts[2];
#Elts[2]:=ii;
#fi;
#fi;
	#####
	#####
#############


RemoveSet(Gens,Identity(G));

ExtendedElts:=List(Gens,g->Position(Elts,g));	
Append(ExtendedElts,[1..N]);			
#Append(ExtendedElts,Reversed([1..N])); #This line added 10/12/2009

if Charact=0 then AlgebraicRed:=AlgebraicReduction;
else
	AlgebraicRed:=function(w);
	return AlgebraicReduction(w,Charact);
	end;
fi;

if Order(G)<5096  then
MT:=MultiplicationTable(Elts);

#####################################################################
Action:=function(g,l);
return [l[1],MT[g][l[2]]];
end;
#####################################################################

else

#####################################################################
Action:=function(g,l);
return [l[1],Position(Elts,Elts[g]*Elts[l[2]])];
end;
#####################################################################

fi;


#####################################################################
ChangeSign:=function(j,b);
if j>0 then return b; else 
return List(b,x->[-x[1],x[2]]); fi;
end;
#####################################################################

#####################################################################
Abs:=function(l)
local r;

r:=ShallowCopy(l);
Apply(r,x->[AbsInt(x[1]),x[2]]);

return r;
end;
#####################################################################

#####################################################################
MaxComplex:=[];
MaxComplex[1]:=[[1..N]];
Apply(MaxComplex[1][1],i->0);
MaxComplex[1][1][1]:=1;
PseudoBoundary:=[];
ContractionMatrix:=[];
ComputedContractions:=[];
#####################################################################

#####################################################################
Dimension:=function(i);
if i<0 then return 0; fi;
if i=0 then return 1; fi;
return Length(PseudoBoundary[i]); 
end;
#####################################################################

#####################################################################
Boundary:=function(i,j);
if i<=0 then return []; else 
return ChangeSign(j,PseudoBoundary[i][AbsInt(j)]); fi;
end;
#####################################################################

#####################################################################
CellValue:=function(i,MC,p)		#MC=MaxComplex[i]
local x,l,v,e,q;

l:=ShallowCopy(Boundary(i,p[1]));
Apply(l,x->Action(p[2],x));
v:=Length(l);				#l is the boundary of cell p
					#where p has dimension i.
for e in l do
v:=v-MC[AbsInt(e[1])][e[2]];
od;

q:=0;
if v = 1 then 
for x in l do
if MC[AbsInt(x[1])][x[2]] = 0 then q:=x; break; fi;
od;
fi;

if tietze then
   if (v = 0)  then Add(Spheres,l);fi;	
fi;					#Spheres is a list of contractible
					#spheres. 
return [v,q];
end;
#####################################################################

#####################################################################
FirstZero:=function(MC,i)
local j,g, temp, temp2, compare, p,q, LengthOfDiff;

temp:=[];

for j in [1..Length(MC)] do
for g in [1..N] do
#for g in ExtendedElts do    #TRIED 11 May 2019, BUT WAS SLOWER!!
if MC[j][g]=0 then Add(temp,[j,g]);
fi;
od;
od;

#return Random(temp);

LengthOfDiff:=function(p);
if DiffLengths[p[1]][p[2]]=0 then 
  DiffLengths[p[1]][p[2]]:=Length(Differential(i,p)); fi;###################
##################################################Modified 5 Nov 2015#######
#DiffLengths[p[1]][p[2]]:=Reversed([Length(Differential(i,p)),        ######
#Length(SSortedList(List(Differential(i,p),a->AbsInt(a[1]))))]); fi;  ######
return DiffLengths[p[1]][p[2]];
end;

temp2:=List(temp, LengthOfDiff );

return temp[Position(temp2,Minimum(temp2))];

end;
#####################################################################

#####################################################################
IsFinished:=function(MC);
if Product(Flat(MC)) = 1 then return true; 
else return false; fi;
end;
#####################################################################

#####################################################################
Contraction:=function(i,x)		#x is an (i-1)-cell and the output
local l,m,b,c,y,z;			#is a collection of i-cells. For technical
					#reasons we need to allow i=0. The sign factors
if i<1 then return [[-1,1]]; else	#in the function were arrived at by trial and
					#error, rather than mathematical reasoning!
					#When Contraction computes a value for (i,x)
					#it stores it in ComputedContractions[i] and uses
					#this computed value subsequently.
					
   if ComputedContractions[i][AbsInt(x[1])][x[2]]=0 then
   z:=[AbsInt(x[1]),x[2]];
      if ContractionMatrix[i][z[1]][z[2]]=1 then return []; 
      else
      m:=ContractionMatrix[i][z[1]][z[2]];
      b:=ShallowCopy(Boundary(i,m[1]));
      Apply(b,y-> Action(m[2],y));
      b:=ChangeSign(-x[1],b);

      c:=ChangeSign(-x[1],[m]);
      for y in b do
         if (not Abs([y])=Abs([x])) 
         then Append(c,Contraction(i,y));
         fi;
      od;

      if i<K or Extendible then ComputedContractions[i][z[1]][z[2]]:=
                  ChangeSign(x[1],c); fi; 	
      return c;
      fi;
   else

   c:=ComputedContractions[i][AbsInt(x[1])][x[2]];
   return ChangeSign(x[1],c);
   fi;

fi;
end;
#####################################################################

#####################################################################
Homotopy:=function(i,p); 	#It is useful to have a second name
if i <0 then return fail; fi;	#for Contraction with the correct indexing!
return ChangeSign(-1,Contraction(i+1,p));	
end;
#####################################################################

#####################################################################
Differential:=function(i,p)  	#The boundary of an i-cell corresponding
local j,k,l,x,Diff;		#to an (i-1)-cell p in MaxComplex. The cell
				#p itself should be in the boundary.
j:=p[1];
k:=p[2];
Diff:=[p];
if i=1 then l:=[1]; else
l:=ShallowCopy(PseudoBoundary[i-1][AbsInt(j)]); 
Apply(l,x->Action(k,x));	#l is the boundary of p.
l:=ChangeSign(j,l);
fi;

for x in l do
Append(Diff,Contraction(i-1,x)); 
od;

return Diff;
end;
#####################################################################

#####################################################################
FindConsequences:=function(i,MC)  	#MC=StructuralCopy(MacComplex[i])
local j,g,c,p,toggle, SignIntLoc,CellValueLoc,iterset;


SignIntLoc:=SignInt;
CellValueLoc:=CellValue;

toggle:=true;

if i<K or Extendible then
iterset:= Concatenation([1..Dimension(i)],Reversed([1..Dimension(i)]));
while toggle do
toggle:=false;
for g in ExtendedElts do
for j in iterset do
#for j in [1..Dimension(i)]  do #Changed this line 10/12/2009
c:=CellValueLoc(i,MC,[j,g]); p:=c[2];
if c[1]=1 then MC[AbsInt(p[1])][p[2]]:=1;
ContractionMatrix[i][AbsInt(p[1])][p[2]]:=[SignInt(p[1])*j,g];
MaxComplex[i+1][j][g]:=1;
toggle:=true;fi;
od;
od;
od;

else

while toggle do
toggle:=false;
for g in ExtendedElts do
for j in [1..Dimension(i)] do
c:=CellValueLoc(i,MC,[j,g]); p:=c[2];
if c[1]=1 then MC[AbsInt(p[1])][p[2]]:=1;
MaxComplex[i+1][j][g]:=1;
toggle:=true;fi;
od;
od;
od;

fi;
end;
#####################################################################

#####################################################################
InitComputedContractions:=function(i)
local ii,j;

ComputedContractions[i]:=[];
for ii in [1..Dimension(i-1)] do
ComputedContractions[i][ii]:=[];
for j in [1..N] do
ComputedContractions[i][ii][j]:=0;
od;
od;

end;
#####################################################################

#####################################################################
NextResTerm:=function(i)
local ii, p, MC,l,j,Diff;
PseudoBoundary[i]:=[];
MaxComplex[i+1]:=[];
ContractionMatrix[i]:=ShallowCopy(MaxComplex[i]);
MC:=StructuralCopy(MaxComplex[i]);
Spheres:=[];
DiffLengths:=[];

for ii in [1..Length(MC)] do
DiffLengths[ii]:=[];
for j in [1..N] do
DiffLengths[ii][j]:=0;
od;
od;

InitComputedContractions(i); 

l:=[1..N];Apply(l,x->0);l[1]:=1;

if i<K or Extendible then
while not IsFinished(MC) do
p:=FirstZero(MC,i);
if tietze then 
   Diff:=TietzeReduction(Spheres,AlgebraicRed((Differential(i,p))));
else
   Diff:=AlgebraicRed(Differential(i,p));
fi;
Add(PseudoBoundary[i],Diff);
Add(MaxComplex[i+1],ShallowCopy(l));
MC[p[1]][p[2]]:=1;
ContractionMatrix[i][p[1]][p[2]]:=[Length(MaxComplex[i+1]),1];
FindConsequences(i,MC);
od;

else

while not IsFinished(MC) do
p:=FirstZero(MC,i);
if tietze then
Diff:=TietzeReduction(Spheres,AlgebraicRed((Differential(i,p))));
else
Diff:=AlgebraicRed(Differential(i,p));
fi;


Add(PseudoBoundary[i],Diff);
Add(MaxComplex[i+1],ShallowCopy(l));
MC[p[1]][p[2]]:=1;
FindConsequences(i,MC);
od;

fi;

DiffLengths:=0; MC:=0;
end;
#####################################################################

R:=Objectify(HapResolution,
		rec(
		dimension:=Dimension,
		boundary:=Boundary,
		homotopy:=Homotopy,
		elts:=Elts,
		group:=G,
                vectorField:=ContractionMatrix,
		properties:=
		   [["length",0],
		    ["reduced",true],
		    ["type","resolution"],
		    ["characteristic",Charact]  ])); 

#####################################################################
ExtendRes:=function()
local i;
i:=R!.properties[1][2]+1;
NextResTerm(i);
R!.properties[1][2]:=i;
MaxComplex[i]:=[];
if i>1 and saveSpace then InitComputedContractions(i-1); fi;
end;
#####################################################################
if Extendible then
R!.extend:=ExtendRes;
fi;

for i in [1..K] do
ExtendRes();
od;


return R;
end);
#####################################################################
