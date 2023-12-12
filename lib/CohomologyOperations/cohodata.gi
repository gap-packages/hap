
if IsPackageMarkedForLoading("singular","06.07.23") then
#######################################################
#######################################################
InstallGlobalFunction(MinimizeRingRelations,
function(P)
local newP, rels, gens, newrels, newgens, H, newH, i;
rels:=SSortedList(P!.relations);
if Length(rels)=0 then return P; fi;
newrels:=1*rels;
H:=HilbertPoincareSeries(P);

for i in [1..Length(rels)] do
RemoveSet(newrels,rels[i]);
#RemoveSet(newgens,gens[i]);

newP:=rec(
relations:=newrels,
degrees:=P!.degrees,
ring:=P!.ring,
);
newP:=Objectify(GradedAlgebraPresentationType,newP);

newH:=HilbertPoincareSeries(newP);
if not H=newH then
AddSet(newrels,rels[i]);
newP:=rec(
relations:=newrels,
degrees:=P!.degrees,
ring:=P!.ring,
);
newP:=Objectify(GradedAlgebraPresentationType,newP);
fi;
od;

return newP;

end);
#######################################################
#######################################################
else
InstallGlobalFunction(MinimizeRingRelations,
function(x) return x; end);
fi;

#########################################################
#########################################################
InstallGlobalFunction(CohomologicalData,
function(arg)
local G,N,file, alpha, alpha1, A, gens, gensletters, gensletters1, gensdegrees,
d, p, tmp, tmpdir, x, relabel, relabeltwo, pres, rels, r, s,i, k, w,
AppendTo, PrintTo; 

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

G:=arg[1];
if not SSortedList(Factors(Order(G)))=[2] then 
Print("This function is only implemented for 2-groups.\n");
return fail;
fi;
N:=arg[2];
if Length(arg)=2 then
tmpdir := DirectoryTemporary();;
file:=Filename( tmpdir , "cdata.txt" );
else
file:=arg[3];
fi;

if IsBound(hap_cr[IdGroup(G)[1]]) then
if IsBound(hap_cr[IdGroup(G)[1]][IdGroup(G)[2]]) then
if N >= Maximum(  hap_cr[IdGroup(G)[1]][IdGroup(G)[2]]  ) then 
Print("Integer argument is large enough to ensure completeness of cohomology ring presentation.\n\n");
else
Print("Integer argument should be at least ",Maximum(  hap_cr[IdGroup(G)[1]][IdGroup(G)[2]]  )," ","to ensure completeness of cohomology ring presentation.\n\n");
fi;
fi;
fi;
AppendTo(file,"Group order: ", IdGroup(G)[1],"\n");
AppendTo(file,"Group number: ",IdGroup(G)[2],"\n");
AppendTo(file,"Group description: ",StructureDescription(G),"\n\n");

alpha:=
['1','a','b','c','d','e','f','g','h','p','q','r','s','t','u','v','w','x'];
alpha1:=List(alpha,i->[i]);

A:=ModPCohomologyRing(G,N);
gens:=ModPRingGenerators(A);

#gens:=ModPRingGenerators(A);
gensletters:=alpha{[1..Length(gens)]};
gensletters1:=alpha1{[1..Length(gens)]};

AppendTo(file,"Cohomology generators\n");
gensdegrees:=List(gens,A!.degree);
for d in SSortedList(gensdegrees) do
if d>0 then
AppendTo(file, "Degree ",d,": ");
tmp:=Filtered([1..Length(gens)],i->A!.degree(gens[i])=d);
tmp:=gensletters1{tmp};
for x in tmp do
if Position(tmp,x)<Length(tmp) then
AppendTo(file,x,", ");
else
AppendTo(file,x,"\n");
fi;
od;
fi;
od;


AppendTo(file,"\n");


##################################
##################################
relabel:=function(ss) 
local i, s, us;
s:=String(ss);
s:=List(s,i->i);
s:=Filtered(s,i->not i='x');
Add(s,' ');
us:=Filtered([1..Length(s)],i->s[i]='_');
for i in us do
s[i+1]:=gensletters[1+EvalString([s[i+1]])];
if not s[i+2] in ['^','+','*'] then s[i+2]:=' '; fi;
od;
s:=Filtered(s,i->not i='_');
s:=Filtered(s,i->not i=' ');
if Length(s)=0 then return 0; fi;
return s;
end;
##################################
##################################



AppendTo(file,"Cohomology relations\n");
pres:=Mod2CohomologyRingPresentation(A);
pres:=MinimizeRingRelations(pres);
rels:=pres!.relations;
for r in [1..Length(rels)] do
s:=relabel(rels[r]);
AppendTo(file,r,": ",s,"\n");
od;
AppendTo(file,"\n");

AppendTo(file,"Poincare series\n");
p:=HilbertPoincareSeries(pres);
p:=String(p);
p:=List(p,i->i);
for i in [1..Length(p)] do
if p[i]='_' then p[i]:=' '; p[i+1]:=' '; fi;
od;
p:=Filtered(p,i->not i=' ');
AppendTo(file,p,"\n\n");

N:=Maximum(List(gens,A!.degree));

A:=Mod2SteenrodAlgebra(G,2*N);
gens:=ModPRingGenerators(A);
gensletters:=alpha{[1..Length(gens)]};
gensletters1:=alpha1{[1..Length(gens)]};

gensdegrees:=List(gens,A!.degree);



##################################
##################################
relabeltwo:=function(ss)
local i, s, us, t, ii, l;
s:=String(ss);
s:=List(s,i->i);
s:=Filtered(s,i->not i='v');
us:=Filtered([1..Length(s)],i->s[i]='.');
Add(s,' ');
for i in us do
ii:=i+1;
l:=[];
while not s[ii] in ['^','+','*', ' '] do
Add(l,s[ii]);
s[ii]:=' ';
ii:=ii+1;
od;
t:=Basis(A)[EvalString(l)];
s[i+1]:=gensletters[Position(gens,t)];
#if not s[i+2] in ['^','+','*'] then s[i+2]:=' '; fi;
od;
s:=Filtered(s,i->not i='.');
s:=Filtered(s,i->not i=' ');
if Length(s)=0 then return 0; fi;
return s;
end;
##################################
##################################



AppendTo(file,"Steenrod squares\n");
for i in [2..Length(gens)] do
for k in [1..A!.degree(gens[i])-1] do
if k=2^Log(k,2)then
AppendTo(file,"Sq^",k,"(",[gensletters[i]],")=");
w:=Sq(A,k,gens[i]);
if IsZero(w) then AppendTo(file,0,"\n");  
else
w:=PrintAlgebraWordAsPolynomial(A,w,1);
w:=relabeltwo(w);
AppendTo(file,w,"\n");
fi;
fi;
od;
od;

#Exec(Concatenation("display ",file));  
i:=InputTextFile(file);
s:=ReadAll(i);
Print(s);
if Length(arg)=2 then RemoveFile(file); fi;

Print("\n");
return A;
end);
#########################################################
#########################################################

