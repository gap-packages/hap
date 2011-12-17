#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(PolytopalGenerators,
function(GG,v)
local
	Points, Vertices, 
	G, D, CG, N, 
	proj, 
	StringToVector,
	input,
	EdgeGenerators,
	tmp, x, w, i, y, p;

if not (IsPermGroup(GG) or IsMatrixGroup(GG)) then 
Print("The group G must be a permutation or matrix group.\n");
return fail;
fi;

G:=Filtered(Elements(GG),x->not x=Identity(GG));
D:=Length(v);
CG:=List([1..D], i->0);   	#CG will eventually be the centre of gravity
Points:=[];			#of the polytope.
Vertices:=[];
EdgeGenerators:=[];

###################### CALCULATE CENTRE OF GRAVITY ##################
if IsPermGroup(GG) then

for x in G do
w:=[];
	for i in [1..D] do
	Append(w,[v[i^x]]);
	od;
Append(Points, [w]); 
      	for i in [1..D] do
        CG[i]:=CG[i]+w[i];
      	od;
od;

for i in [1..D] do
CG[i]:=CG[i]/Order(GG);
od;

else

for x in G do
	w:=x*v;
	Append(Points, [w]);
        for i in [1..D] do
        CG[i]:=CG[i]+w[i];
        od;
od;

for i in [1..D] do
CG[i]:=CG[i]/Order(GG);
od;

fi;
##################### CENTRE OF GRAVITY DONE ########################


	################# PROJECTION ################################
	N:=CG-v;
	proj:=function(w)
	local k, m  ;
	m:=(w-v)*N;
	k:= ((CG-v)*N)/m;
	return v+(k*(w-v));
	end;
	#############################################################
	
################## CALCULATE HULL OF PROJECTED POINTS ###############
for i in [1..Length(Points)] do   
Points[i]:=proj(Points[i]);	  
od;

AppendTo("tmpIn.log","POINTS","\n");
for p in Points do
AppendTo("tmpIn.log",1);
	for i in [1..Length(p)] do
	AppendTo("tmpIn.log"," ",p[i]);
	od;
AppendTo("tmpIn.log","\n");
od;

Exec("beneath_beyond tmpIn.log POINTS");
################# HULL CALCULATED ###################################

#####################################################################
StringToVector:=function(s)
local V,x,y;

V:=[]; y:=[];

for x in s do
if not x=' ' then Append(y,[x]);
else
Append(V,[Rat(Chomp(y))]); y:=[]; 
fi;
od;
Append(V,[Rat(Chomp(y))]);


V[1]:='G';
V:=Filtered(V,x->not x='G');
return V;
end;
#####################################################################

################# READ VERTICES #####################################
input:=InputTextFile("tmpIn.log");
tmp:="hello";
while not tmp="VERTICES\n" do
tmp:=ReadLine(input);
od;

tmp:=ReadLine(input);
i:=1;
while tmp[1]='1' do
Vertices[i]:=tmp;
tmp:=ReadLine(input);
i:=i+1;
od;

Exec("rm tmpIn.log");
Apply(Vertices, s->StringToVector(s));
################ VERTICES READ ######################################

################ RECOVER THE EDGE GENERATORS ########################
for w in Points do
if proj(w) in Vertices then
Append(EdgeGenerators,[G[Position(Points,w)]]);
fi;
od;
################ EDGE GENERATORS RECOVERED ##########################
return EdgeGenerators;

end);
#####################################################################
