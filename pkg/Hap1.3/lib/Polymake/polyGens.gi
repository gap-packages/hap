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
	Faces, FacesFinal, FacesFn,
	Index, IndexFn,
	Count,
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

Apply(Vertices, s->StringToVector(s));
################ VERTICES READ ######################################

################ RECOVER THE EDGE GENERATORS ########################
for w in Points do
p:=proj(w);
if p in Vertices then
x:=G[Position(Points,w)];
EdgeGenerators[Position(Vertices,p)]:=x;
fi;
od;
################ EDGE GENERATORS RECOVERED ##########################

################ READ HASSE DIAGRAM #################################
Exec("polymake tmpIn.log HASSE_DIAGRAM >tmp2In.log");
Exec("rm tmpIn.log");

Faces:=[];
Index:=[];

	#############################################################
	IndexFn:=function(tmp)
	local V,x,y;
	V:=[];
	y:=[];
	
	for x in tmp do
	if (not x=' ') and (not x='\n') then Append(y,[x]);
	else
	Append(V,[Rat(Chomp(y))]); y:=[]; fi;
	od;

	return V;
	end;
	############################################################

	############################################################
	FacesFn:=function(tmp)
	local V,x,y;
	V:=[]; y:=[];
	
	for x in tmp do
	if (not x='<') and (not x='(') and (not x='{') then
	if (not x=' ') and (not x='}') then Append(y,[x]);
	else
	if Length(y)>0 then
        Append(V,[Rat(Chomp(y))]); y:=[]; fi; fi;
	fi;
	if x='}' then break; fi;
	od;
	
	if Length(V)>0 then
	return [V];
	else
	return []; fi;
	end;
	############################################################

input:=InputTextFile("tmp2In.log");
tmp:=ReadLine(input);
tmp:=ReadLine(input);
Index:= IndexFn(tmp);

tmp:=ReadLine(input);

while Length(tmp)>2 do
Append(Faces, FacesFn(tmp));
tmp:=ReadLine(input);
od;
Exec("rm tmp2In.log");

if Length(Faces[1])=1 then

	FacesFinal:=[];
	for i in [1..Length(Index)-1] do
	Append(FacesFinal,[[Index[i]..Index[i+1]-1]]);
	od;
	Append(FacesFinal,[[Index[i+1]]]);
	FacesFinal:=(List(FacesFinal,x->List(x,i->Faces[i])));

else

	FacesFinal:=[[1..Index[1]]];
	for i in [2..Length(Index)-1] do
	Append(FacesFinal,[[Index[i]+1..Index[i+1]]]);
	od;
	FacesFinal:=Reversed(List(FacesFinal,x->List(x,i->Faces[i])));
	
fi;


############### HASSE DIAGRAM READ ##################################
return rec(
             generators:=EdgeGenerators,
             hasseDiagram:=FacesFinal,
	     vector:=v);
	     

end);
#####################################################################
