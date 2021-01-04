#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(PolytopalGenerators,
function(GG,v)
local
	Points, Vertices, 
        Action,
	G, D, CG, N, 
	proj, 
	StringToVector,
	input,
	EdgeGenerators,
	Faces, FacesFinal, FacesFn,
	Index, IndexFn, tmplst,
	Count,
	tmp, x, w, i, y, p,
        AppendTo, PrintTo,
	tmpInlog, tmp2Inlog,tmp3Inlog, tmp4Inlog, tmpDir;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();
tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmp2Inlog:=Filename(tmpDir,"tmp2In.log");
tmp3Inlog:=Filename(tmpDir,"tmp3In.log");
tmp4Inlog:=Filename(tmpDir,"tmp4In.log");

if not (IsPermGroup(GG) or IsMatrixGroup(GG)) then 
Print("The group G must be a permutation or matrix group.\n");
return fail;
fi;

if IsPermGroup(GG) then
GG:=Image(PermToMatrixGroup(GG));
fi;

#####################################################################
Action:=function(g,V) ;
return g*V;;   #This actually works!
end;
#####################################################################

G:=[];
tmplst:=[];

for x in GG do
w:=Action(x,v);
if not w in tmplst and not w=v then
Add(G,x);
Add(tmplst,w);
fi;
od;


D:=Length(v);
CG:=List([1..D], i->0);   	#CG will eventually be the centre of gravity
Points:=[];			#of the polytope.
Vertices:=[];
EdgeGenerators:=[];

###################### CALCULATE CENTRE OF GRAVITY ##################
for x in G do
Add(Points, Action(x,v)); 
od;

Points:=SSortedList(Points);

G:=[];
for w in Points do
for x in GG do
if Action(x,v)=w then Add(G,x); break; fi;
##if Action(x,v) in Points then Add(G,x);  fi;

od;
od;

for w in Points do
    CG:=CG+w;
od;

CG:=CG/Size(Points);
##################### CENTRE OF GRAVITY DONE ########################


	################# PROJECTION ################################
	N:=CG-v;    #This might be parallel to an edge!!! 
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

AppendTo(tmpInlog,"POINTS","\n");
for p in Points do
AppendTo(tmpInlog,1);
	for i in [1..Length(p)] do
	AppendTo(tmpInlog," ",p[i]);
	od;
AppendTo(tmpInlog,"\n");
od;

#Exec(Concatenation(POLYMAKE_PATH,tmpInlog,"  beneath_beyond VERTICES > ",tmp3Inlog));   #Added "beneath_beyond"
Exec(Concatenation(POLYMAKE_PATH,tmpInlog," VERTICES > ",tmp3Inlog));
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

return V{[2..Length(V)]};
end;
#####################################################################

################# READ VERTICES #####################################
input:=InputTextFile(tmp3Inlog);

tmp:=ReadLine(input);
tmp:=ReadLine(input);
i:=1;
while tmp[1]='1' do
Vertices[i]:=tmp;
tmp:=ReadLine(input);
i:=i+1;
od;

Apply(Vertices, s->StringToVector(s));
RemoveFile(tmp3Inlog);

################ VERTICES READ ######################################

################ RECOVER THE EDGE GENERATORS ########################
for w in Points do
p:=w;
if p in Vertices then
x:=G[Position(Points,w)];
EdgeGenerators[Position(Vertices,p)]:=x;
fi;
od;
################ EDGE GENERATORS RECOVERED ##########################

################ READ HASSE DIAGRAM #################################
Exec(Concatenation(POLYMAKE_PATH,tmpInlog," \"HASSE_DIAGRAM -> FACES\" >",tmp2Inlog));
#Exec(Concatenation(POLYMAKE_PATH,tmpInlog," \"HASSE_DIAGRAM -> DIMS\" >",tmp4Inlog));
Exec(Concatenation(POLYMAKE_PATH,tmpInlog," \"F_VECTOR\" >",tmp4Inlog));
RemoveFile(tmpInlog);
Filename(tmpDir,tmp);
RemoveFile(tmp);

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

input:=InputTextFile(tmp4Inlog);
tmp:=ReadLine(input);
tmp:=ReadLine(input);
#Index:=IndexFn(tmp);
tmp:= IndexFn(tmp);           #Modified July 2019
Index:=[1];                   #because Polymake has
for i in [1..Length(tmp)] do  #discontinued the DIMS
Add(Index,Index[i]+tmp[i]);   #property.
od;                           #
RemoveFile(tmp4Inlog);



input:=InputTextFile(tmp2Inlog);
tmp:=ReadLine(input);
tmp:=ReadLine(input);

while Length(tmp)>2 do
Append(Faces, FacesFn(tmp));
tmp:=ReadLine(input);
od;
RemoveFile(tmp2Inlog);

if Length(Faces[1])=1 then

	FacesFinal:=[];
	for i in [1..Length(Index)-1] do
	Append(FacesFinal,[[Index[i]..Index[i+1]-1]]);
	od;
	Append(FacesFinal,[[Index[i+1]]]);
	FacesFinal:=(List(FacesFinal,x->List(x,i->Faces[i])));

else

	FacesFinal:=[[1..Index[1]]];
	for i in [1..Length(Index)-1] do
	Append(FacesFinal,[[Index[i]+1..Index[i+1]]]);
	od;
	FacesFinal:=Reversed(List(FacesFinal,x->List(x,i->Faces[i])));
	
fi;

RemoveFile(tmpInlog);
RemoveFile(tmp2Inlog);
RemoveFile(tmp4Inlog);
RemoveFile(tmp3Inlog);
RemoveFile(Filename(tmpDir," "));
############### HASSE DIAGRAM READ ##################################
return rec(
             generators:=EdgeGenerators,
             hasseDiagram:=FacesFinal,
	     vector:=v);
	     
end);
#####################################################################
