#####################################################################
#####################################################################
InstallGlobalFunction(OrbitPolytope,
function(G,v,Props)
local
	Points,p,x,i,w,Dim,tmp,input, tmpdir, tmpIn, tmp2dir, tmp2In,
	a,b,c,d,U,V,W,AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

Dim:=Length(v);
Points:=[];

tmpdir := DirectoryTemporary();;
tmpIn:=Filename( tmpdir , "tmpIn.log" );
tmp2dir := DirectoryTemporary();;
tmp2In:=Filename( tmpdir , "tmp2In.log" );


############################CREATE POINTS###############
if IsPermGroup(G) then

for x in G do
w:=[];
    for i in [1..Dim] do
    Append(w,[v[i^x]]);
    od;
Append(Points, [w]);
od;
						
else

for x in G do
    w:=v*x;
    Append(Points, [w]);
    od;
fi;
######################POINTS CREATED###################


################## CALCULATE HULL POINTS ###############
AppendTo(tmpIn,"POINTS","\n");
for p in Points do
AppendTo(tmpIn,1);
    for i in [1..Dim] do
    AppendTo(tmpIn," ",p[i]);
    od;
AppendTo(tmpIn,"\n");
od;


################# HULL CALCULATED ###################################


if "DIMENSION" in Props or "dimension" in Props then
Exec(Concatenation(POLYMAKE_PATH, tmpIn," DIM > ", tmp2In));
input:=InputTextFile(tmp2In);
tmp:=ReadLine(input);
tmp:=ReadLine(input);
Print("Dimension of orbit polytope is: ", tmp, "\n");
Exec(Concatenation("rm ", tmp2In));
fi;

if "VERTEX_DEGREES" in Props or "vertex_degrees" in Props then
Exec(Concatenation(POLYMAKE_PATH,tmpIn,  " VERTEX_DEGREES > ",tmp2In));
input:=InputTextFile(tmp2In);
tmp:=ReadLine(input);
tmp:=ReadLine(input);
Print("Vertex degree in graph of polytope is: ", Rat(Concatenation([tmp{[1..Position(tmp,' ')-1]}])), "\n");
#Exec(Concatenation("rm ",tmp2In));
fi;



if "VISUAL_GRAPH" in Props or "visual_graph" in Props then
Exec(Concatenation(POLYMAKE_PATH,tmpIn ," VISUAL_GRAPH"));
fi;

if "SCHLEGEL" in Props or "schlegel" in Props then
Exec(Concatenation(POLYMAKE_PATH,tmpIn ," SCHLEGEL"));
fi;


if "VISUAL" in Props or "visual" in Props then

if IsPermGroup(G) and Length(v)=4 then
	Exec(Concatenation("rm ",tmpIn));
	AppendTo(tmpIn,"POINTS","\n");

	for x in G do
	a:=v[1^x]-v[1];
	b:=v[2^x]-v[2];
	c:=v[3^x]-v[3];
	d:=v[4^x]-v[4];
	U:=2*a-2*b;
	V:=2*c-2*d;
	W:=a+b-c-d;
	AppendTo (tmpIn,1," ",U, " ",V," ",W,  "\n");
	od;
fi;
Exec(Concatenation(POLYMAKE_PATH,tmpIn ," VISUAL"));
fi;

#Exec(Concatenation("rm ",tmpIn));



end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod(Display,
"method for displaying convex hulls of sets of points",
[IsPolymakeObject],
function(F)
local
        Points,p,x,i,w,Dim,tmp,input, tmpdir, tmpIn, tmp2dir, tmp2In,
        a,b,c,d,U,V,W,AppendTo, PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

Dim:=Length(Polymake(F,"F_VECTOR"));
Points:=Polymake(F,"VERTICES");

tmpdir := DirectoryTemporary();;
tmpIn:=Filename( tmpdir , "tmpIn.log" );
tmp2dir := DirectoryTemporary();;
tmp2In:=Filename( tmpdir , "tmp2In.log" );




################## CALCULATE HULL POINTS ###############
AppendTo(tmpIn,"POINTS","\n");
for p in Points do
AppendTo(tmpIn,1);
    for i in [1..Dim] do
    AppendTo(tmpIn," ",p[i]);
    od;
AppendTo(tmpIn,"\n");
od;


################# HULL CALCULATED ###################################


Exec(Concatenation(POLYMAKE_PATH,tmpIn ," VISUAL"));

#Exec(Concatenation("rm ",tmpIn));

end);
#####################################################################
#####################################################################

