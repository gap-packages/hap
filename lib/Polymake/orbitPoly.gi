#####################################################################
#####################################################################
InstallGlobalFunction(OrbitPolytope,
function(G,v,Props)
local
	Points,p,x,i,w,Dim,tmp,input,
	a,b,c,d,U,V,W;

Dim:=Length(v);
Points:=[];

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
    w:=x*v;
    Append(Points, [w]);
    od;
fi;
######################POINTS CREATED###################


################## CALCULATE HULL POINTS ###############
AppendTo("/tmp/tmpIn.log","POINTS","\n");
for p in Points do
AppendTo("/tmp/tmpIn.log",1);
    for i in [1..Dim] do
    AppendTo("/tmp/tmpIn.log"," ",p[i]);
    od;
AppendTo("/tmp/tmpIn.log","\n");
od;


################# HULL CALCULATED ###################################


if "DIMENSION" in Props or "dimension" in Props then
Exec("polymake /tmp/tmpIn.log DIM > /tmp/tmp2In.log");
input:=InputTextFile("/tmp/tmp2In.log");
tmp:=ReadLine(input);
tmp:=ReadLine(input);
Print("Dimension of orbit polytope is: ", tmp, "\n");
Exec("rm /tmp/tmp2In.log");
fi;

if "VERTREX_DEGREE" in Props or "vertex_degree" in Props then
Exec("polymake /tmp/tmpIn.log VERTEX_DEGREES > /tmp/tmp2In.log");
input:=InputTextFile("/tmp/tmp2In.log");
tmp:=ReadLine(input);
tmp:=ReadLine(input);
Print("Vertex degree in graph of polytope is: ", Rat([tmp[1]]), "\n");
Exec("rm /tmp/tmp2In.log");
fi;



if "VISUAL_GRAPH" in Props or "visual_graph" in Props then
Exec("polymake /tmp/tmpIn.log VISUAL_GRAPH");
fi;

if "SCHLEGEL" in Props or "schlegel" in Props then
Exec("polymake /tmp/tmpIn.log SCHLEGEL");
fi;


if "VISUAL" in Props or "visual" in Props then

if IsPermGroup(G) and Length(v)=4 then
	Exec("rm /tmp/tmpIn.log");
	AppendTo("/tmp/tmpIn.log","POINTS","\n");

	for x in G do
	a:=v[1^x]-v[1];
	b:=v[2^x]-v[2];
	c:=v[3^x]-v[3];
	d:=v[4^x]-v[4];
	U:=2*a-2*b;
	V:=2*c-2*d;
	W:=a+b-c-d;
	AppendTo ("/tmp/tmpIn.log",1," ",U, " ",V," ",W,  "\n");
	od;
fi;
Exec("polymake /tmp/tmpIn.log VISUAL");
fi;

Exec("rm /tmp/tmpIn.log");



end);
#####################################################################
#####################################################################
