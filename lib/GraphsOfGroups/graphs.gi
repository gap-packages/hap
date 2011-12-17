#(C) 2005-2006 Graham Ellis

#A:=Group([(1,2),(1,2,3,4,5)]);SetName(A,"S5");
#B:=Group([(1,2),(1,2,3,4)]);SetName(B,"S4");
#C:=SymmetricGroup(3);SetName(C,"S3");
#CA:=GroupHomomorphismByFunction(C,A,x->x);
#CB:=GroupHomomorphismByFunction(C,B,x->x);
#D:=[A,B,[CA,CB]];

#####################################################################
#####################################################################
InstallGlobalFunction(GraphOfGroupsTest,
function(D)
local
	Boole, x, VertexNames, EdgeNames, RangeNames;

VertexNames:=[];
EdgeNames:=[];
RangeNames:=[];
for x in D do
Boole:=false;
if IsGroup(x) and HasName(x) then 
Boole:=true; 
Append(VertexNames,[Name(x)]);fi;
if IsList(x) then
        if Length(x)=2 then 
        	if IsGroupHomomorphism(x[1])
		and IsGroupHomomorphism(x[2]) then
			if HasName(Range(x[1])) 
	       		and HasName(Source(x[1]))
			and HasName(Source(x[2]))
        		and HasName(Range(x[2])) then
			  if Name(Source(x[1]))=Name(Source(x[2]))
			  then Boole:=true;
			  Append(EdgeNames,[Name(Source(x[1]))]); 
			  Append(RangeNames,[Name(Range(x[1]))]);
			  Append(RangeNames,[Name(Range(x[2]))]);
			  fi;
			fi;
		fi;
	fi;
fi;

if not Boole then return false; fi;
od;

RangeNames:=SSortedList(RangeNames);

if Length(SSortedList(VertexNames))=Length(VertexNames)
and Length(RangeNames)=Length(Intersection(VertexNames,RangeNames))
and  Length(SSortedList(EdgeNames))=Length(EdgeNames) then
return true;
else return false;
fi;

end);
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallGlobalFunction(GraphOfGroupsDisplay,
function(arg)
local

	PositionName,
	D, Boole, Vertices, Edges,x;

#####################################################################
PositionName:=function(L,x);
return PositionProperty(L,n->Name(n)=Name(x));
end;
#####################################################################

D:=arg[1];
Vertices:=[];
Edges:=[];

if not GraphOfGroupsTest(D) then
Print("The list D does not represent a Graph of Groups \n");
return fail; fi;

for x in D do
if IsGroup(x) then Append(Vertices,[x]); fi;
if IsList(x) then Append(Edges,[x]); fi;
od;

################ WRITE TO TMPIN.LOG #################################

AppendTo("/tmp/tmpIn.log"," graph G { \n size=\"4,4\" \n node [shape=circle, style=filled, color=blue] \n edge [style=\"setlinewidth(5)\"] \n");


for x in Vertices do
AppendTo("/tmp/tmpIn.log",PositionName(Vertices,x), "[label=\" ", Name(x), "\",fontsize=10]\n");
od;

for x in Edges do
AppendTo("/tmp/tmpIn.log",
PositionName(Vertices,Range(x[1]))," -- ",
PositionName(Vertices,Range(x[2])), "[label=\" ",Name(Source(x[1])),"\",fontsize=10,color=brown] \n");
od;

AppendTo("/tmp/tmpIn.log","} \n");
################ WRITTEN ############################################

Exec("neato -Tgif /tmp/tmpIn.log > /tmp/basic.gif");

if Length(arg)=1 then
Exec("mozilla /tmp/basic.gif");
Sleep(2);
Exec("rm /tmp/tmpIn.log; rm /tmp/basic.gif");

else
AppendTo("/tmp/tmpIn2.log", "Browser=",arg[2],"\n");
AppendTo("/tmp/tmpIn2.log","$Browser /tmp/basic.gif");
Exec("chmod a+x /tmp/tmpIn2.log; /tmp/tmpIn2.log");
Exec("rm /tmp/tmpIn.log; rm /tmp/basic.gif; rm /tmp/tmpIn2.log;");
fi;

end);
#####################################################################
#####################################################################
