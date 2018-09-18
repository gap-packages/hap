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
        	#if IsGroupHomomorphism(x[1])
		#and IsGroupHomomorphism(x[2]) then
                if HasSource(x[1])
                and HasSource(x[2]) then

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

	PositionName, tmpDir, tmpInlog, tmpIn2log, basicgif,
	D, Boole, Vertices, Edges,x,AppendTo,PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;

tmpDir:=DirectoryTemporary();
tmpInlog:=Filename(tmpDir,"tmpIn.log");
tmpIn2log:=Filename(tmpDir,"tmpIn2.log");
basicgif:=Filename(tmpDir,"basic.gif");

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

AppendTo(tmpInlog," graph G { \n size=\"10,10\" \n node [shape=circle, style=filled,  color=gray] \n edge [style=\"setlinewidth(2)\"] \n");


for x in Vertices do
AppendTo(tmpInlog,PositionName(Vertices,x), "[label=\" ", Name(x), "\",fontsize=8]\n");
od;

for x in Edges do
AppendTo(tmpInlog,
PositionName(Vertices,Range(x[1]))," -- ",
PositionName(Vertices,Range(x[2])), "[label=\" ",Name(Source(x[1])),"\",fontsize=8,color=black] \n");
od;

AppendTo(tmpInlog,"} \n");
################ WRITTEN ############################################

Exec(Concatenation(NEATO_PATH,"-Tgif ",tmpInlog ," > ",basicgif));

if Length(arg)=1 then
Exec(Concatenation(DISPLAY_PATH, basicgif));
Sleep(2);
Exec(Concatenation("rm ",tmpInlog, "; rm ",basicgif));

else
AppendTo(tmpIn2log, "Browser=",arg[2],"\n");
AppendTo(tmpIn2log,"$Browser ",basicgif);
Exec(Concatenation("chmod a+x ",tmpIn2log,"; ",tmpIn2log));
Exec(Concatenation("rm ",tmpInlog,"; rm ",basicgif,"; rm ",tmpIn2log,";"));
fi;

end);
#####################################################################
#####################################################################
