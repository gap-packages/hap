#To run this file just something like the following (the -A is important). 
#    gap -l '/home/graham/;' -A -r update.gi

root:=DirectoriesPackageLibrary("HAP");
rootversion:=Filename(root,"../version");
root:=Filename(root,".");
root:=Concatenation(root{[1..Length(root)-5]},"doc/");

MakeGAPDocDoc(root,"newHapMan.xml",[],"HAP","MathJax");
i:=InputTextFile(rootversion);;
version:=ReadAll(i);;
version:=NormalizedWhitespace(version);;

hapvariables:= PackageVariablesInfo("HAP",version);;
h:=List(hapvariables,x->x[2]);;
hapvariables:=[];;
for i in [1..Length(h)] do
for j in [1..Length(h[i])] do
Add(hapvariables,h[i][j][1]);
od;od;


input:=InputTextFile(Concatenation(root,"chapInd.txt"));
hapindex:=ReadAll(input);;
hapindex:=ReplacedString(hapindex,"[2X","@");;
hapindex:=SplitString(hapindex,['@']);;
hapindex:=hapindex{[2..Length(hapindex)]};;
hapindex:=List(hapindex,s->SplitString(s,['\033']));;
hapindex:=List(hapindex,x->x[1]);;

hapundocumented:=Filtered(hapvariables,x->not x[1] in hapindex);;

s:="<Chapter><Heading>HAP variables that are not yet documented</Heading>\n<Section><Heading> &nbsp;</Heading><Br/>";

PrintTo(Concatenation(root,"Undocumented.xml"),s);

for x in hapundocumented do
s:=ReplacedString(x[1],"<","\<");
s:=Concatenation("<C>",s,"</C><Br/>\n");

AppendTo(Concatenation(root,"Undocumented.xml"),s);
od;

s:="</Section>\n</Chapter>\n";
AppendTo(Concatenation(root,"Undocumented.xml"),s);

MakeGAPDocDoc(root,"newHapMan.xml",[],"HAP","MathJax");
QUIT;
