HAPchildren:=[];


####################################################
####################################################
InstallGlobalFunction(ChildProcess,
function(arg)
local d,f,s,Name,Remote,Computer,ST;
d:= DirectoryCurrent();

if Length(arg)=0 then 
Remote:=false;
Computer:="localhost";
f := Filename(DirectoriesSystemPrograms(), "gap");
s := InputOutputLocalProcess(d,f,["-T","-q","-b"]);
fi;
if Length(arg)=1 then
if not IsString(arg[1]) then
Remote:=false;
Computer:="localhost";
f := Filename(DirectoriesSystemPrograms(), "gap");

s := InputOutputLocalProcess(d,f,Concatenation(["-T","-q","-b"],arg[1]));
else
Remote:=true;
Computer:=arg[1];
f:=Filename(DirectoriesSystemPrograms(),"ssh");
s := InputOutputLocalProcess(d,f,[arg[1],"gap","-f","-T","-q","-b"]);
fi;fi;

if Length(arg)=2 then
Remote:=true;
Computer:=arg[1];
f:=Filename(DirectoriesSystemPrograms(),"ssh");
s := InputOutputLocalProcess(d,f,Concatenation([arg[1],"gap","-f","-T","-q","-b"],arg[2]));
fi;

if not Remote then
Name:=Filename(DirectoryTemporary(),"name");
else
ST:=String(Random([1..100000]));
Name:=Concatenation("/tmp/",ST);
while IsExistingFile(Name) do
ST:=String(Random([1..100000]));
Name:=Concatenation("/tmp/",ST);
od;
Name:=Filename(Directory("/tmp"),ST);
fi;

Add(HAPchildren,Name);
PrintTo(Name,"HAPchildToggle\:=true;");

return rec(
	stream:=s,
	name:=Name,
	number:=Length(HAPchildren),
	remote:=Remote,
	computer:=Computer);
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildClose,
function(s)
local i;
CloseStream(s.stream);
if not s!.remote then
i:=Concatenation("rm -r ",s!.name{[1..Length(s!.name)-4]});
else
i:=Concatenation(["rm ",s!.name]);
fi;

Exec(i);

end);
####################################################
####################################################

HAPchildFunctionToggle:=true;
####################################################
####################################################
InstallGlobalFunction(ChildFunction,
function(func,s)
local i;

PrintTo(s.name,"HAPchildToggle\:=false;");
WriteLine(s.stream,"\"STOP\";");;
i:=0;
while true do 			#FLUSH THE CHILD
if not i="\"STOP\"\n" then
i:=ReadAllLine(s.stream,true);
else break; fi;
od;

i:=Concatenation(["\"START\";xyxy:=",func,"\"STOP\";"]);
WriteLine(s.stream,i);;

i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;
HAPchildFunctionToggle:=false;
end);
####################################################
####################################################


####################################################
####################################################
InstallGlobalFunction(ChildCommand,
function(command,s)
local i;

PrintTo(s.name,"HAPchildToggle\:=false;");
Append(command,";");

if HAPchildFunctionToggle then
while not ReadAllLine(s.stream)=fail do  #Flush stream
od;
fi;

WriteLine(s.stream,command);;
i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;

#if HAPchildFunctionToggle then
#while not ReadAllLine(s.stream)=fail do  #Flush stream
#od;
#fi;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildRead,
function(s)
local i,output;

while true do
i:=ReadAllLine(s.stream,true); 
if  i="\"START\"\n" then 
fi;break; 
od;

if not s!.remote then i:=(ReadAllLine(s.stream,true));fi;
i:=(ReadAllLine(s.stream,true));

output:="";
while true do
if not i="\"STOP\"\n" then 
Append(output,i);
i:=(ReadAllLine(s.stream,true)); 
else break; fi;
od;

HAPchildFunctionToggle:=true;
return output;

end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(ChildReadEval,
function(s);
return EvalString(ChildRead(s));
end);
#####################################################
#####################################################


#####################################################
#####################################################
InstallGlobalFunction(ChildPut,
function(X,Name,s)
local i,fle,flebatch,fleRemote;

NextAvailableChild([s]); #Don't start if s is busy.

PrintTo(s.name,"HAPchildToggle\:=false;");

if not s!.remote then
fle:=Filename(DirectoryTemporary(),"fle");
else

fle:=Concatenation("/tmp/",String(Random([1..100000])));
while IsExistingFile(fle) do
fle:=Concatenation("/tmp/",String(Random([1..100000])));
od;
fleRemote:=Concatenation(fle,"Remote");
flebatch:=Concatenation(fle,"Batch");
PrintTo(flebatch,Concatenation(["put ",fle, " ",fleRemote ]) );
fi;


PrintTo(fle,Concatenation(Name,":="));
AppendTo(fle,X);
AppendTo(fle,";");


if s!.remote then
        i:=Concatenation(["sftp -b ",flebatch," ",s!.computer]);
        Exec(i);
	i:=Concatenation(["Exec(\"mv ",fleRemote," ",fle," \")"]);
        ChildCommand(i,s);
        fi;

i:=Concatenation("Read(\"",fle,"\");");
ChildCommand(i,s);
if not s!.remote then
i:=Concatenation("Exec(\"rm -r ",fle{[1..Length(fle)-4]},"\");");
else 
i:=Concatenation("Exec(\"rm ",fle,"\");");
fi;
ChildCommand(i,s);
if s!.remote then 
i:=Concatenation("Exec(\"rm ",flebatch,"\");");
ChildCommand(i,s);
fi;


end);
#####################################################
#####################################################

HAP_XYXYXYXY:=0;
#####################################################
#####################################################
InstallGlobalFunction(ChildGet,
function(X,s)
local i,ST, fle, fleLocal, batchfile;

NextAvailableChild([s]); #Don't start if s is busy.

PrintTo(s.name,"HAPchildToggle\:=false;");


if not s!.remote then
fle:=Filename(DirectoryTemporary(),"fle");
else

ST:=String(Random([1..100000]));
fle:=Concatenation("/tmp/",ST);
while IsExistingFile(fle) do
ST:=String(Random([1..100000]));
fle:=Concatenation("/tmp/",ST);
od;
#fle:=Filename(Directory("/tmp"),ST);
fleLocal:=Concatenation(fle,"Local");
batchfile:=Concatenation(fle,"Batch");
PrintTo(batchfile,Concatenation(["get ",fle, " ",fleLocal," ; rm", fle]) );

fi;


i:=Concatenation("PrintTo(\"",fle,"\", \"HAP_XYXYXYXY:=\");");
WriteLine(s.stream,i);

i:=Concatenation("AppendTo(\"",fle,"\"," ,  X,");");

ChildCommand(i,s);
NextAvailableChild([s]);

PrintTo(s.name,"HAPchildToggle\:=false;");

if s!.remote then 
i:=Concatenation(["sftp -b ",batchfile," ",s!.computer]);
Exec(i);
i:=Concatenation(["rm ",batchfile]);
Exec(i);
i:=Concatenation(["Exec(\"mv ",fleLocal," ",fle," \");"]);
ChildCommand(i,s);
fi;

AppendTo(fle,";");

Read(fle);
i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;
if s!.remote then
Exec(Concatenation(["rm ",fle]));
else
Exec(Concatenation("rm -r ",fle{[1..Length(fle)-4]}));
fi;


return HAP_XYXYXYXY;

end);
#####################################################
#####################################################






HAPchildToggle:=false;

#####################################################
#####################################################
InstallGlobalFunction(NextAvailableChild,
function(L)
local
	s,i,localname;

while true do

	for s in L do
	if s!.remote then 
        localname:=Concatenation(s!.name,"Local");
	i:=Concatenation(["ssh ",s!.computer," less ",s!.name," > ",localname]);
	Exec(i);
        Exec(Concatenation("mv ",localname," ",s!.name));
	fi;
	Read(s!.name);
	if HAPchildToggle=true then
	return s; break; fi;
	od;

	Exec("sleep 0.01");
od;
end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(IsAvailableChild,
function(s)
local x;

Read(s!.name);

if HAPchildToggle then return true;
else return false;
fi;

end);
#####################################################
#####################################################

#####################################################
#####################################################
InstallGlobalFunction(ParallelList,
function(lst,func,processes)
local
	ind,s,i,j,kk,answer,final,cnt,f,prm1,prm2;

answer:=[];
prm1:=[];
prm2:=[];
cnt:=1;

for i in lst{[1..Length(processes)]} do
f:=Concatenation([func,"(",String(i),");"]);
ChildFunction(f,processes[cnt]);
Add(prm1,processes[cnt].number);
cnt:=cnt+1;
od;

for i in lst{[1+Length(processes)..Length(lst)]} do
s:=NextAvailableChild(processes);
Add(answer,ChildReadEval(s));
Add(prm2,s.number);
f:=Concatenation([func,"(",String(i),");"]);
ChildFunction(f,s);
Add(prm1,s.number);

od;

for s in processes do
Add(answer,ChildReadEval(s));
Add(prm2,s.number);
od;

#Reorder the answer
ind:=SSortedList(prm1);
final:=[];
cnt:=[];
for i in ind do
cnt[i]:=0;
od;

for i in [1..Length(prm1)] do
cnt[prm1[i]]:=cnt[prm1[i]] + 1;
final[i]:=answer[PositionNthOccurrence(prm2,prm1[i],cnt[prm1[i]])];

od;

return final;
end);
#####################################################
#####################################################
