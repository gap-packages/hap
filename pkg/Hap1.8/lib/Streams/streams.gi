HAPchildren:=[];


####################################################
####################################################
InstallGlobalFunction(ChildProcess,
function(arg)
local d,f,s,Name;
d:= DirectoryCurrent();
f := Filename(DirectoriesSystemPrograms(), "gap");
if Length(arg)>0 then
###
###
###
#s := InputOutputLocalProcess(d,f,["-T","-q","-b",'-L arg[1] "]);
###This does not work!!!!!!!!!
###
###
else
s := InputOutputLocalProcess(d,f,["-T","-q","-b"]);
fi;
WriteLine(s,"SizeScreen([255,24]);");

Name:=Concatenation(["/tmp/HAPchild",String(1+Length(HAPchildren))]);
Add(HAPchildren,Name);
PrintTo(Name,"HAPchildToggle\:=true;");

return rec(
	stream:=s,
	name:=Name,
	number:=Length(HAPchildren));
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildClose,
function(s);
CloseStream(s.stream);
PrintTo(s.name,"HAPchildToggle\:=false;");
end);
####################################################
####################################################

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
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildCommand,
function(command,s)
local i;
PrintTo(s.name,"HAPchildToggle\:=false;");

WriteLine(s.stream,command);;
i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;

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

i:=(ReadAllLine(s.stream,true));
i:=(ReadAllLine(s.stream,true));

output:="";
while true do
if not i="\"STOP\"\n" then 
Append(output,i);
i:=(ReadAllLine(s.stream,true)); 
else break; fi;
od;

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
HAPchildToggle:=false;

#####################################################
#####################################################
InstallGlobalFunction(NextAvailableChild,
function(L)
local
	s;

while true do

	for s in L do
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
