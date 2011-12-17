####################################################
####################################################
InstallGlobalFunction(ChildProcess,
function(arg)
local d,f,s;
d:= DirectoryCurrent();
f := Filename(DirectoriesSystemPrograms(), "gap");
if Length(arg)>0 then
###
###
###
s := InputOutputLocalProcess(d,f,["-T -q -b -L arg[1] "]);
###This does not work!!!!!!!!!
###
###
else
s := InputOutputLocalProcess(d,f,["-T -q -b  "]);
fi;
WriteLine(s,"SizeScreen([255,24]);");
return s;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildClose,
function(s);
CloseStream(s);
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildFunction,
function(func,s)
local i;
WriteLine(s,"\"STOP\";");;
i:=0;
while true do
if not i="\"STOP\"\n" then
i:=ReadAllLine(s,true);
else break; fi;
od;

i:=Concatenation(["\"START\";xyxy:=",func,"\"STOP\";"]);
WriteLine(s,i);;
#WriteAll(s,i);;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildCommand,
function(command,s)
local i;
WriteLine(s,command);;
#WriteAll(s,command);;

end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildRead,
function(s)
local i,output;

while true do
i:=ReadAllLine(s,true); 
if  i="\"START\"\n" then 
fi;break; 
od;

i:=(ReadAllLine(s,true));
i:=(ReadAllLine(s,true));

output:="";
while true do
if not i="\"STOP\"\n" then 
Append(output,i);
i:=(ReadAllLine(s,true)); 
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


