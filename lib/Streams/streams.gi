HAPchildren:=[];
HAP_XYXYXYXY:=0;

####################################################
####################################################
InstallGlobalFunction(ChildProcess,
function(arg)
local d,f,s,Name,Remote,Computer,ST,AppendTo,PrintTo,tmpdir;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;
tmpdir := DirectoryTemporary();;

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
#Name:=Concatenation("/tmp/",ST);
Name:=Concatenation(tmpdir,ST);
while IsExistingFile(Name) do
ST:=String(Random([1..100000]));
#Name:=Concatenation("/tmp/",ST);
Name:=Concatenation(tmpdir,ST);
od;
#Name:=Filename(Directory("/tmp"),ST);
Name:=Filename(Directory(tmpdir),ST);
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

ChildCreate:=ChildProcess;
MakeReadOnlyGlobal("ChildCreate");

####################################################
####################################################
InstallGlobalFunction(ChildClose,
function(s)
local i;
CloseStream(s.stream);
if not s!.remote then
#i:=Concatenation("rm -r ",s!.name{[1..Length(s!.name)-4]});
RemoveFile(s!.name{[1..Length(s!.name)-4]});
else
#i:=Concatenation(["rm ",s!.name]);
RemoveFile(s!.name);
fi;

#Exec(i);

end);
####################################################
####################################################

HAPchildFunctionToggle:=true;
####################################################
####################################################
InstallGlobalFunction(ChildFunction,
function(func,s)
local i,AppendTo,PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;


PrintTo(s.name,"HAPchildToggle\:=false;");
WriteLine(s.stream,"\"STOP\";");;
i:=0;     
while true do 			#FLUSH THE CHILD
if not i="\"STOP\"\n" then
i:=ReadAllLine(s.stream,true);
i:=StripEscapeSequences(i);
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

ChildKill:=ChildClose;
MakeReadOnlyGlobal("ChildKill");

####################################################
####################################################
InstallGlobalFunction(ChildRestart,
function(s)
local r;
ChildKill(s);
r:=ChildCreate();
s.name:=r.name; s.stream:=r.stream;
end);
####################################################
####################################################


####################################################
####################################################
InstallGlobalFunction(ChildCommand,
function(arg)
local command, s, i,AppendTo,PrintTo;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;


command:=arg[1];
s:=arg[2];
if Length(arg)=2 then
NextAvailableChild([s]);
fi;

PrintTo(s.name,"HAPchildToggle\:=false;");
Append(command,";");

WriteLine(s.stream,"true;");
if HAPchildFunctionToggle then
while not ReadAllLine(s.stream)=fail do  #Flush stream
od;
fi;

WriteLine(s.stream,command);;
i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;

return true;
end);
####################################################
####################################################

####################################################
####################################################
InstallGlobalFunction(ChildRead,
function(s)
local i,output;

while true do
i:=StripEscapeSequences(ReadAllLine(s.stream,true)); 
if  i="\"START\"\n" then 
break; fi;
od;

if not s!.remote then i:=(ReadAllLine(s.stream,true));fi;
#i:=(ReadAllLine(s.stream,true));

output:="";
while true do
if not i="\"STOP\"\n" then 
Append(output,i);
i:=(ReadAllLine(s.stream,true)); 
i:=StripEscapeSequences(i);
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
local i,fle,flebatch,fleRemote,AppendTo,PrintTo, tmpdir;

if IsHapResolution(X) then
ChildPutObj(X,Name,s);
return true;
fi;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;
tmpdir := DirectoryTemporary();;


NextAvailableChild([s]); #Don't start if s is busy.

PrintTo(s.name,"HAPchildToggle\:=false;");

if not s!.remote then
fle:=Filename(DirectoryTemporary(),"fle");
else

#fle:=Concatenation("/tmp/",String(Random([1..100000])));
fle:=Concatenation(tmpdir,String(Random([1..100000])));
while IsExistingFile(fle) do
#fle:=Concatenation("/tmp/",String(Random([1..100000])));
fle:=Concatenation(tmpdir,String(Random([1..100000])));
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
        ChildCommand(i,s,true);
        fi;

i:=Concatenation("Read(\"",fle,"\");");
ChildCommand(i,s,true);
if not s!.remote then
#i:=Concatenation("Exec(\"rm -r ",fle{[1..Length(fle)-4]},"\");");
i:=Concatenation("RemoveFile(\"",fle{[1..Length(fle)-4]},"\");");
else 
#i:=Concatenation("Exec(\"rm ",fle,"\");");
i:=Concatenation("RemoveFile(\"",fle{[1..Length(fle)]},"\");");
fi;
ChildCommand(i,s,true);
if s!.remote then 
#i:=Concatenation("Exec(\"rm ",flebatch,"\");");
i:=Concatenation("RemoveFile(\"",flebatch,"\");");
ChildCommand(i,s,true);
fi;


end);
#####################################################
#####################################################

HAP_XYXYXYXY:=0;
#####################################################
#####################################################
InstallGlobalFunction(ChildGet,
function(X,s)
local i,ST, fle, fleLocal, batchfile, AppendTo, PrintTo, tmpdir;

AppendTo:=HAP_AppendTo;
PrintTo:=HAP_PrintTo;
tmpdir:=DirectoryTemporary();
NextAvailableChild([s]); #Don't start if s is busy.

PrintTo(s.name,"HAPchildToggle\:=false;");


if not s!.remote then
fle:=Filename(DirectoryTemporary(),"fle");
else

ST:=String(Random([1..100000]));
#fle:=Concatenation("/tmp/",ST);
fle:=Concatenation(tmpdir,ST);
while IsExistingFile(fle) do
ST:=String(Random([1..100000]));
#fle:=Concatenation("/tmp/",ST);
fle:=Concatenation(tmpdir,ST);
od;
#fle:=Filename(Directory("/tmp"),ST);
fleLocal:=Concatenation(fle,"Local");
batchfile:=Concatenation(fle,"Batch");
PrintTo(batchfile,Concatenation(["get ",fle, " ",fleLocal," ; rm", fle]) );

fi;


i:=Concatenation("PrintTo(\"",fle,"\", \"HAP_XYXYXYXY:=\");");
WriteLine(s.stream,i);

i:=Concatenation("AppendTo(\"",fle,"\"," ,  X,");");

ChildCommand(i,s,true);
NextAvailableChild([s]);

PrintTo(s.name,"HAPchildToggle\:=false;");

if s!.remote then 
i:=Concatenation(["sftp -b ",batchfile," ",s!.computer]);
Exec(i);
i:=Concatenation(["rm ",batchfile]);
Exec(i);
i:=Concatenation(["Exec(\"mv ",fleLocal," ",fle," \");"]);
ChildCommand(i,s,true);
fi;

AppendTo(fle,";");

Read(fle);
i:=Concatenation(["PrintTo(\"",String(s.name),"\",\"HAPchildToggle\:=true;\");"]);
WriteLine(s.stream,i);;
if s!.remote then
#Exec(Concatenation(["rm ",fle]));
RemoveFile(fle);
else
#Exec(Concatenation("rm -r ",fle{[1..Length(fle)-4]}));
RemoveFile(fle{[1..Length(fle)-4]});
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
#Print(s.number,"\n");
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
#Print(prm2,"\n");
return final;
end);
#####################################################
#####################################################

########################################################
########################################################
InstallGlobalFunction(ChildTransfer,
function(arg)
local X,Y,S,T,fle,i;

X:=arg[1]; if IsString(X) then X:=[X]; fi;
Y:=arg[2]; if IsString(Y) then Y:=[Y]; fi;
S:=arg[3]; if IsRecord(S) then S:=[S]; fi;
T:=arg[4]; if IsRecord(T) then T:=[T]; fi;

if Length(X)<>Length(Y) or Length(X)<>Length(S) or Length(X)<>Length(T)
then Print("Error: input arguments must be lists of equal length.\n");
return fail;
fi;
#Transfer object X[i] on child S[i] to child T[i] and call it Y[i]. 
#Here X and Y are lists of strings.
##THIS ONLY WORKS IF THE TWO CHILDREN ARE ON THE SAME MACHINE
fle:=[];
#################################
for i in [1..Length(X)] do
fle[i]:=Filename(DirectoryTemporary(),"fle");
ChildCommand(Concatenation("HAPPrintTo(\"",fle[i],"\",",String(X[i]),"\);")  ,S[i]);
od;
List(S,s->NextAvailableChild([s]));
for i in [1..Length(X)] do
ChildCommand(Concatenation(String(Y[i]),":=HAPRead(\"",fle[i],"\");"), T[i]);
od;
#NextAvailableChild([t]);
#RemoveFile(fle);
#od;
#################################

end);
########################################################
########################################################

########################################################
########################################################
InstallMethod(ChildPutObj,
"Put a resolution onto a child process",
[IsHapResolution,IsString,IsRecord],
function(X,Y,s)
local fle;

fle:=Filename(DirectoryTemporary(),"fle");

HAPPrintTo(fle,X);
ChildCommand(Concatenation(String(Y),":=HAPRead(\"",fle,"\");"), s);
NextAvailableChild([s]);
RemoveFile(fle);
end);
########################################################
########################################################


########################################################
########################################################
InstallOtherMethod(ChildPutObj,
"Put a resolution onto a child process",
[IsHapResolution,IsString,IsList],
function(X,Y,S)
local fle,s;

fle:=Filename(DirectoryTemporary(),"fle");

HAPPrintTo(fle,X);
for s in S do
ChildCommand(Concatenation(String(Y),":=HAPRead(\"",fle,"\");"), s);
od;
for s in S do
NextAvailableChild([s]);
od;
RemoveFile(fle);
end);
########################################################
########################################################

########################################################
########################################################
InstallOtherMethod(ChildPutObj,
"Put a matrix onto a child process",
[IsTable,IsString,IsRecord],
#[IsMatrix,IsString,IsRecord],
function(X,Y,s)
local fle;

fle:=Filename(DirectoryTemporary(),"fle");

PrintTo(fle,Concatenation(Y,":="));
AppendTo(fle,X);
AppendTo(fle,";");
ChildCommand(Concatenation("Read(\"",fle,"\");"), s);
NextAvailableChild([s]);
RemoveFile(fle);
end);
########################################################
########################################################

########################################################
########################################################
InstallOtherMethod(ChildPutObj,
"Put a matrix onto a child process",
[IsTable,IsString,IsList],
#[IsMatrix,IsString,IsList],
function(X,Y,S)
local fle,s;

fle:=Filename(DirectoryTemporary(),"fle");

PrintTo(fle,Concatenation(Y,":="));
AppendTo(fle,X);
AppendTo(fle,";");
for s in S do
ChildCommand(Concatenation("Read(\"",fle,"\");"), s);
od;
for s in S do
NextAvailableChild([s]);
od;
RemoveFile(fle);
end);
########################################################
########################################################

########################################################
########################################################
InstallOtherMethod(ChildPutObj,
"Put a filtered pure cubical complex onto a child process",
[IsHapFilteredPureCubicalComplex,IsString,IsRecord],
function(X,Y,s)
local fle, cmd, HAP_AA, HAP_BB;

fle:=Filename(DirectoryTemporary(),"fle");
HAP_AA:=Random([1..1000000]);
HAP_AA:=Concatenation("HAP_AA",String(HAP_AA));
HAP_BB:=Concatenation("HAP_BB",String(HAP_AA));
ChildPutObj(X!.binaryArray,HAP_AA,s);
ChildPutObj(X!.filtration,HAP_BB,s);
cmd:=Concatenation(Y,":=FilteredPureCubicalComplex(",HAP_AA,",",HAP_BB,");");
ChildCommand(cmd,s);
NextAvailableChild([s]);

RemoveFile(fle);
end);
########################################################
########################################################

########################################################
########################################################
InstallOtherMethod(ChildPutObj,
"Put a filtered pure cubical complex onto a child process",
[IsHapFilteredPureCubicalComplex,IsString,IsList],
function(X,Y,S)
local fle, cmd, HAP_AA, HAP_BB, s;

fle:=Filename(DirectoryTemporary(),"fle");
HAP_AA:=Random([1..1000000]);
HAP_AA:=Concatenation("HAP_AA",String(HAP_AA));
HAP_BB:=Concatenation("HAP_BB",String(HAP_AA));
for s in S do
ChildPutObj(X!.binaryArray,HAP_AA,s);
ChildPutObj(X!.filtration,HAP_BB,s);
cmd:=Concatenation(Y,":=FilteredPureCubicalComplex(",HAP_AA,",",HAP_BB,");");
ChildCommand(cmd,s);
od;

RemoveFile(fle);
end);
########################################################
########################################################


#########################################################
########################################################
InstallMethod(ChildGetObj,
"Get a resolution from a child process",
[IsString,IsRecord],
function(X,s)
local fle,Y;

fle:=Filename(DirectoryTemporary(),"fle");

ChildCommand(Concatenation("HAPPrintTo(\"",fle,"\",",X,");"),s);
NextAvailableChild([s]);
Y:=HAPRead(fle);
RemoveFile(fle);
return Y;
end);
########################################################
########################################################

##########################################################################
##########################################################################
InstallGlobalFunction(ParallelPersistentBettiNumbers,
function(arg)
local FF,NN,p,S, F,N,B,terms,children,i,j,s,cmd1,cmd2,cmd3,cmd4,cmd5,
      cmd6,cmd7,cmd8,cmd, maps,m,T,TT,lst,P;
FF:=arg[1];
NN:=arg[2];
if Length(arg)=3 then p:=2; S:=arg[3]; else p:=arg[3]; S:=arg[4];fi;
F:=ContractedComplex(FF);
N:=NN;
N:=Filtered(N,x->x>0);
ChildPutObj(F,"F",S);;
terms:=SSortedList(Flat(F!.filtration));;
RemoveSet(terms,0);;
children:=[];;

for i in terms do
  s:=NextAvailableChild(S);;
  children[i]:=s;;
  cmd1:="F:=FilteredPureCubicalComplexToCubicalComplex(F);";
  cmd2:="F:=FilteredCubicalComplexToFilteredRegularCWComplex(F);";
  cmd3:=Concatenation("Y",String(i),":=FiltrationTerm(F,",String(i),");");
  cmd4:=Concatenation("CriticalCells(Y",String(i),");");
  cmd:=Concatenation(cmd1,cmd2,cmd3,cmd4);
  ChildCommand(cmd,s);
od;

for i in [1..Length(terms)-1] do
  ChildTransfer(Concatenation("Y",String(terms[i+1])),Concatenation("YY",String(terms[i])),children[terms[i+1]],children[terms[i]]);
od;

for i in [1..Length(terms)-1] do
  cmd1:="map:=function(n,j); return j; end;\n";
  cmd2:=Concatenation("cwmap:=Objectify( HapRegularCWMap, rec( source:=Y",String(terms[i]),", target:=YY",String(terms[i]),", mapping:=map));\n");
  cmd3:="L:=ChainMapOfRegularCWMap(cwmap);\n";
  cmd4:=Concatenation("Eq:=ChainComplexEquivalenceOfRegularCWComplex(Y",String(terms[i]),");\n");
  cmd5:=Concatenation("EEq:=ChainComplexEquivalenceOfRegularCWComplex(YY",String(terms[i]),");\n");
  cmd6:="f:=Compose( EEq[1], Compose(L,Eq[2]) );";
  cmd7:=Concatenation("ff",String(terms[i]),":=[];");
  cmd8:=Concatenation("for n in ",String(N),"do\n ff",String(terms[i]),"[n+1]:=HomologyVectorSpace(TensorWithIntegersModP(f,",String(p),"),n);\n od;\n");
  cmd:=Concatenation(cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8);
  ChildCommand(cmd,children[terms[i]]);
if 0 in NN then
cmd:=Concatenation("B",String(terms[i]),":=Length(Homology(Y",String(terms[i]),",0));");
ChildCommand(cmd,children[terms[i]]);
fi;
od;

maps:=[];
B:=[];
for i in terms{[1..Length(terms)-1]} do
  maps[i]:=ChildGet(Concatenation("ff",String(i)),children[i]);
  if 0 in NN then
  B[i]:=ChildGet(Concatenation("B",String(i)),children[i]);
  fi;
od;

T:=[];
for i in terms do
  T[i]:=i;
od;

for i in [1..Length(T)] do
  if not IsBound(T[i]) then T[i]:=T[i-1]; B[i]:=B[i-1];fi;
od;

TT:=[];
lst:=1;
for i in [1..Length(T)-1] do
  if not T[i]=T[i+1] then TT[i]:=maps[lst]; lst:=i+1;
else
  TT[i]:=List(maps[lst],x->IdentityMapping(Source(x)));
  fi;
od;

P:=[];;
for i in N do
  P[i+1]:=LinearHomomorphismsPersistenceMat(List(TT,t->t[i+1]));
od;
if 0 in NN then
P[1]:=NullMat(Length(B),Length(B));
for i in [1..Length(B)] do
for j in [i..Length(B)] do
  P[1][i][j]:=B[j];
od;od;
fi;
return P;
end);
######################################################################
######################################################################

