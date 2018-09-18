InstallGlobalFunction(PolymakeFaceLattice,
function(arg)
local  polygon, callPolymake, F,D, L, i,j,s,t, toggle, DD,x,y;

polygon:=arg[1];
if Length(arg)=2 then toggle:=false;
else toggle:=true;
fi;

    callPolymake:=function(object,splitoption)
        local   returnedstring,  stdout,  stdin,  dir,  exitstatus;

        returnedstring:="";
        stdout:=OutputTextString(returnedstring,false);
        stdin:=InputTextNone();;
        dir:=DirectoryOfPolymakeObject(object);
        if dir=fail
           then
            dir:=DirectoryCurrent();
        fi;
        exitstatus:=Process( dir, POLYMAKE_COMMAND, stdin, stdout,
                            Concatenation([FullFilenameOfPolymakeObject(object)],
                                    splitoption)
                            );;
        CloseStream(stdout);
        CloseStream(stdin);
        return rec(status:=exitstatus,stdout:=stdout,string:=returnedstring);
    end;


if toggle then
F:=callPolymake(polygon,["HASSE_DIAGRAM->FACES"]);
else
F:=callPolymake(polygon,["HASSE_DIAGRAM->ADJACENCY"]);
fi;
F:=F.string;
F:=SplitString(F,"\n");
F:=F{[2..Length(F)]};
F:=ConvertPolymakeListOfSetsToGAP(F);
F:=F+1;

D:=callPolymake(polygon,["F_VECTOR"]);
D:=SplitString(D.string,"\n")[2];
D:=SplitString(D," ");
Apply(D,EvalString);
D:=Reversed(D);
DD:=[1];;
x:=1;;
for y in D do
x:=x+y;
Add(DD,x);
od;
D:=DD;

Add(D,Length(F)-1);

L:=List([1..Length(D)-1],i->[]);

s:=1;
for i in [1..Length(D)-1] do
t:=D[i+1]-1;
for j in [s+1..t+1] do
Add(L[i],F[j]);
od;
s:=t+1;
od;

if Length(L[1][1])=1 then

L:=Reversed(L);
Remove(L,1);
Add(L,[[]]);
fi;


return L;
end);
