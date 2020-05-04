################################################
################################################
InstallGlobalFunction(ResolutionPSL2QuadraticIntegers,
function(d,n)
local K, PK, R, S, D, name, ints,x, Rgroup, Kgroup;

ints:=[ -1, -2, -3, -5, -6, -7, -10, -11, -13, -14, -15, -17, -19, -21, -22, -23, -26, -43, -67, -163, "-26+I", "-22+I", "-21+I2", "-21+I3", "-21+I4", "-17+I", "-15+I", "-14+I", "-13+I", "-10+I", "-6+I", "-5+I"  ];
if not d in ints then
Print("PSL(2,Sqrt(d)) is implemented for d= -1, -2, -3, -5, -6, -7, -10, -11, -13, -14, -15, -17, -19, -21, -22, -23, -26, -43, -67, -163.\n");
Print("PSL(O-d) is implemented for d= \"-26+I\", \"-22+I\", \"-21+I2\", \"-21+I3\", \"-21+I4\", \"-17+I\", \"-15+I\", \"-14+I\", \"-13+I\", \"-10+I\", \"-6+I\", \"-5+I\" .\n");
return fail;
fi;

if IsInt(d) then
name:=Concatenation("SL(2,O", String(d), ")");
else
name:=Concatenation("SLO",d,")");
fi;
if name="SL(2,O-2)" then 
K:=ContractibleGcomplex("SL2O-2_a");
Kgroup:=K!.group; Kgroup!.bianchiInteger:=-2;
else
K:=ContractibleGcomplex(name);
Kgroup:=K!.group;
fi;
D:=Group( -Identity(K!.group) );;
PK:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(PK,n);
Rgroup:=R!.group;
if IsInt(d) then
Rgroup!.bianchiInteger:=d;
else
Rgroup!.bianchiInteger:=fail;
fi;
R!.group:=Rgroup;

if not '(' in name then
name:=SplitString(name,['O']);
name:=Concatenation(name[1],"(2,O",name[2],")");
fi;
name:=Concatenation("P",name);
SetName(R!.group,name);

return R;
end);
################################################
################################################

################################################
################################################
InstallGlobalFunction(ResolutionSL2QuadraticIntegers,
function(arg)
local d,n,K, PK, R, S, D, name, ints,x, gens, Q,OQ,I,G,i,k;

d:=arg[1];
n:=arg[2];

ints:=[ -1, -2, -3, -5, -6, -7, -10, -11, -13, -14, -15, -17, -19, -21, -22, -23, -26, -43, -67, -163, "-26+I", "-22+I", "-21+I2", "-21+I3", "-21+I4", "-17+I", "-15+I", "-14+I", "-13+I", "-10+I", "-6+I", "-5+I"  ];
if not d in ints then
Print("PSL(2,Sqrt(d)) is implemented for d= -1, -2, -3, -5, -6, -7, -10, -11, -13, -14, -15, -17, -19, -21, -22, -23, -26, -43, -67, -163.\n");
Print("PSL(O-d) is implemented for d= \"-26+I\", \"-22+I\", \"-21+I2\", \"-21+I3\", \"-21+I4\", \"-17+I\", \"-15+I\", \"-14+I\", \"-13+I\", \"-10+I\", \"-6+I\", \"-5+I\" .\n");
return fail;
fi;

if IsInt(d) then
name:=Concatenation("SL(2,O", String(d), ")");
else
name:=Concatenation("SLO",d,")");
fi;
if name="SL(2,O-2)" then 
K:=ContractibleGcomplex("SL(2,O-2)_a"); ;
else
K:=ContractibleGcomplex(name);
fi;
#D:=Group( -Identity(K!.group) );;
#K:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(K,n);

if not '(' in name then
name:=SplitString(name,['O']);
name:=Concatenation(name[1],"(2,O",name[2],")");
fi;
#name:=Concatenation("P",name);
SetName(R!.group,name);

if Length(arg)=3 then
  if arg[3]=true then
  for k in [1..n] do
    for i in [1..R!.dimension(k)] do
    R!.boundary(k,i);
    od; 
  od;
  Q:=QuadraticNumberField(d);;OQ:=RingOfIntegers(Q);;I:=QuadraticIdeal(OQ,1);;
  Apply(R!.elts,x->HAP_4x4MatTo2x2Mat(x,d));
  G:=HAP_CongruenceSubgroupGamma0(I);;
  G!.tree:=true;
  R!.group:=G;
  fi;
fi;
return R;
end);
################################################
################################################


################################################
################################################
InstallGlobalFunction(ResolutionGL2QuadraticIntegers,
function(d,n)
local K, PK, R, S, D, name, ints,x;

ints:=[ -1, -10, -11, 11, -13, 13, -14, 14, -15, -17, -19, -2, -21, -22, -23, -26, 2, -3, 3, -43, -5, 5, -6, 6, -7, 7, "GL2O10Steinitz1HAP",
"GL2O10Steinitz2HAP", "GL2O15Steinitz1HAP", "GL2O15Steinitz2HAP"  ];
if not d in ints then
Print("GL(2,Sqrt(d)) is implemented for d= -43, -26, -23, -22, -21, -19, -17, -15, -14, -13, -11, -10, -7, -6, -5, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 14, \"GL2O10Steinitz1HAP\", \"GL2O10Steinitz2HAP\", \"GL2O15Steinitz1HAP\", \"GL2O15Steinitz2HAP\" .\n");
return fail;
fi;

if IsInt(d) then
   if d<0 then
   name:=Concatenation("GL(2,O", String(d), ")");
   else
   name:=Concatenation("GL(2,O", String(d), "HAP)");
   fi;
else
name:=d;
fi;
K:=ContractibleGcomplex(name);
#D:=Group( -Identity(K!.group) );;
#K:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(K,n);

if IsInt(d) then
  if d>0 then
  name:=name{[1..Length(name)-4]}; name:=Concatenation(name,")");
  fi;
fi;
#name:=Concatenation("P",name);
SetName(R!.group,name);

return R;
end);
################################################
################################################


################################################
################################################
InstallGlobalFunction(ResolutionPGL2QuadraticIntegers,
function(d,n)
local K, PK, R, S, D, name, ints,x;

ints:=[ -1, -10, -11, 11, -13, 13, -14, 14, -15, -17, -19, -2, -21, -22, -23, -26, 2, -3, 3, -43, -5, 5, -6, 6, -7, 7, "GL2O10Steinitz1HAP",
"GL2O10Steinitz2HAP", "GL2O15Steinitz1HAP", "GL2O15Steinitz2HAP"  ];
if not d in ints then
Print("PGL(2,Sqrt(d)) is implemented for d= -43, -26, -23, -22, -21, -19, -17, -15, -14, -13, -11, -10, -7, -6, -5, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 14, \"GL2O10Steinitz1HAP\", \"GL2O10Steinitz2HAP\", \"GL2O15Steinitz1HAP\", \"GL2O15Steinitz2HAP\" .\n");
return fail;
fi;

if IsInt(d) then
   if d<0 then
   name:=Concatenation("GL(2,O", String(d), ")");
   else
   name:=Concatenation("GL(2,O", String(d), "HAP)");
   fi;
else
name:=d;
fi;
K:=ContractibleGcomplex(name);
D:=Group( -Identity(K!.group) );;
PK:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(PK,n);

if IsInt(d) then
  if d>0 then
  name:=name{[1..Length(name)-4]}; name:=Concatenation(name,")");
  fi;
fi;
name:=Concatenation("P",name);
SetName(R!.group,name);

return R;
end);
################################################
################################################

################################################
################################################
InstallGlobalFunction(ResolutionGL3QuadraticIntegers,
function(d,n)
local K, PK, R, S, D, name, ints,x;

ints:=[ -1, -2, -3, -5, -7, -11, -15 ];
if not d in ints then
Print("GL(3,Sqrt(d)) is implemented for d=  -1, -2, -3, -5, -7, -11, -15  .\n");
return fail;
fi;

if IsInt(d) then
   name:=Concatenation("GL(3,O", String(d), "HAP)");
else
name:=d;
fi;
K:=ContractibleGcomplex(name);
#D:=Group( -Identity(K!.group) );;
#K:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(K,n);

if IsInt(d) then
  name:=name{[1..Length(name)-4]}; name:=Concatenation(name,")");
fi;
#name:=Concatenation("P",name);
SetName(R!.group,name);

return R;
end);
################################################
################################################


################################################
################################################
InstallGlobalFunction(ResolutionPGL3QuadraticIntegers,
function(d,n)
local K, PK, R, S, D, name, ints,x;

ints:=[ -1, -2, -3, -5, -7, -11, -15 ];
if not d in ints then
Print("GL(3,Sqrt(d)) is implemented for d=  -1, -2, -3, -5, -7, -11, -15  .\n");
return fail;
fi;

if IsInt(d) then
   name:=Concatenation("GL(3,O", String(d), "HAP)");
else
name:=d;
fi;
K:=ContractibleGcomplex(name);
D:=Group( -Identity(K!.group) );;
PK:=QuotientOfContractibleGcomplex(K,D);;
R:=FreeGResolution(PK,n);

if IsInt(d) then
  name:=name{[1..Length(name)-4]}; name:=Concatenation(name,")");
fi;
name:=Concatenation("P",name);
SetName(R!.group,name);

return R;
end);
################################################
################################################



