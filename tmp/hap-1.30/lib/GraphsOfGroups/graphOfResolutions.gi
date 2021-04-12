#(C) 2012 Graham Ellis

##############################################
InstallGlobalFunction(GraphOfResolutions,
function(D,N)
local v,e,r,s,t,u,R,i,j,names;

R:=[];
names:=[];

for v in D do
if IsGroup(v) then
r:=ResolutionGenericGroup(v,N);
r!.group:=v;
Add(R,r);
Add(names,Name(v));
fi;
od;

for e in D do
if IsList(e) then
u:=[];
r:=ResolutionGenericGroup(Source(e[1]),N);
r!.group:=Source(e[1]);
t:=Range(e[1]);
j:=Position(names,Name(t));
Add(u,EquivariantChainMap(r,R[j],e[1]));
r:=ResolutionGenericGroup(Source(e[2]),N);
t:=Range(e[2]);
j:=Position(names,Name(t));
Add(u,EquivariantChainMap(r,R[j],e[2]));
Add(R,u);
fi;
od;

return R;
end);
##############################################

##############################################
InstallGlobalFunction(GraphOfResolutionsToGroups,
function(D)
local G,v,e;

G:=[];

for v in D do
if IsHapResolution(v) then
Add(G,v!.group);
fi;
if IsList(v) then
e:=[v[1]!.originalHom,v[2]!.originalHom];
Add(G,e);
fi;
od;

return G;
end);
##############################################

##############################################
InstallGlobalFunction(GraphOfResolutionsTest,
function(D)
local v;

if not IsList(D) then return false; fi;

for v in D do
if not 
(IsHapResolution(v) or IsList(v))
then return false;
fi;
od;

return GraphOfGroupsTest(GraphOfResolutionsToGroups(D));

end);
##############################################

##############################################
InstallGlobalFunction(GraphOfResolutionsDisplay,
function(D);

if not GraphOfResolutionsTest(D) then 
Print("Input is not a graph of resolutions.\n");
return fail;
fi;

GraphOfGroupsDisplay(GraphOfResolutionsToGroups(D));

end);
###############################################

##############################################
InstallGlobalFunction(TreeOfGroupsToContractibleGcomplex,
function(D,G)
local Vertices, VerticesNames, Edges, Dimension, 
PseudoBoundary, Boundary, Elts, Action, Stabilizer, ID;

Vertices:=Filtered(D,x->IsGroup(x));
VerticesNames:=List(Vertices,x->Name(x));
Edges:=Filtered(D,x->IsList(x));
ID:=Group(One(G));
Elts:=[One(G)];

#######################
Dimension:=function(n);
if not n in [0,1] then return 0; fi;
if n=0 then return Length(Vertices); fi;
if n=1 then return Length(Edges); fi;
end;
#######################

PseudoBoundary:=[];
#######################
Boundary:=function(n,k);
if not n=1 then return []; fi;
if not IsBound(PseudoBoundary[k]) then
PseudoBoundary[k]:=
[[Position(VerticesNames,Name(Range(Edges[k][1]))),1],
[-Position(VerticesNames,Name(Range(Edges[k][2]))),1]];
fi;
if k>0 then 
return PseudoBoundary[k];
else
return NegateWord(PseudoBoundary[k]);
fi;
end;
#######################

#######################
Action:=function(n,k,l);
return 1;
end;
#######################

#######################
Stabilizer:=function(n,k);
if n=0 then return 
Vertices[k]; fi;
if n=1 then return
Source(Edges[k][1]); fi;
return ID;
end;
#######################


return Objectify(HapNonFreeResolution,
            rec(
            dimension:=Dimension,
            boundary:=Boundary,
            homotopy:=fail,
            elts:=Elts,
            group:=G,
            stabilizer:=Stabilizer,
            action:=Action,
            properties:=
            [["length",10000],
             ["characteristic",0],
             ["type","resolution"]]  ));

end);
###############################################

###############################################
InstallGlobalFunction(TreeOfResolutionsToContractibleGcomplex,
function(R,G)
local D, C, v,s, NamesOfGroups, Resolutions;

D:=GraphOfResolutionsToGroups(R);
C:=TreeOfGroupsToContractibleGcomplex(D,G);

NamesOfGroups:=[];
Resolutions:=[];
for v in R do
if IsHapResolution(v) then
Add(Resolutions,v);
Add(NamesOfGroups,Name(v!.group));
else
s:=Source(v[1]);
Add(Resolutions,s);
Add(NamesOfGroups,Name(s!.group));
fi;
od;

C!.resolutions:=[Resolutions,NamesOfGroups];
return C;
end);
###############################################

############################################
InstallGlobalFunction(ConjugatedResolution,
function(R,p)
local newelts, eltlst, update,G,pinv;

pinv:=p^-1;
newelts:=[];
eltlst:=[];
newelts:=ListToPseudoList(newelts);
SetIsPseudoListWithFunction(newelts,true);
newelts!.eltfun:= function(i) local g, el;
              if not IsBound(eltlst[i]) then
                  el:=p*R!.elts[i]*pinv;
                  eltlst[i]:=el;
              fi;
              return eltlst[i];
              end;
newelts!.addfun:=function(g);
                 Add(R!.elts,pinv*g*p);
                 end;
newelts!.lnthfun:=function(); return Length(R!.elts); end;
newelts!.posfun:=function(g) local pos; 
		 pos:= Position(R!.elts,pinv*g*p); 
		 return pos; 
                 end;

G:=R!.group^(pinv);
SetName(G,String(Random([1..1000000])) );
return Objectify(HapResolution,
                rec(
                dimension:=R!.dimension,
                boundary:=R!.boundary,
                homotopy:=R!.homotopy,
                elts:=newelts,
                group:=G,
                properties:=R!.properties));

end);
############################################
