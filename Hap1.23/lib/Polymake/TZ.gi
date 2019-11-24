#DeclareGlobalFunction("HAP_TzPair");
#################################################################
#################################################################
InstallGlobalFunction( HAP_TzPair, function ( arg )

    local geni, genj, gens, k, m, n, num, pairs, T, tietze;

    # check the first argument to be a Tietze record.
    T := arg[1];
    TzCheckRecord( T );

    # get the second argument.
    n := 1000000;

    # intialize the local variables.
    tietze := T!.tietze;
    gens := tietze[TZ_GENERATORS];

    # determine the most frequently occurring pair.
    pairs := TzMostFrequentPairs( T, n );
    pairs:=Filtered(pairs,x->x[1]>1);
    pairs:=SSortedList(pairs);
if Length(pairs)=0 then return fail; fi;
    # output it.
        num := pairs[1][1];
        k := pairs[1][4];
        geni := gens[pairs[1][2]];
        if k > 1 then  geni := geni^-1;  fi;
        genj := gens[pairs[1][3]];
        if k mod 2 = 1 then  genj := genj^-1;  fi;
        
        return   [num,  geni*genj];
     
end );
#################################################################
#################################################################

#DeclareGlobalFunction("HAP_AddGenerator");
#################################################################
#################################################################
InstallGlobalFunction(HAP_AddGenerator,
function(G)
local FG,F,P,g,w,rels,newrels,sub,by,rr,r,s,ns,nt,bool,mx  ;
P:=PresentationFpGroup(G);
w:=HAP_TzPair(P);
if w=fail then return G; fi;

if w[1]=1 then return G; fi;
w:=w[2];
AddGenerator(P);
g:=GeneratorsOfPresentation(P);
g:=g[Length(g)];
rr:=g^-1*w;
#AddRelator(P,r);
FG:=FpGroupPresentation(P);
rels:=RelatorsOfFpGroup(FG);
F:=FreeGroupOfFpGroup(FG);
w:=MappedWord(w,GeneratorsOfPresentation(P),GeneratorsOfGroup(F));
g:=MappedWord(g,GeneratorsOfPresentation(P),GeneratorsOfGroup(F));
rr:=MappedWord(rr,GeneratorsOfPresentation(P),GeneratorsOfGroup(F));
sub:=w; by:=g;
newrels:=[];
for r in rels do
s:=r;
bool:=true;
while bool do
ns:=SubstitutedWord(s,sub,1,by); 
if not ns=fail then s:=ns; else bool:=false; fi;
ns:=SubstitutedWord(s,sub^-1,1,by^-1);
if not ns=fail then s:=ns; bool:=true; fi;
od;
Add(newrels,s);
od;
Add(newrels,rr);
return F/newrels;
return [newrels,sub,by];
end);
#################################################################
#################################################################

#DeclareGlobalFunction("SmoothedFpGroup");
#################################################################
#################################################################
InstallGlobalFunction(SmoothedFpGroup,
function(G)
local mx, m, bool, GG, H;
GG:=G;
mx:=Reversed(SortedList((List(RelatorsOfFpGroup(GG),Length))));

bool:=true;
while bool do
H:=HAP_AddGenerator(GG);
m:=Reversed(SortedList((List(RelatorsOfFpGroup(H),Length))));
if m<mx then GG:=H; mx:=m;
else bool:=false; fi; 

od;
return GG;
end);
#################################################################
#################################################################

