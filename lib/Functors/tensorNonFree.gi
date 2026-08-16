#####################################################################
#####################################################################
InstallGlobalFunction(TensorWithIntegersMod2Torsion,
function(R)
local
        Dimension,
        BoundaryRec,
        Boundary,
        LengthC,
        TFSummands,
        returnvec,
        bound,
        M, n, k, g, i, x, tmp, pos;

LengthC:=EvaluateProperty(R,"length");
M:=[1..LengthC];
TFSummands:=[];
BoundaryRec:=[];

TFSummands[1]:=[1..R!.dimension(0)];
for n in [1..LengthC] do
BoundaryRec[n]:=[];
TFSummands[n+1]:=[];
for k in [1..R!.dimension(n)] do
tmp:=true;
for g in Elements(R!.stabilizer(n,k)) do
pos := Position(R!.elts,g);
if pos=fail then Add(R!.elts,g); pos:=Length(R!.elts); fi;
if R!.action(n,k,pos)=-1 then tmp:=false; break; fi;
od;
if tmp then Add(TFSummands[n+1],k); fi;
od;
od;

##########################
Dimension:=function(n);
return Length(TFSummands[n+1]);
end;
##########################

for n in [1..LengthC] do
for k in TFSummands[n+1] do
#####################################################################

returnvec:=0*[1..R!.dimension(n-1)];
bound:=R!.boundary(n,k);
for x in [1..Size(bound)] do
i:=AbsInt(bound[x][1]);
returnvec[i]:=returnvec[i]+SignInt(bound[x][1]);
od;

returnvec:=returnvec{TFSummands[n]};
Add(BoundaryRec[n],returnvec);

#####################################################################
od;
od;


##########################
Boundary:=function(n,k);
if n<0 then return false; fi;
if n=0 then return [0]; fi;
return BoundaryRec[n][k];
end;
##########################

return Objectify(HapChainComplex,
                rec(
                TFSummands:=TFSummands,
                dimension:=Dimension,
                boundary:=Boundary,
                properties:=
                [["length",LengthC],
                ["connected",true],
                ["type", "chainComplex"],
                ["characteristic",
                EvaluateProperty(R,"characteristic")] ]));
end);
#####################################################################
#####################################################################

