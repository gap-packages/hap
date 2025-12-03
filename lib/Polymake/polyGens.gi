#(C) Graham Ellis, 2005-2006

#####################################################################
#####################################################################
InstallGlobalFunction(PolytopalGenerators,
function(GG,v)
local G, Points, vertices, CG, x, w, N, proj, poly, EdgeGenerators, tmp, index, i, Faces, FacesFinal;


if not (IsPermGroup(GG) or IsMatrixGroup(GG)) then 
    Print("The group G must be a permutation or matrix group.\n");
    return fail;
fi;

if IsPermGroup(GG) then
    GG:=Image(PermToMatrixGroup(GG));
fi;

#####################################################################

Points:=[];

G:=[];
for x in GG do
    w:=x*v;  # warning: left action!
    if not w in Points and not w=v then
        Add(Points,w);
        Add(G,x);
    fi;
od;
Points:=Set(Points);

###################### CALCULATE CENTRE OF GRAVITY ##################

CG:=Sum(Points)/Size(Points);

##################### CENTRE OF GRAVITY DONE ########################


	################# PROJECTION ################################
	N:=CG-v;    #This might be parallel to an edge!!! 
	proj:=function(w)
        local k, m  ;
        m:=(w-v)*N;
        k:= ((CG-v)*N)/m;
        return v+(k*(w-v));
	end;
	#############################################################
	
################## CALCULATE HULL OF PROJECTED POINTS ###############
Points := List(Points, proj);

poly:=CreatePolymakeObject();
AppendPointlistToPolymakeObject(poly,Points);

################# HULL CALCULATED ###################################

################# READ VERTICES #####################################
vertices := Polymake(poly,"VERTICES");
################ VERTICES READ ######################################

################ RECOVER THE EDGE GENERATORS ########################
EdgeGenerators:=[];
for w in Points do
    if w in vertices then
        x:=G[Position(Points,w)];
        EdgeGenerators[Position(vertices,w)]:=x;
    fi;
od;
################ EDGE GENERATORS RECOVERED ##########################

################ READ HASSE DIAGRAM #################################

tmp := Polymake(poly,"F_VECTOR");
index:=[1];                   #because Polymake has
for i in [1..Length(tmp)] do  #discontinued the DIMS
Add(index,index[i]+tmp[i]);   #property.
od;                           #

Faces := Polymake(poly,"FACES");
if Length(Faces[1])=0 then
    Remove(Faces,1);
fi;
Faces:=List(Faces, f -> f-1);

if Length(Faces[1])=1 then

	FacesFinal:=[];
	for i in [1..Length(index)-1] do
    	Append(FacesFinal,[[index[i]..index[i+1]-1]]);
	od;
	Append(FacesFinal,[[index[i+1]]]);
	FacesFinal:=(List(FacesFinal,x->List(x,i->Faces[i])));

else

	FacesFinal:=[[1..index[1]]];
	for i in [1..Length(index)-1] do
	    Append(FacesFinal,[[index[i]+1..index[i+1]]]);
	od;
	FacesFinal:=Reversed(List(FacesFinal,x->List(x,i->Faces[i])));
	
fi;

############### HASSE DIAGRAM READ ##################################
return rec(
             generators:=EdgeGenerators,
             hasseDiagram:=FacesFinal,
	     vector:=v);
	     
end);
#####################################################################
