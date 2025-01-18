#(C) Graham Ellis 2005-2006

#####################################################################
InstallGlobalFunction(ResolutionExtension,
function(arg)
local
	EEhomGG, RN, RG, TestFinite,PreImRep,
	N,E,G,
	NhomE,
	EhomG,
	GmapE,
	NEhomN,
	NEhomNrecord,
	EltsE,
	MultE,
	InvE,
	PreimagesRecordG,PreimagesRecordE,
	NisFinite,GisFinite,EisFinite,
	Lngth,T,
	AppendToElts,
	gn,i,j,x,y;

T:=0;

EEhomGG:=arg[1];
RN:=arg[2];
RG:=arg[3];
TestFinite:=false;
if Length(arg)>3 then
if arg[4]="TestFiniteness" then TestFinite:=true; fi;
fi;
if Length(arg)>4 then 
PreImRep:=arg[5];
else
	#############################################################
	PreImRep:=function(x);
	return PreImagesRepresentative(EEhomGG,x);
	end;
	#############################################################
fi;

N:=RN!.group;
E:=Source(EEhomGG);
#G:=Image(EEhomGG);
G:=Range(EEhomGG);

NisFinite:=false;
GisFinite:=false;
EisFinite:=false; 
if TestFinite then
if IsFinite(N) then
	if Order(N)<=Length(RN!.elts) then NisFinite:=true; fi;
fi;

if IsFinite(G) then
        if Order(G)<=Length(RG!.elts) then GisFinite:=true; fi;
fi;

EisFinite:=IsFinite(E);
fi;

if EisFinite then EltsE:=Elements(E);
else
EltsE:=[Identity(E)];
for gn in GeneratorsOfGroup(E) do
Append(EltsE,[gn,gn^-1]);
od;
fi;

EltsE:=SSortedList(EltsE);
gn:=Position(EltsE,Identity(E));
EltsE[gn]:=EltsE[1];
EltsE[1]:=Identity(E);


	########################################################
	AppendToElts:=function(x);
	Append(EltsE,[x]);
	end;
	########################################################

if GisFinite then
	#########################################
	EhomG:=function(x);
	return Position(RG!.elts,ImageElm(EEhomGG,EltsE[x]));
	end;
	#########################################
else
        #########################################
        EhomG:=function(x)
	local g,Eltg;
	Eltg:=ImageElm(EEhomGG,EltsE[x]); 
        g:=Position(RG!.elts,Eltg); 
	if g=fail then 
	RG!.appendToElts(Eltg); 
	Append(RG!.elts,[Eltg]);
	g:=Length(RG!.elts); fi;
	#if Position(RG!.elts,Eltg^-1)=fail then 
	#Append(RG!.elts,[Eltg^-1]);fi;
	return g;
        end;
        #########################################
fi;

if EisFinite then
	#########################################
	NhomE:=function(x);
	return Position(EltsE,RN!.elts[x]);
	end;
	#########################################
else	
        #########################################
        NhomE:=function(x)
	local e,Elte;
	Elte:=RN!.elts[x];
        e:=Position(EltsE,Elte);
	if e=fail then AppendToElts(Elte);
	e:=Length(EltsE); fi;
	#if Position(EltsE,Elte^-1)=fail then
	#Append(EltsE,[Elte^-1]); fi;
	return e;
        end;
        #########################################
fi;

if GisFinite and EisFinite then
PreimagesRecordE:=List([1..Order(G)],x->
	Position(EltsE,PreImRep(RG!.elts[x])));
	
	#########################################
	GmapE:=function(x);
	return PreimagesRecordE[x];
	end;
	#########################################
else	
PreimagesRecordG:=[];
PreimagesRecordE:=[];
	#########################################
	GmapE:=function(x)
	local e,Elte,Eltg,pos;
	Eltg:=RG!.elts[x];
	pos:=Position(PreimagesRecordG,Eltg); 
	if not pos=fail then 
	return PreimagesRecordE[pos]; fi;
	
	Elte:=PreImRep(Eltg);
	e:=Position(EltsE,Elte);
	if e=fail then AppendToElts(Elte);
        e:=Length(EltsE); fi;
	if Position(EltsE,Elte^-1)=fail then
        AppendToElts(Elte^-1); fi;
	Append(PreimagesRecordG,[Eltg]);
	Append(PreimagesRecordE,[e]);
        return e;
	end;
	#########################################
fi;

if NisFinite then
	#########################################
	NEhomN:=function(x);
	return Position(RN!.elts,EltsE[x]);
	end;
	#########################################
else
	#########################################
	NEhomN:=function(x)
	local p,Eltp;
	Eltp:= EltsE[x];
	p:=Position(RN!.elts,Eltp);
	if p=fail then RN!.appendToElts(Eltp);
	Append(RN!.elts,[Eltp]);
	p:=Length(RN!.elts); fi;
	return p;
	end;
	#########################################
fi;

if EisFinite then
	#########################################
	MultE:=function(x,y);
	return Position(EltsE,EltsE[x]*EltsE[y]);
	end;
	#########################################
else
	#########################################
	MultE:=function(x,y)
	local p,Eltp;
	Eltp:=EltsE[x]*EltsE[y];
	p:= Position(EltsE,Eltp);
	if p=fail then AppendToElts(Eltp);
	p:=Length(EltsE); fi;
	#if Position(EltsE,Eltp^-1)=fail then
        #Append(EltsE,[Eltp^-1]); fi;
	return p;
	end;
	#########################################
fi;

if EisFinite then
	#########################################
	InvE:=function(x);
	return Position(EltsE,EltsE[x]^-1);
	end;
	#########################################
else
	#########################################
	InvE:=function(x)
	local p;
	p:=(Position(EltsE,EltsE[x]^-1));
	if p=fail then AppendToElts(EltsE[x]^-1); 
	p:=Length(EltsE);fi;
	return p;
	end;
	#########################################
fi;

if (not EisFinite ) and (not Length(RN!.elts)=infinity) and HAPconstant<50 then

for x in RN!.elts do
for y in RG!.elts do
AppendToElts(x*PreImRep(y));
od;
od;

fi;

#Print("\n",[NisFinite, GisFinite, EisFinite],"\n");


T:=TwistedTensorProduct(RG,RN,EhomG,GmapE,NhomE,NEhomN,EltsE,MultE,InvE,E);

        ########################################################
        AppendToElts:=function(x);
        Append(T!.elts,[x]); 
        end;
        ########################################################


T!.appendToElts:=AppendToElts;

return T;
end);
#####################################################################


