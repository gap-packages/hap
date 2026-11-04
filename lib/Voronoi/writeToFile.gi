#This code was written by Sebastian Schoennenbeck

#####################################################################
#####################################################################
WriteContractibleGComplexToFile:=function(R,file)
local output, length, x, i, j, Stab, StabList, Rot, ListIFace, ListSign, 
ListElt, xx;

output:=OutputTextFile(file,true);

WriteAll(output, "HAP_GCOMPLEX_SETUP:=[false]; \n");

WriteAll(output, "HAP_GCOMPLEX_LIST:=[");

length:=0;
x:=false;

while not x do 
	if R!.dimension(length)=0 then 
		x:=true;
		break;
	fi;
	length:=length+1;
od;

length:=length-1;

for i in [0..length] do
	WriteAll(output, "[");
	for j in [1..R!.dimension(i)] do
		Stab:=R!.stabilizer(i,j);
		StabList:=List(Elements(Stab),x->x);
		Rot:=[];
		for x in StabList do 
			xx:=Position(R!.elts,x);
			if R!.action(i,j,xx)=1 then
				Add(Rot,x);
			fi;
		od;
		WriteAll(output,"rec(");
		WriteAll(output,"TheMatrixStab:=Group(");
		PrintTo(output, StabList);
		WriteAll(output,"),");
		WriteAll(output,"TheRotSubgroup:=Group(");
		PrintTo(output,Rot);
		WriteAll(output,"),");
		ListIFace:=List(R!.boundary(i,j),x->AbsInt(x[1]));
		ListSign:=List(R!.boundary(i,j),x->SignInt(x[1]));
		ListElt:=List(R!.boundary(i,j),x->R!.elts[x[2]]);
		WriteAll(output,"BoundaryImage:=rec(");
		WriteAll(output,"ListIFace:=");
		PrintTo(output,ListIFace);
		WriteAll(output,",ListSign:=");
		PrintTo(output,ListSign);
		WriteAll(output,",ListElt:=");
		PrintTo(output,ListElt);
		WriteAll(output,"))");
		if j <> R!.dimension(i) then
			WriteAll(output,",");
		fi;
	od;
	WriteAll(output,"]");
	if i <> length then
		WriteAll(output,",");
	fi;
od;
WriteAll(output,"];");
CloseStream(output);
end;
#####################################################################
#####################################################################
