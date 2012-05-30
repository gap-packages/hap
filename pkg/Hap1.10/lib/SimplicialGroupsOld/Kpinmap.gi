InstallGlobalFunction(EilenbergMacLaneSimplicialGroupMap,
function(f,n,len)
	local 
		H,G,KH,KG,
		Mapf,Mapping,HH,GG,
		GensH,Prj,Emb,ImgsH,
		i,j,k,t,lengens,
		h,g;
		
H:=f!.Source;
G:=f!.Range;
KH:=EilenbergMacLaneSimplicialGroup(H,n,len);
KG:=EilenbergMacLaneSimplicialGroup(G,n,len);
Mapf:=[];
for i in [0..n-2] do
	Mapf[i+1]:=GroupHomomorphismByImages(Group(Identity(H)),Group(Identity(G)),[],[]);
od;
Mapf[n]:=f;  ##n-1
for i in [n..len] do
	HH:=KH!.groupsList(i);
	GG:=KG!.groupsList(i);
	GensH:=GeneratorsOfGroup(HH);
	Prj:=[];
	Emb:=[];
	k:=Length(HH!.DirectProductInfo!.groups);	
	for j in [1..k] do
		Prj[j]:=Projection(HH,j);
		Emb[j]:=Embedding(GG,j);
	od;
	ImgsH:=[];
	lengens:=Length(GensH);
    for j in [1..lengens] do
		h:=GensH[j];
		g:=Identity(GG);
		for t in [1..k] do
			g:=g*Image(Emb[t],Image(f,Image(Prj[t],h)));
		od;
		ImgsH[j]:=g;
	od;
	Mapf[i+1]:=GroupHomomorphismByImages(HH,GG,GensH,ImgsH);
od;
###################
Mapping:=function(i)
	return Mapf[i+1];
end;
###################
return Objectify(HapSimplicialGroupMap,
       rec(
            source:=KH,
            target:=KG,
            mapping:=Mapping,
            properties:=[["length",len]]
          ));
end);		
