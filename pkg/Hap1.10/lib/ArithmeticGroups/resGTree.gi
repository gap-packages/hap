InstallGlobalFunction(ResolutionGTree,
function(arg)
local 
	R,n,EltsG,G,i,L,
	StabRes,StabGrps,Triple2Pair,Quad2One,GrpsRes,Pair2Quad,Quad2Pair,
	Dimension,Hmap,pos,Mult,AlgRed,CorrectList,Boundary,
	Homotopy,FinalHomotopy,
	StRes,hmap,Action,PseudoBoundary,PseudoHomotopy,ZeroDimensionHmap,ZeroDimensionHtpy,
	HmapRec,p,q,r,s,HtpyRec,k,g;
R:=arg[1];
n:=arg[2];
EltsG:=R!.elts;
G:=R!.group;
Action:=function(p,r,g)
if not IsBound(R!.action) then return 1;
else return R!.action(p,r,g);
fi;
end;
#############################################
pos:=function(g)   #return the position of gth - element in EltsG, 
		    # if not then add to EltsG
if Position(R!.elts,g)=fail then 
	Add(R!.elts,g); 
fi;
return Position(R!.elts,g);
end;
#############################################
AlgRed:=function(g) #Algebraic Reduction
local l,x;
l:=[];
for x in g do
	if Position(l,[-x[1],x[2]])=fail then Add(l,x); 
	else 
		Remove(l,Position(l,[-x[1],x[2]]));
	fi;
od;
return l;
end;
#############################################
Mult:=function(g,w) # Multiply gth-element with a word
local 
	l,x;
l:=[];
if R!.elts[g]=[] then return [];fi;
Append(l,List(w,y->[y[1],pos(R!.elts[g]*R!.elts[y[2]])]));
return AlgRed(l);
end;
##############################################
GrpsRes:=function(G,n) # Resolutions of Group
local 
	iso,Q,res,x;
if IsBound(R!.resolutions) and HasName(G) then 
x:=Position(R!.resolutions[2], Name(G)); 
if not x=fail then return R!.resolutions[1][x]; fi;
fi;
iso:=RegularActionHomomorphism(G);
Q:=Image(iso);
res:=ResolutionFiniteGroup(Q,n);
res!.group:=G;
res!.elts:=List(res!.elts,x->PreImagesRepresentative(iso,x));

return res;
end;
#############################################
#Create list of stabilizer groups and resolutions
StabGrps:= List([0..Length(R)],n->
           List([1..R!.dimension(n)], k->R!.stabilizer(n,k))); 
StabRes:=[];
for L in StabGrps do
Add(StabRes,List(L,	
g->ExtendScalars(GrpsRes(g,n),G,EltsG))
);

od;

#############################################
CorrectList:=function(list)
local 
	l,i;
if list=[] then return [];fi;
l:=StructuralCopy(list[1]);
for i in [2..Length(list)] do
		Append(l,StructuralCopy(list[i]));
	od;
return l;
end;
############################################
# return nth-generator of F_(p,q) from (r,s)th-generator of 
# stabilizer 
Quad2One:=function(p,q,r,s)
local
	n,d,i,j;
n:=0;
i:=SignInt(s);
s:=AbsInt(s);
d:=List([1..R!.dimension(p)],x->StabRes[p+1][x]!.dimension(q));
for j in [1..r-1] do
	n:=n+d[j];
od;
n:=n+s;
if q=0 and n>R!.dimension(p) then n:=R!.dimension(p);fi;
return i*n;
end;
############################################
Triple2Pair:=function(p,q,n)
local 
	r,s,d,i;
r:=0;
d:=List([1..R!.dimension(p)],x->StabRes[p+1][x]!.dimension(q));
i:=SignInt(n);
n:=AbsInt(n);
while n>0 do
	r:=r+1;
	s:=n;
	n:=n-d[r];
od;
return [r,i*s];
end;
############################################
HmapRec:=[];
for p in [1..2] do
  HmapRec[p]:=[];
  for q in [1..n+1] do
    HmapRec[p][q]:=[];
    for r in [1..R!.dimension(p-1)] do
      HmapRec[p][q][r]:=[];
    od;
  od;
od;
############################################
ZeroDimensionHmap:=function(k)
local i,j,pk;
pk:=AbsInt(k);
j:=0;
for i in [1..pk-1] do
j:=j+StabRes[1][i]!.dimension(0);
od;
j:=j+1;
if k>0 then return j;
else return -j;fi;
end;
############################################
Hmap:=function(p,q,r,s)     #Horiziontal map Hmap:A(p,q)->A(p-1,q), acts on the (r,s) th-generator of A(p,q)
local 
	i,l,d0,m,bdr,ps,d1d0,w;
ps:=AbsInt(s);
if p<>1 then return [];
else
if not IsBound(HmapRec[p+1][q+1][r][ps]) then
	#if q=0 then bdr:=R!.boundary(1,Quad2One(p,q,r,ps)); 
	if q=0 then bdr:=StructuralCopy(R!.boundary(1,1));
		#Print("bdr",bdr);
		Apply(bdr,w->[ZeroDimensionHmap(w[1]),w[2]]);
		#Print("bdr",bdr);
		HmapRec[p+1][q+1][r][ps]:=bdr;
	else
		l:=[];m:=[];
		d0:=StructuralCopy(List(StabRes[p+1][r]!.boundary(q,ps),x->[Action(p,r,x[2])*x[1],x[2]]));
		#d1d0:=AlgRed(CorrectList(List(d0,x->Mult(x[2],Hmap(p,q-1,r,x[1])))));
		for w in d0 do
			Append(m,Mult(w[2],Hmap(p,q-1,r,w[1])));
		od;
		Apply(m,x->[Triple2Pair(p-1,q-1,x[1]),x[2]]);
		for w in m do
			Append(l,List(StabRes[p][w[1][1]]!.homotopy(q-1,[w[1][2],w[2]]),y->[Quad2One(p-1,q,w[1][1],y[1]),y[2]]));
		od;
		HmapRec[p+1][q+1][r][ps]:=AlgRed(l);
	fi;
fi;
fi;
if SignInt(s)=1 then return HmapRec[p+1][q+1][r][ps];
else return NegateWord(HmapRec[p+1][q+1][r][ps]);fi;
end;
##############################################
Pair2Quad:=function(k,n)
local 
	p,q,r,s,i,temp,j1,j2;
temp:=0;
for j1 in [0..k] do
	for j2 in [1..R!.dimension(j1)] do
		temp:=temp+StabRes[j1+1][j2]!.dimension(k-j1);
	od;
od;
if n>temp then return "generator does not exist";fi;
p:=-1;
i:=SignInt(n);
n:=AbsInt(n);
while n>0 do;
	p:=p+1;
	r:=0;
	while (n>0 and r<R!.dimension(p)) do
		r:=r+1;
		s:=n;
		n:=n-StabRes[p+1][r]!.dimension(k-p);
	od;
od;
q:=k-p;
return [p,q,r,i*s];
end;
##############################################
Quad2Pair:=function(p,q,r,s)
local
	k,n,i,j;
k:=p+q;
n:=0;
for i in [0..p-1] do
	for j in [1..R!.dimension(i)] do
		n:=n+StabRes[i+1][j]!.dimension(k-i);
	od;
od;	
for i in [1..r-1] do
	n:=n+StabRes[p+1][i]!.dimension(k-p);
od;
n:=n+AbsInt(s);
return [k,SignInt(s)*n];
end;
#############################################
PseudoBoundary:=[];
for k in [1..n+1] do
    PseudoBoundary[k]:=[];
od;
##############################################
Boundary:=function(k,n)
local
	d,l,p,q,r,s,w,pn;
pn:=AbsInt(n);
if not IsBound(PseudoBoundary[k+1][pn]) then
w:=Pair2Quad(k,pn);
p:=w[1];q:=w[2];r:=w[3];s:=w[4];
d:=[];
if q<>0 then 
	l:=StructuralCopy(List(StabRes[p+1][r]!.boundary(q,s),x->[Quad2Pair(p,q-1,r,Action(p,r,x[2])*x[1])[2],x[2]]));
	Append(d,StructuralCopy(l));
fi;
if IsEvenInt(q) then 
	Append(d,StructuralCopy(Hmap(p,q,r,s)));
else
	Append(d,StructuralCopy(NegateWord(Hmap(p,q,r,s))));
fi;
PseudoBoundary[k+1][pn]:=AlgRed(d);	
fi;
if SignInt(n)=1 then 
    return PseudoBoundary[k+1][pn];
else return NegateWord(PseudoBoundary[k+1][pn]);
fi;
end;
##############################################
Dimension:=function(n)
local
	dim,p,i;
dim:=0;
for p in [0..n] do
	for i in [1..R!.dimension(p)] do
		dim:=dim+StabRes[p+1][i]!.dimension(n-p);
	od;
od;
return dim;
end;
##############################################
HtpyRec:=[];
for k in [1..n] do
  HtpyRec[k]:=[];
  for s in [1..Dimension(k-1)] do
    HtpyRec[k][s]:=[];
  od;
od;
##############################################
ZeroDimensionHtpy:=function(k)
local i,j,r;
i:=0;
#r:=i;
#k:=k-StabRes[1][i]!.dimension(0);
while k>0 do
  i:=i+1;
  k:=k-StabRes[1][i]!.dimension(0);
  r:=i;
od;
return r;
end;
##############################################
Homotopy:=function(n,w)
local 
	t,g,h0,h11,e,h,dh,
	p,q,r,s,v,m,pt,ppt,
	h1,d1h1,x,k,y,ps;
t:=w[1];
g:=w[2];
e:=[];
h:=[];
dh:=[];
pt:=AbsInt(t);
v:=Pair2Quad(n,pt);#Print(v);
p:=v[1];q:=v[2];r:=v[3];s:=v[4];
if not IsBound(HtpyRec[n+1][pt][g]) then
if n=0 then
	ppt:=ZeroDimensionHtpy(pt);
	h1:=StructuralCopy(R!.homotopy(n,[ppt,g]));
	#Print("h1=",h1,"\n");
	#Apply(h1,w->[Action(1,1,w[2])*w[1],w[2]]);
	#d1h1:=AlgRed(CorrectList(List(h1,x->Mult(x[2],Hmap(p+1,q,r,x[1])))));   #need to fix 'r'
	d1h1:=StructuralCopy(AlgRed(CorrectList(List(h1,x->Mult(x[2],Hmap(p+1,q,1,x[1]))))));
	for x in d1h1 do
	    k:=Pair2Quad(n,x[1]);
	    y:=StructuralCopy(StabRes[k[1]+1][k[3]]!.homotopy(q,[k[4],x[2]]));
	    Apply(y,w->[Quad2Pair(k[1],k[2]+1,k[3],w[1])[2],w[2]]);
	    Append(e,y);
	od;
	h0:=StructuralCopy(StabRes[p+1][r]!.homotopy(0,[s,g]));
	Apply(h0,w->[Quad2Pair(p,q+1,r,w[1])[2],w[2]]);
	h11:=List(h1,x->[Quad2Pair(p+1,q,Triple2Pair(p+1,q,x[1])[1],Triple2Pair(p+1,q,x[1])[2])[2],x[2]]);
	#Print("e=   ",e,"\n");
	#Print("h0=   ",h0,"\n");
	#Print("h11=   ",h11,"\n");
	Append(h,NegateWord(e));
	Append(h,h0);
	Append(h,h11);
	HtpyRec[n+1][pt][g]:=AlgRed(h);
else
	if p=0 then 
	  h0:=StructuralCopy(StabRes[p+1][r]!.homotopy(q,[s,g]));
	  Apply(h0,w->[Quad2Pair(p,q+1,r,w[1])[2],w[2]]);
	  Append(h,h0);
	else
	ps:=Action(1,1,g)*s;
	m:=StructuralCopy(StabRes[p+1][r]!.homotopy(q,[ps,g]));
	#Print("m=",m,"\n");
	Apply(m,x->[Action(1,1,x[2])*x[1],x[2]]);
	h0:=List(m,x->[Quad2Pair(p,q+1,r,x[1])[2],x[2]]);
	#Print("h=",h0,"\n");
	Append(h,h0);
	#dh:=AlgRed(CorrectList(List(m,x->Mult(x[2],Hmap(p,q+1,r,x[1])))));
	dh:=AlgRed(CorrectList(List(m,x->Mult(x[2],Hmap(p,q+1,1,x[1])))));
	for x in dh do
	   k:=Pair2Quad(n,x[1]);
	   y:=StructuralCopy(StabRes[k[1]+1][k[3]]!.homotopy(k[2],[k[4],x[2]]));
	   Apply(y,w->[Quad2Pair(k[1],k[2]+1,k[3],w[1])[2],w[2]]);
	   Append(e,y);
	   #k:=Pair2Quad(n,x[1]);
	   #Append(e,StabRes[k[1]+1][k[3]]!.homotopy(k[2],[k[4],x[2]]));
	od;
	#Print("e=",e,"\n");
	if IsEvenInt(q) then 
	  Append(h,e);
	else Append(h,NegateWord(e));fi;
	fi;
	HtpyRec[n+1][pt][g]:=AlgRed(h);
fi;
fi;
if SignInt(t)=1 then return HtpyRec[n+1][pt][g];
else return NegateWord(HtpyRec[n+1][pt][g]);
fi;
end;

##############################################
FinalHomotopy:=function(n,g)
if R!.homotopy=fail then
  return fail;
else return Homotopy(n,g);
fi;
end;
####ADDED MAY##############################################
StRes:=function(n,k)
return StabRes[n+1][k];
end;
##############################################
return Objectify(HapResolution,
                rec(
                dimension:=Dimension,
                boundary:=Boundary,
                homotopy:=FinalHomotopy,
                elts:=R!.elts,
                group:=R!.group,
                stabres:=StRes,
                properties:=
                   [["length",n],
                    ["initial_inclusion",true],
                    ["type","resolution"],
                    ["characteristic",EvaluateProperty(R,"characteristic")]  ]	));
end);