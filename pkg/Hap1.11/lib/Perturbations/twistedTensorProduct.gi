#(C) Graham Ellis, 2005-2006

#####################################################################
InstallGlobalFunction(TwistedTensorProduct,
function(R,S,EhomG,GmapE,NhomE,NEhomN,EltsE,Mult,InvE)
local   
		DimensionR,BoundaryR,HomotopyR,
		DimensionS,BoundaryS,HomotopyS,
		Dimension,Boundary,Homotopy,
                FilteredLength, FilteredDimension, 
		DimPQ,DimPQrec,
		Int2Pair, Pair2Int,
		Htpy, HtpyRecord, CompHtpy,
		Del, CompDel, DelRecord,
                srtfn,
		PseudoBoundary,
		Charact,
		AddWrds, AddLst,
		i,j,k,l,n,p,q,r,rr,grp,dl, 
				######Remaining variables are concerned
				######with the contracting homotopy.
		Int2Vector,
		Vector2Int,
		HorizontalBoundaryGen,
		HorizontalBoundaryWord,
		HomotopyGradedGen,
		EmapN,
		HomotopyRec,
		Homtpy,
		HomotopyOfWord,
		FinalHomotopy,
		HorizontalPseudoBoundary,
		SMALL,SizeE,BoolE,MT,
#################################
AbsInt,                         #
SignInt;                        #
                                #
AbsInt:=AbsInt_HAP;             #
SignInt:=SignInt_HAP;           #
#################################
		
SMALL:=4096;;
if Order(S!.group) >SMALL  then SizeE:=Order(S!.group); fi;
if Order(R!.group) >SMALL
 then SizeE:=Order(R!.group); fi;
if not IsBound(SizeE) then
SizeE:=Order(R!.group)*Order(S!.group);
fi;

BoolE:=SizeE<=SMALL and Size(R!.group)=Size(R!.elts)
and Size(S!.group)=Size(S!.elts);
#################################
if BoolE then
MT:=[];
#for i in [1..SizeE] do
for i in [1..Length(EltsE)] do
MT[i]:=[];
#for j in [1..SizeE] do
for j in [1..Length(EltsE)] do
MT[i][j]:=Mult(i,j);
od;
od;
fi;
#################################

DimensionR:=R!.dimension;
DimensionS:=S!.dimension;
BoundaryR:=R!.boundary;
BoundaryS:=S!.boundary;
HomotopyR:=R!.homotopy;
HomotopyS:=S!.homotopy;
n:=Minimum(EvaluateProperty(R,"length"),EvaluateProperty(S,"length"));

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")=0
and EvaluateProperty(S,"characteristic")>0
then Charact:=EvaluateProperty(S,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")=0
then Charact:=EvaluateProperty(R,"characteristic");
fi;

if EvaluateProperty(R,"characteristic")>0
and EvaluateProperty(S,"characteristic")>0
then Charact:=Product(Intersection([
DivisorsInt(EvaluateProperty(R,"characteristic")),
DivisorsInt(EvaluateProperty(S,"characteristic"))
]));
fi;


if Charact=0 then AddWrds:=AddFreeWords; else
	AddWrds:=function(v,w);
	return AddFreeWordsModP(v,w,Charact);
	end;
fi;

#####################################################################
AddLst:=function(v,SM)
local x,ab;

for x in v do
if not IsBound(SM[x[2]]) then SM[x[2]]:=[]; fi;
ab:=AbsInt(x[1]);
if not IsBound(SM[x[2]][ab]) then
        SM[x[2]][ab]:=[ab,x[2]];
else
Unbind(  SM[x[2]][ab]);
fi;
od;

end;
#####################################################################


#####################################################################
Dimension:=function(i)
local D,j;

if i=0 then return 1; fi;

D:=0;

for j in [0..i] do
D:=D+DimensionR(j)*DimensionS(i-j);
od;

return D; 


end;
#####################################################################

DimPQrec:=List([1..n+1],i->[]);
#####################################################################
DimPQ:=function(p,q)
local D,j;

if (p<0) or (q<0) then return 0; fi;
if not IsBound(DimPQrec[p+1][q+1]) then 
D:=0;
for j in [0..q] do
D:=D+DimensionR(p+q-j)*DimensionS(j);
od;
DimPQrec[p+1][q+1]:=D;
fi;
return DimPQrec[p+1][q+1];
end;
#####################################################################

#####################################################################
Int2Pair:=function(i,p,q)       #Assume that x<=DimR(p)*DimS(q).
local s,r,x;
                           	#The idea is that the generator f_i in F
				#corresponds to a tensor (e_r x e_s)
x:=AbsInt(i)-DimPQ(p+1,q-1);     #with e_r in R_p, e_s in S_q. If we
s:= x mod DimensionS(q);                #input i we get output [r,s].
r:=(x-s)/DimensionS(q);

if s=0 then return [SignInt(i)*r,DimensionS(q)];
else return [SignInt(i)*(r+1),s]; fi;

end;
#####################################################################

#####################################################################
Pair2Int:=function(x,p,q)
local y;                        #Pair2Int is the inverse of Int2Pair.

y:=[AbsInt(x[1]),AbsInt(x[2])];
return
SignInt(x[1])*SignInt(x[2])*((y[1]-1)*DimensionS(q)+y[2]+DimPQ(p+1,q-1));


end;
#####################################################################

HtpyRecord:=[];
for p in [0..n] do
HtpyRecord[p+1]:=[];
for q in [0..n-p] do
HtpyRecord[p+1][q+1]:=[];
for j in [1..DimPQ(p,q)] do
HtpyRecord[p+1][q+1][j]:=[];
od;
od;
od;

if not BoolE then
#####################################################################
Htpy:=function(p,q,x)
local tensor, t,g, r, s,AB;

AB:=AbsInt(x[1]);

if not IsBound(HtpyRecord[p+1][q+1][AB][x[2]]) then

tensor:=Int2Pair(AB,p,q);

g:=NEhomN(Mult(InvE(GmapE(EhomG(x[2]))),x[2] ));
t:=GmapE(EhomG(x[2]));
r:=ShallowCopy(HomotopyS(q,[tensor[2],g]));
Apply(r,y->[y[1],NhomE(y[2])]);
Apply(r,y->[Pair2Int([tensor[1],y[1]],p,q+1),Mult(t,y[2])]);
HtpyRecord[p+1][q+1][AB][x[2]]:=r;
fi;

if SignInt(x[1])>0 then
return HtpyRecord[p+1][q+1][AB][x[2]];
else
return NegateWord(HtpyRecord[p+1][q+1][AB][x[2]]);
fi;

end;
#####################################################################
else
#####################################################################
Htpy:=function(p,q,x)
local tensor, t,g, r, s,AB;

AB:=AbsInt(x[1]);

if not IsBound(HtpyRecord[p+1][q+1][AB][x[2]]) then

tensor:=Int2Pair(AB,p,q);

g:=NEhomN(MT[InvE(GmapE(EhomG(x[2])))][x[2]] );
t:=GmapE(EhomG(x[2]));
r:=ShallowCopy(HomotopyS(q,[tensor[2],g]));
Apply(r,y->[y[1],NhomE(y[2])]);
Apply(r,y->[Pair2Int([tensor[1],y[1]],p,q+1),MT[t][y[2]]]);
HtpyRecord[p+1][q+1][AB][x[2]]:=r;
fi;

if SignInt(x[1])>0 then
return HtpyRecord[p+1][q+1][AB][x[2]];
else
return NegateWord(HtpyRecord[p+1][q+1][AB][x[2]]);
fi;

end;
#####################################################################
fi;




if not Charact=2 then
#####################################################################
CompHtpy:=function(p,q,b)
local w, r;

r:=[];
for w in b do
r:=AddWrds(Htpy(p,q,w),r);
od;

return r;
end;
#####################################################################
else
#####################################################################
CompHtpy:=function(p,q,b)
local w, x, r, SM,j;

SM:=[];
for w in b do
AddLst(Htpy(p,q,w),SM);
od;

r:=[];
for x in SM do
for j in x do
if not j=0 then
Add(r,j);fi;
od;
od;

return r;
end;
#####################################################################
fi;
#####################################################################
Del:=function(k,p,q,x)
local b,i,r,v,w,tensor,Record,Ab,j,SM,y;  #Assume that 1 <= x <= DimR(p)*DimS(q)

Ab:=AbsInt(x);
if not DelRecord[k+1][p+1][q+1][Ab] = 0 then 
    if SignInt(x)=1 then return DelRecord[k+1][p+1][q+1][Ab]; 
    else return NegateWord( DelRecord[k+1][p+1][q+1][Ab] );
    fi;
fi;

	#############################################################
	Record:=function();
	if SignInt(x)=1 then
	DelRecord[k+1][p+1][q+1][Ab]:=v;
	else
	DelRecord[k+1][p+1][q+1][Ab]:=NegateWord(v);
	fi;
	end;
	#############################################################

tensor:=Int2Pair(x,p,q);

if k=0 then
   b:=BoundaryS(q,tensor[2]);
   v:= List(b,v->[Pair2Int([tensor[1],v[1]],p,q-1),NhomE(v[2])]);
   Record();
   return v;
fi;

if k=1 then
   if q=0 then
         if p>0 then
            b:=ShallowCopy(BoundaryR(p,-tensor[1]));
            v:=List(b,v->[Pair2Int([v[1],tensor[2]],p-1,q),GmapE(v[2])]);
	    Record();
	    return v;
	 else return []; 
	 fi;
   else
         if p>0 then
            v:=CompHtpy(p-1,q-1,CompDel(1,p,q-1,Del(0,p,q,-x)));
	    Record();
	    return v;
         else return []; 
	 fi;
   fi;
fi;

if k>1 then
    if p>(k-1) then

       r:=[];
       for i in [1..k] do
###
### THE NEXT LINE TAKES UP ALL THE TIME!
       r:=AddWrds(CompDel(i,p-k+i,q+k-i-1,Del(k-i,p,q,-x)),r);
###
###
       od;
       v:= CompHtpy(p-k,q+k-2,r);
       Record();
       return v;
    else return [];
    fi;
fi;

end;
#####################################################################


if not BoolE then
#####################################################################
CompDel:=function(k,p,q,b)
local r,v,w,x, map,SM,j,y;

map:=function(x);
return Del(k,p,q,x);
end;

###############
if not Charact=2 then
r:=[];
for v in b do
w:=ShallowCopy(map(v[1])); 
Apply(w,y->[y[1],Mult(v[2],y[2])]);
r:=AddWrds(w,r);
od;

else
SM:=[];
for v in b do
w:=ShallowCopy(map(v[1])); 
Apply(w,y->[y[1],Mult(v[2],y[2])]);
AddLst(w,SM);
od;

r:=[];
for y in SM do
for j in y do
#if not j=0 then
Add(r,j);
#fi;
od;
od;
fi;
###############

return r;
end;
#####################################################################
else
#####################################################################
CompDel:=function(k,p,q,b)
local r,v,w,x, map,SM,j,y;

map:=function(x);
return Del(k,p,q,x);
end;

###############
if not Charact=2 then
r:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],MT[v[2]][y[2]]]);
r:=AddWrds(w,r);
od;

else
SM:=[];
for v in b do
w:=ShallowCopy(map(v[1]));
Apply(w,y->[y[1],MT[v[2]][y[2]]]);
AddLst(w,SM);
od;

r:=[];
for y in SM do
for j in y do
#if not j=0 then
Add(r,j);
#fi;
od;
od;
fi;
###############

return r;
end;
#####################################################################
fi;

DelRecord:=[];
for l in [0..n] do
DelRecord[l+1]:=[];
for p in [0..n] do
DelRecord[l+1][p+1]:=[];
for q in [0..n-p] do
DelRecord[l+1][p+1][q+1]:=[];
for j in [DimPQ(p+1,q-1)+1..DimPQ(p,q)] do
DelRecord[l+1][p+1][q+1][j]:=0;
od;
od;
od;
od;


PseudoBoundary:=[];
HorizontalPseudoBoundary:=[];
for k in [1..n] do
PseudoBoundary[k]:=[];
HorizontalPseudoBoundary[k]:=[];
   for q in [0..k] do
   p:=k-q;
      for j in [DimPQ(p+1,q-1)+1..DimPQ(p,q)] do

      r:=[];rr:=[];

         for l in [0..p] do
	 dl:=Del(l,p,q,j);    
         r:=AddWrds(dl,r);
	 if l>0 then rr:=AddWrds(dl,rr);fi;
         od;

         Add(PseudoBoundary[k],[r,p]);    #I'm pretty sure it's p and not q
	 Add(HorizontalPseudoBoundary[k],rr);
      od;
  od;
od;


#####################################################################
Boundary:=function(k,j);
if k=0 then return []; 
else
	if SignInt(j)=1 then return PseudoBoundary[k][j][1]; 
	else return NegateWord(PseudoBoundary[k][-j][1]);
	fi;
fi;
end;
#####################################################################










#######START WORKING ON THE CONTRACTING HOMOTOPY####################

#####################################################################
EmapN:=function(x);

#return NEhomN(Mult(x,InvE(GmapE(EhomG(x))))); 
return NEhomN(Mult(InvE(GmapE(EhomG(x))),x)); #Added 9th March 2011

end;
#####################################################################



#####################################################################
Int2Vector:=function(k,j)
local tmp,p,q;

p:=k;q:=0;
while j>=DimPQ(p,q)+1 do
p:=p-1;q:=q+1;
od;                             #p,q are now computed from k,j

tmp:=Int2Pair(j,p,q);

return [p,q,tmp[1],tmp[2]];
end;
#####################################################################

#####################################################################
Vector2Int:=function(p,q,r,s);
return Pair2Int([r,s],p,q);
end;
#####################################################################


#####################################################################
HorizontalBoundaryGen:=function(k,y)
local horizontal;

if k=0 then return [];
else
        if SignInt(y[1])=1 then horizontal:=ShallowCopy(
	HorizontalPseudoBoundary[k][y[1]]);
        else horizontal:=NegateWord(ShallowCopy(
	HorizontalPseudoBoundary[k][-y[1]]));
        fi;
fi;

Apply(horizontal,x->[x[1],Mult(y[2],x[2])]);

return horizontal;
end;
#####################################################################

#####################################################################
HorizontalBoundaryWord:=function(n,w)
local x, bnd;

bnd:=[];
for x in w do
Append(bnd,HorizontalBoundaryGen(n,x));
od;
return bnd;

end;
#####################################################################

#####################################################################
HomotopyGradedGen:=function(g,p,q,r,s,bool)    #Assume EltsE[g] exists!
local aa,hty, hty1, Eg, Eg1, Eg2, g1, g2;       #bool=true for vertical homotopy
#This function seems to work! But I should really check the maths again!!

#Eg:=EltsE[g];
#Eg1:=Image(EhomG,Eg);
#Eg2:=Image(EmapN,Eg);
#g2:=Position(S!.elts,Eg2);
#g1:=Position(R!.elts,Eg1);
#Eg1:=Image(GmapE,Eg1);
#Eg2:=Image(NhomE,Eg2);

g1:=EhomG(g);
g2:=EmapN(g);
Eg1:=GmapE(g1);
Eg2:=NhomE(g2);



hty:=HomotopyS(q,[s,g2]);
#Apply(hty,x->[ Vector2Int(p,q+1,r,x[1]), Image(NhomE,S!.elts[x[2]])]);
#Apply(hty,x->[ x[1], Elts2Int(Eg1*x[2])]);
Apply(hty,x->[ Vector2Int(p,q+1,r,x[1]), NhomE(x[2])]);
Apply(hty,x->[ x[1], Mult(Eg1,x[2])]);


if (p=0 and q>0) or bool then return hty; fi;

if p>0 then
hty1:=HomotopyOfWord(p+q,ShallowCopy(HorizontalBoundaryWord(p+q+1,hty)),false);
Append(hty, NegateWord(hty1));

fi;

if q>0 then return hty; fi;


hty1:=HomotopyR(p,[r,g1]);
#Apply(hty1,x->[ Vector2Int(p+1,q,x[1],s), Image(GmapE,R!.elts[x[2]])]);
#Apply(hty1,x->[ x[1], Elts2Int(x[2])]); #Here
Apply(hty1,x->[ Vector2Int(p+1,q,x[1],s), GmapE(x[2])]);


hty1:=NegateWord(hty1); ####added
Append(hty,hty1);

hty1:=HomotopyOfWord(p+q,ShallowCopy(HorizontalBoundaryWord(p+q+1,hty1)),true);

Append(hty,NegateWord(hty1));

hty1:=HomotopyOfWord(p+q,ShallowCopy(HorizontalBoundaryWord(p+q+1,hty1)),false);

Append(hty,hty1);       #I think this perturbation term is always zero and
                        #thus not necessary.

return hty;
end;
#####################################################################

#####################################################################
Homtpy:=function(n,x,bool)
local vec,a;


a:=AbsInt(x[1]);
vec:=Int2Vector(n,a);
if SignInt(x[1])=1 then
return HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool);
else
return NegateWord(HomotopyGradedGen(x[2],vec[1],vec[2],vec[3],vec[4],bool));
fi;
end;
#####################################################################

#####################################################################
HomotopyOfWord:=function(n,w,bool)
local x, hty;

hty:=[];
for x in w do
Append(hty,Homtpy(n,x,bool));
od;
return hty;

end;
#####################################################################

#HomotopyRec:=[];
#for i in [1..n] do  #n=Length
#HomotopyRec[i]:=[];
#for j in [1..Dimension(i-1)] do
#HomotopyRec[i][j]:=[];
#od;od;


#####################################################################
FinalHomotopy:=function(n,x);
return Homtpy(n,x,false);
end;
#####################################################################

if HomotopyR=fail or HomotopyS=fail then
FinalHomotopy:=fail;
fi;

#########FINISHED WORKING ON THE CONTRACTING HOMOTOPY##############
grp:=Group(EltsE);



################spectral sequence requirements##################

FilteredLength:=Length(R);

##################################################
FilteredDimension:=function(r,i);
return Length(Filtered(List(PseudoBoundary[i],x->x[2]),y->y<=r));
end;
##################################################


return Objectify(HapResolution,
	    rec(
	    dimension:=Dimension, 
            filteredDimension:=FilteredDimension,
	    boundary:=Boundary, 
	    homotopy:=FinalHomotopy, 
	    elts:=EltsE, 
	    group:=grp,
	    vectorToInt:=Vector2Int,
	    intToVector:=Int2Vector,
            pseudoBoundary:=PseudoBoundary,
	    properties:=
	    [["type","resolution"],
	     ["length",n],
             ["filtration_length",FilteredLength],
             ["initial_inclusion",false],
	     ["characteristic",Charact],
	     ["isTwistedTensorProduct",true]
	      ]));
	     
end);
#####################################################################

