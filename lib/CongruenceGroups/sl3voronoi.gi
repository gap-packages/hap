# Graham Ellis 2005-2026 
#This file was written with the help of ChatGPT


#############################################################################
#
# Voronoi reduction for (currently) G = SL(3,Z)
#
# Convention:
#
#        Q -> g^T Q g
#
# and matrices are multiplied on the RIGHT during reduction.


#############################################################################
# Voronoi generators for group G such as G=SL(3,Z)

InstallGlobalFunction(VoronoiGenerators,
function(G);

if G=SL(2,Integers) then return
  [
  [ [ 0, -1 ], [ 1, 0 ]],
  [ [ 1, 1 ], [ 0, 1 ]]
  ];
fi;

if G=SL(3,Integers) then return 
  [
  [ [ -1, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ],
  [ [  0, 0, 1 ], [ 1, 0,-1 ], [ 0, 1, 1 ] ],
  [ [  0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0,-1 ] ],
  [ [  0, 0, 1 ], [ 1, 0, 0 ], [ 0, 1, 0 ] ],
  [ [  0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ],
  [ [  0, 1, 0 ], [ 0, 0, 1 ], [ 1, 1, 0 ] ]
];
fi;

return fail;

end);
#
############################################################################


############################################################################
# Matrix in inerior of fundamental Voronoi domain

InstallGlobalFunction(VoronoiBaseMatrix,
function(G);

if G=SL(3,Integers) then
    return 
    [ [3,2,1],
      [2,4,2],
      [1,2,3] ];
fi;

return fail;

end);
#
############################################################################


#############################################################################
# Matrix corresponding to a word

InstallGlobalFunction(MatrixFromWord, 
function(gens,w)
local RightMove,ginv, A,a;

ginv := List(gens,x->x^-1);

   ###########
   RightMove := function(M,a)
   if a > 0 then
       return M*gens[a];
   else
       return M*ginv[-a];
   fi;
   end;
   ###########

A := IdentityMat(Length(gens[1]));

for a in w do
    A := RightMove(A,a);
od;

return A;

end);
#
############################################################################


#############################################################################
# Factorize a matrix as a word in voronoi generators

InstallGlobalFunction(VoronoiFactorize,
function(G,M)

local w,gens,CheckWord, elts,nums;

if Order(M)=1 then return []; fi;

####################################
if G=SL(2,Integers) then 
elts:=VoronoiGenerators(G);
elts:=[elts[1],elts[1]^2,elts[1]^3,elts[2],elts[2]^-1];
nums:=[1,[1,1],-1,2,-2];
w:= List(AsWordInSL2Z(M), m->Position(elts,m));   
w:=Filtered(w,x-> not x=fail);
return Flat(List(w,i->nums[i]));
fi;
####################################

w := VoronoiWord(G,M);
gens:=VoronoiGenerators(SL(3,Integers));

##############
CheckWord := function(M,w)
return MatrixFromWord(gens,w) = M;
end;
##############

if not CheckWord(M,w) then
Print("Check = ",CheckWord(M,w),"\n");
return fail;
fi;

return w;

end);
#
############################################################################


#############################################################################
# Voronoi coordinates
# Returns the coordinates c_i with Q = Sum_i {c_i*R_i} with R_i=v_i*v_i^t
# where v_i are the vectors yielding minimum values of the perfect form.

InstallGlobalFunction(VoronoiCoordinates,
function(G,Q)
local a,b,c,d,e,f;

if G=SL(3,Integers) then
a := Q[1][1];
b := Q[2][2];
c := Q[3][3];

d := Q[1][2];
e := Q[1][3];
f := Q[2][3];

return [
    a-d,
    b-d-f+e,
    c-f,
    d-e,
    f-e,
    e
];
fi;

return fail;
end);
#
#############################################################################


#############################################################################
# Generate the 24-element stabilizer and find a word for every element.
#
# We only need the stabilizer generators above, so this is very small.

InstallGlobalFunction(VoronoiBuildStabilizer,
function(G)
local stab,words,frontier,newfrontier,InverseWord,
      A,w,B,u,key,gens,StabGenerators;

if not G=SL(3,Integers) then return fail; fi;

StabGenerators := [
    [ 1,-3 ],
    [ 2, 2 ],
    [ 1,-2,-3 ]
];

stab := rec();
words := rec();
gens:=VoronoiGenerators(SL(3,Integers));

######################
InverseWord := function(w)
return List(Reversed(w),x -> -x);
end;
######################

A := IdentityMat(3);
key := String(A);

stab.(key) := A;
words.(key) := [];

frontier := [A];

while Length(frontier) > 0 do

    newfrontier := [];

    for A in frontier do

        w := words.(String(A));

        for u in Concatenation(
                    StabGenerators,
                    List(StabGenerators,InverseWord)) do

            B := A*MatrixFromWord(gens,u);
            key := String(B);

            if not IsBound(stab.(key)) then

                stab.(key) := B;
                words.(key) := Concatenation(w,u);

                Add(newfrontier,B);

            fi;

        od;

    od;

    frontier := newfrontier;

od;

return [stab,words];

end);
#
#############################################################################


#############################################################################
# Voronoi reduction.
#
# Starting with M, repeatedly choose a negative Voronoi coordinate c_i
# and replace M by M*g_i.
#
# Returns:
#
#   [ Mred, word ]
#
# where
#
#       Mred = M * MatrixFromWord(gens,word)
#
# and Mred lies in the fundamental cone.

InstallGlobalFunction(VoronoiReduce,
function(G,M,VoronoiCoordinates)

local A,Anew,VoronoiHeight,MCoordinates,c,i,H,Hnew,w,candidates,j,S0,gens;

    A := M;
    w := [];

    S0 := VoronoiBaseMatrix(G);
    gens:=VoronoiGenerators(G);

########################
VoronoiHeight := function(M)
local Q,c;
Q := TransposedMat(M)*S0*M;
c := VoronoiCoordinates(G,Q);
return Sum(c);
end;
########################

########################
MCoordinates := function(M)
return VoronoiCoordinates(G,TransposedMat(M)*S0*M);
end;
########################
    while true do

        c := MCoordinates(A);
        H := Sum(c);

        candidates := [];

        for j in [1..Length(VoronoiCoordinates(G,S0))] do

            if c[j] < 0 then

                Anew := A * gens[j];
                Hnew := VoronoiHeight(Anew);

                if Hnew < H then
                    Add(candidates,[j,Hnew]);
                fi;

            fi;

        od;

        if Length(candidates) = 0 then
            break;
        fi;

        Sort(candidates,
             function(x,y)
                 return x[2] < y[2];
             end);

        i := candidates[1][1];

        A := A * gens[i];
        Add(w,i);

    od;

    return [A,w];

end);
#
#############################################################################


#############################################################################
# Full Voronoi expression.
#
# If
#
#       A = M * w
#
# with A in the stabilizer, then
#
#       M = A * w^-1.

InstallGlobalFunction(VoronoiWord,
function(G,M)

local R,A,w,ws,gens,StabilizerData,StabilizerWords,WordForStabilizerElement;

gens:=VoronoiGenerators(SL(3,Integers));

R := VoronoiReduce(G,M,VoronoiCoordinates);

A := R[1];
w := R[2];

StabilizerData := VoronoiBuildStabilizer(SL(3,Integers));
StabilizerWords := StabilizerData[2];

#####################################
WordForStabilizerElement := function(A)
local key;
key := String(A);
if not IsBound(StabilizerWords.(key)) then
    Error("Matrix is not in the stabilizer");
fi;
return StabilizerWords.(key);
end;
#####################################


ws := WordForStabilizerElement(A);

return Concatenation(ws,List(Reversed(w), x->-x));

end);
#
#############################################################################

