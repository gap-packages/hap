
#####################################################################
#####################################################################
InstallGlobalFunction (QuadraticNumber,
function(x,y,BI)
local num;
#if y=0 then return x; fi;
num  := rec( rational:=x,
             irrational:=y,
             bianchiInteger:=BI);

return Objectify(HapQuadraticNumber, num);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallGlobalFunction (QuadraticNumberConjugate,
function(x);
if IsRat(x) then return x; fi;
return QuadraticNumber(x!.rational, -x!.irrational, x!.bianchiInteger);
end);
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( NiceMonomorphism,
    "for small group of 2x2 matrices of quadratic numbers",
    [IsHap2x2matrixGroup],
100000000,
    function( G )
    local T, mono;
    if HasNiceMonomorphism(G) then return G!.NiceMonomorphism; fi;
    T:=List(MultiplicationTable(G),PermList);;
    mono:= GroupHomomorphismByImagesNC(G,Group(T),Elements(G),T);;
    SetIsInjective(mono,true);
    SetIsSurjective(mono,true);
    return mono;
    end );
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod( \+,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    if not x!.bianchiInteger=y!.bianchiInteger then return fail; fi;
    #if x!.irrational=-y!.irrational then return x!.rational + y!.rational; fi;
    return QuadraticNumber( x!.rational + y!.rational, x!.irrational + y!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \+,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    return QuadraticNumber( x + y!.rational, y!.irrational, y!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \+,
    "for quadratic number",
    [IsHapQuadraticNumber, IsRat],
    function( x, y )
    return QuadraticNumber( y + x!.rational, x!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \-,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    return QuadraticNumber( x - y!.rational, -y!.irrational, y!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \-,
    "for quadratic number",
    [IsHapQuadraticNumber, IsRat],
    function( x, y )
    return QuadraticNumber( -y + x!.rational, x!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod( \-,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    if not x!.bianchiInteger=y!.bianchiInteger then return fail; fi;
    #if x!.irrational=y!.irrational then return x!.rational - y!.rational; fi;
    return QuadraticNumber( x!.rational - y!.rational, x!.irrational - y!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( AdditiveInverseOp,
    "for quadratic number",
    [IsHapQuadraticNumber],
    function( x )
    return QuadraticNumber( -x!.rational, -x!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( InverseOp,
    "for quadratic number",
    [IsHapQuadraticNumber],
    function( y )
local BI,D;
BI:=y!.bianchiInteger;
D:=y!.rational^2 -BI*y!.irrational^2;
return  QuadraticNumber( (y!.rational )/D,
     - (y!.irrational)/D     , y!.bianchiInteger );
    end );
#####################################################################
#####################################################################


#####################################################################
#####################################################################
InstallOtherMethod( \*,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    #if x=0 then return 0; fi;
    return QuadraticNumber( x * y!.rational, x * y!.irrational, y!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \*,
    "for quadratic number",
    [IsHapQuadraticNumber, IsRat],
    function( x, y )
    #if y=0 then return 0; fi;
    return QuadraticNumber( y * x!.rational, y * x!.irrational, x!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \*,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    local p;
    if not x!.bianchiInteger=y!.bianchiInteger then return fail; fi;
    p:= QuadraticNumber( x!.rational * y!.rational + x!.bianchiInteger * x!.irrational * y!.irrational, x!.rational * y!.irrational + x!.irrational * y!.rational,  x!.bianchiInteger );
    #if IsRat(p) then return p; fi;
    #if p!.irrational = 0 then return p!.rational; fi;
    return p;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \/,
    "for quadratic number",
    [IsHapQuadraticNumber, IsRat],
    function( y, x )
    return QuadraticNumber( y!.rational/x,  y!.irrational/x, y!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \/,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    local D, BI,p;
    if not x!.bianchiInteger=y!.bianchiInteger then return fail; fi;
    BI:=x!.bianchiInteger;
    D:=y!.rational^2 -BI*y!.irrational^2;
    p:=QuadraticNumber( (x!.rational*y!.rational -BI*x!.irrational*y!.irrational)/D,
    (x!.irrational*y!.rational - x!.rational*y!.irrational)/D     , x!.bianchiInteger );
    #if p!.irrational=0 then return p!.rational; fi;
    return p;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \/,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    local D, BI;
    BI:=y!.bianchiInteger;
    D:=y!.rational^2 -BI*y!.irrational^2;
    return QuadraticNumber( (x*y!.rational )/D,
     - (x*y!.irrational)/D     , y!.bianchiInteger );
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( Sqrt,
    "for quadratic number",
    [IsHapQuadraticNumber],
    function( x );
    if x!.irrational=0 then return Sqrt(x!.rational); fi;
    return fail;
    end );
#####################################################################
#####################################################################
#
#####################################################################
#####################################################################
InstallOtherMethod( POW,
    "for quadratic number",
    [IsHapQuadraticNumber, IsInt],
    function( x, n )
    if n=0 then return QuadraticNumber(1,0,x!.bianchiInteger); fi;
    if n>0 then return x*POW(x,n-1); fi;
    if n<0 then return 1/POW(x,AbsInt(n)); fi;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \<,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    if x!.bianchiInteger <> y!.bianchiInteger then return fail; fi;
    if x!.rational <y!.rational then return true; fi;
    if x!.rational =y!.rational and x!.irrational < y!.irrational then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \<,
    "for quadratic number",
    [IsHapQuadraticNumber, IsRat],
    function( x, y )
    if x!.rational <y then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \<,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    if not y=x and not y<x  then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \=,
    "for quadratic number",
    [IsHapQuadraticNumber, IsHapQuadraticNumber],
    function( x, y )
    if x!.bianchiInteger <> y!.bianchiInteger then return false; fi;
    if x!.rational = y!.rational  and x!.irrational = y!.irrational then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \=,
    "for quadratic number",
    [IsRat, IsHapQuadraticNumber],
    function( x, y )
    if x = y!.rational  and y!.irrational=0 then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( \=,
    "for quadratic number",
    [IsHapQuadraticNumber,IsRat],
    function( x, y )
    if y = x!.rational  and x!.irrational=0 then return true; fi;
    return false;
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( One,
    "for quadratic number",
    [IsHapQuadraticNumber],
    function( x )
    return QuadraticNumber(1,0,x!.bianchiInteger);
    end );
#####################################################################
#####################################################################

#####################################################################
#####################################################################
InstallOtherMethod( Zero,
    "for quadratic number",
    [IsHapQuadraticNumber],
    function( x )
    return QuadraticNumber(0,0,x!.bianchiInteger);
    end );
#####################################################################
#####################################################################





#####################################################################
#####################################################################
HAPQuadratic:=function(x)
local Q;

Q:=rec();
if IsRat(x) then
Q.a:=x;
Q.b:=0;
Q.d:=1;
else
Q.a:=x!.rational;
Q.b:=x!.irrational;
Q.d:=1;
fi;

return Q;
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HAPNorm:=function(OQ,x);
if IsRat(x) then return x^2; fi;
if not OQ!.bianchiInteger=x!.bianchiInteger then return fail; fi;
return x!.rational^2+AbsInt(x!.bianchiInteger)*x!.irrational^2;
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
HAPSqrt:=function(n);
if IsRat(Sqrt(n)) then return Sqrt(n); fi;
return QuadraticNumber(0,1,n);
end;
#####################################################################
#####################################################################

#####################################################################
#####################################################################
QuadraticToCyclotomic:=function(a);
if IsRat(a) then return a; fi;
return a!.rational+a!.irrational*Sqrt(a!.bianchiInteger);
end;
#####################################################################
#####################################################################

