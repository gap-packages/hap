InstallGlobalFunction(MakeHAPManual,
function()
local  grp, HAPDOC, HAPTUT, cn;

cn:=Concatenation;

grp:=HAP_ROOT;

############################################################
#IF NECESSARY, CHANGE "HAPDOC" TO THE CORRECT PATH 
#
HAPDOC:=cn(grp{[1..Length(grp)-4]},"doc/");
HAPTUT:=cn(grp{[1..Length(grp)-4]},"tutorial/");
#
############################################################

MakeGAPDocDoc(HAPDOC,"newHapMan.xml",[],"HAP","MathJax");

MakeGAPDocDoc(HAPTUT,"HapTutorial.xml",[],"HAP","MathJax");

end) ;
