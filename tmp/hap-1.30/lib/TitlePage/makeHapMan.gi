InstallGlobalFunction(MakeHAPManual,
function()
local  grp, HAPDOC, cn;

cn:=Concatenation;

grp:=HAP_ROOT;

############################################################
#IF NECESSARY, CHANGE "HAPDOC" TO THE CORRECT PATH 
#
HAPDOC:=cn(grp{[1..Length(grp)-4]},"doc/");
#
############################################################

MakeGAPDocDoc(HAPDOC,"newHapMan.xml",[],"HAP","MathJax");

HAPDOC:=Concatenation(HAPDOC,"tutorial/");
MakeGAPDocDoc(HAPDOC,"HapTutorial.xml",[],"HAP","MathJax");

end) ;
