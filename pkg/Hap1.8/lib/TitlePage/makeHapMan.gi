InstallGlobalFunction(MakeHAPManual,
function()
local  VISCRIPT, HAPDOC, cn;

cn:=Concatenation;

############################################################
#IF NECESSARY, CHANGE "VISCRIPT" AND "HAPDOC" TO THE CORRECT PATHS 
#
VISCRIPT:=cn(GAP_ROOT_PATHS[1],"pkg/Hap1.8/lib/TitlePage/viscript ");
HAPDOC:=cn(GAP_ROOT_PATHS[1],"pkg/Hap1.8/doc/");
#
############################################################

MakeGAPDocDoc(HAPDOC,"HapMan.xml",[],"HAP");
Exec( cn(VISCRIPT, cn(HAPDOC,"chap0.html")) ); 
Exec( cn(VISCRIPT, cn(HAPDOC,"chap1.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap2.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap3.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap4.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap5.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap6.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap7.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap8.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap9.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap10.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap11.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap12.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap13.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap14.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap15.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap16.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap17.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap18.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chapInd.html")) );
SetHelpViewer("mozilla");

end) ;
