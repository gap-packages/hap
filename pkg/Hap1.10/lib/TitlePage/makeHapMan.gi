InstallGlobalFunction(MakeHAPManual,
function()
local  grp, VISCRIPT,  VIS, HAPDOC, HAPPDFDOC, cn;

cn:=Concatenation;

grp:=HAP_ROOT;

############################################################
#IF NECESSARY, CHANGE "VISCRIPT" AND "HAPDOC" TO THE CORRECT PATHS 
#
VISCRIPT:=cn(grp,"TitlePage/viscript ");
VIS:=cn(grp,"TitlePage/vis ");
HAPDOC:=cn(grp{[1..Length(grp)-4]},"doc/");
HAPPDFDOC:=cn(grp{[1..Length(grp)-4]},"pdfdoc/");
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
Exec( cn(VISCRIPT, cn(HAPDOC,"chap19.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap20.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap21.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap22.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap23.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap24.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap25.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap26.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap27.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap28.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chap29.html")) );
Exec( cn(VISCRIPT, cn(HAPDOC,"chapInd.html")) );

Exec( cn("cp ",HAPDOC,"/*.xml ",HAPPDFDOC)); 
Exec( cn(VIS, cn(HAPPDFDOC,"Resolutions.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"ModuleResolutions.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"InducedChainMaps.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Functors.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"ChainComplexes.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Homology.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Poincare.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Rings.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"HAPprime.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Nonabelian.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Lie.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Presentations.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Orbits.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Cocycles.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Words.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Fpgmodules.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Meataxe.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Gouter.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Cat1groups.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"SimplicialGroups.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Coxeter.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Simplicial.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Cubical.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Categories.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Pseudolists.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Parallel.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Access.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Miscellaneous.xml")) );
Exec( cn(VIS, cn(HAPPDFDOC,"Sparse.xml")) );
MakeGAPDocDoc(HAPPDFDOC,"HapMan.xml",[],"HAP");
Exec( cn("rm ",HAPPDFDOC,"*.xml "));
Exec( cn("rm ",HAPPDFDOC,"*.html "));
Exec( cn("mv ",HAPPDFDOC,"*.txt ",HAPDOC));
Exec( cn("mv ",HAPPDFDOC,"*.pdf ",HAPDOC));
Exec( cn("mv ",HAPPDFDOC,"*.six ",HAPDOC));
Exec( cn("mv ",HAPPDFDOC,"HapMan.* ",HAPDOC));




#SetHelpViewer("BROWSER_PATH");

end) ;
