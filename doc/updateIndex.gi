root:=HAP_ROOT{[1..Length(HAP_ROOT)-4]};
rootwww:=Concatenation(root,"www/SideLinks/About/");
rootdoc:=Concatenation(root,"doc/");
roottut:=Concatenation(root,"tutorial/");

L:=[
"chap1.html", "chap2.html", "chap3.html", "chap4.html", "chap5.html",
"chap6.html", "chap7.html", "chap8.html", "chap9.html", "chap10.html",
"chap11.html", 
"aboutAbelianCategories.html",    "aboutLinks.html",
"aboutArithmetic.html",           "aboutMetrics.html",
"aboutArtinGroups.html",          "aboutModPRings.html",
"aboutAspherical.html",           "aboutNonabelian.html",
"aboutNoncrossing.html",
"aboutBogomolov.html",            "aboutParallel.html",
"aboutBredon.html",               "aboutPerformance.html",
"aboutCocycles.html",             "aboutPeriodic.html",
"aboutCoefficientSequence.html",  "aboutPeripheral.html",
"aboutCohomologyRings.html",      "aboutPersistent.html",
"aboutPoincareSeries.html",
"aboutCoveringSpaces.html",       "aboutPoincareSeriesII.html",
"aboutCoverinSpaces.html",        "aboutPolytopes.html",
"aboutCoxeter.html",              "aboutQuandles2.html",
"aboutQuandles.html",
"aboutCrossedMods.html",          "aboutquasi.html",
"aboutCubical.html",              "aboutRandomComplexes.html",
"aboutRosenbergerMonster.html",
"aboutDavisComplex.html",         "aboutSchurMultiplier.html",
"aboutDefinitions.html",          "aboutSimplicialGroups.html",
"aboutExtensions.html",           "aboutSpaceGroup.html",
"aboutFunctorial.html",           "aboutSuperperfect.html",
"aboutGouter.html",               "aboutSurvey.html",
"aboutGraphsOfGroups.html",       "aboutTDA.html",
"aboutIntro.html",                
"aboutKnots.html",                "aboutTensorSquare.html",
"aboutKnotsQuandles.html",        "aboutTopology.html",
"aboutLieCovers.html",            "aboutTorAndExt.html",
"aboutLie.html",                  "aboutTwistedCoefficients.html",
];;

##########################################
IllustratedInFile:=function(fle,fn)
local file, input, bool;

if StartsWith(fle,"about") then
file:=Concatenation(rootwww,fle);
input:=InputTextFile(file);
file:=ReadAll(input);
CloseStream(input);
file:=ReplacedString(file,fn,"@");
bool:= '@' in file;
if not bool then return bool; fi;
return Concatenation("../www/SideLinks/About/",fle);
fi;

if StartsWith(fle,"chap") then
file:=Concatenation(roottut,fle);
input:=InputTextFile(file);
file:=ReadAll(input);
CloseStream(input);
file:=ReplacedString(file,fn,"@");
bool:= '@' in file;
if not bool then return bool; fi;
return Concatenation("../tutorial/",fle);
fi;


end;
##########################################

##########################################
IllustratedInFiles:=function(fn)
local files, fle, ans;
files:=[];
if Length(fn)=0 then return files; fi;
for fle in L do
ans:=IllustratedInFile(fle,fn);
if not ans=false then Add(files,ans); fi;
od;
return files;
end;
##########################################

##########################################
AddExamplesToFile:=function(fle)
local fn, file, str, input, out, i, m ;

file:=Concatenation(rootdoc,fle);
input:=InputTextFile(file);
str:=ReadAll(input);
CloseStream(input);

str:=ReplacedString(str,"<P/><B>Examples","@gyntaf");
str:=ReplacedString(str,"</Description>","@</Description>");
str:=SplitString(str,['@']);
str:=Filtered(str,s->not StartsWith(s,"gyntaf"));
str:=Concatenation(str);

str:=ReplacedString(str,"<ManSection>","@");
str:=SplitString(str,['@']);

for i in [2..Length(str)] do
m:=StructuralCopy(str[i]);
#adjust man section here ################

fn:=ReplacedString(m,"<Func Name=\"","@");

#if not '@' in fn then
fn:=ReplacedString(fn,"<Oper Name=\"","@");
#fi;

#if not '@' in fn then
fn:=ReplacedString(fn,"<Var Name=\"","@");
#fi;

#if not '@' in fn then
fn:=ReplacedString(fn,"\" Arg","@");
#fi;


fn:=SplitString(fn,['@']);
fn:=fn[2];

fn:=IllustratedInFiles(fn);
fn:=List(fn,a->Concatenation("<URL><Link>",a,"</Link><LinkText>",String(Position(fn,a)),"</LinkText></URL>&nbsp;,\n")  );
if Length(fn)>0 then
fn[Length(fn)][Length(fn[Length(fn)])-1]:=' ';
fi;
fn:=Concatenation(fn);
m:=ReplacedString(m,"</Description>",Concatenation("<P/><B>Examples:</B> ",fn,"\n</Description>\n"));

#adjust man section here ################
str[i]:=Concatenation("<ManSection>",m);
od;

str:=Concatenation(str);
out := OutputTextFile( file, false );;
SetPrintFormattingStatus(out, false);
PrintTo(out,str);
CloseStream(out);

return true;
end;
##########################################

##########################################
AddExamplesToUndocumentedFile:=function(fle)
local fn, file, str, input, out, i, m ;

file:=Concatenation(rootdoc,fle);
input:=InputTextFile(file);
str:=ReadAll(input);
CloseStream(input);

#str:=ReplacedString(str,"<Br/>\n<B>Examples","@gyntaf");
str:=ReplacedString(str,"&nbsp;&nbsp;&nbsp;&nbsp;<B>Examples","@gyntaf");

str:=SplitString(str,['@']);
str:=Filtered(str,s->not StartsWith(s,"gyntaf"));
str:=Concatenation(str);

str:=ReplacedString(str,"<C>","@");
str:=SplitString(str,['@']);

for i in [2..Length(str)] do
m:=StructuralCopy(str[i]);
#adjust man section here ################

fn:=ReplacedString(m,"</C>","@");
fn:=SplitString(fn,['@']);
fn:=fn[1];
if Length(fn)>0 then
fn:=IllustratedInFiles(fn);
fi;

fn:=List(fn,a->Concatenation("<URL><Link>",a,"</Link><LinkText>",String(Position(fn,a)),"</LinkText></URL>&nbsp;,\n")  );
if Length(fn)>0 then
fn[Length(fn)][Length(fn[Length(fn)])-1]:=' ';
fi;
fn:=Concatenation(fn);

m:=ReplacedString(m,"</C>",Concatenation("</C>&nbsp;&nbsp;&nbsp;&nbsp;<B>Examples:</B> ",fn,"<Br/>\n"));

#adjust man section here ################
str[i]:=Concatenation("<C>",m);
od;

str:=Concatenation(str);
if not EndsWith(str,"</Chapter>\n") then
Append(str,"</Section>\n</Chapter>\n");
fi;
out := OutputTextFile( file, false );;
SetPrintFormattingStatus(out, false);
PrintTo(out,str);
CloseStream(out);

return true;
end;
##########################################


M:=["newNewCellComplexes.xml",
"newNewResolutions.xml",
"newNewGroups.xml",
"newNewParallel.xml",
"newResolutions.xml",
"newModuleResolutions.xml",
"newInducedChainMaps.xml",
"newFunctors.xml",
"newChainComplexes.xml",
"newSparse.xml",
"newHomology.xml",
"newPoincare.xml",
"newRings.xml",
"newHAPprime.xml",
"newNonabelian.xml",
"newLie.xml",
"newPresentations.xml",
"newOrbits.xml",
"newCocycles.xml",
"newWords.xml",
"newFpgmodules.xml",
"newMeataxe.xml",
"newGouter.xml",
"newCat1groups.xml",
"newSimplicialGroups.xml",
"newCoxeter.xml",
"newTorsionSubcomplexes.xml",
"newSimplicial.xml",
"newCubical.xml",
"newCW.xml",
"newKnots.xml",
"newFunctionsKnotsQuandles.xml",
"newMetrics.xml",
"newCategories.xml",
"newPseudolists.xml",
"newParallel.xml",
"newAccess.xml",
"newMiscellaneous.xml"
];;

for fle in M do
AddExamplesToFile(fle);
od;
AddExamplesToUndocumentedFile("Undocumented.xml");
MakeGAPDocDoc(rootdoc,"newHapMan.xml",[],"HAP","MathJax");
QUIT;
