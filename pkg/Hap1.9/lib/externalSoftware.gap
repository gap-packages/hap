#############################################################
if not IsBound(POLYMAKE_PATH)
   then
    POLYMAKE_PATH:=Filename( DirectoriesSystemPrograms( ), "polymake" );

if POLYMAKE_PATH=fail
   or not IsExecutableFile(POLYMAKE_PATH)
   then
    Info(InfoWarning,1,"polymake command not found. Please set POLYMAKE_PATH by hand if you need to use HAP functions that rely on polymake.");
else
POLYMAKE_PATH:=Concatenation(POLYMAKE_PATH," ");
    MakeReadOnlyGlobal("POLYMAKE_PATH");
fi;
fi;
#############################################################

#############################################################
if not IsBound(NEATO_PATH)
   then
    NEATO_PATH:=Filename( DirectoriesSystemPrograms( ), "neato" );

if NEATO_PATH=fail
   or not IsExecutableFile(NEATO_PATH)
   then
    Info(InfoWarning,1,"neato command not found. Please set NEATO_PATH by hand if you need to use HAP functions that rely on the GraphViz function neato.");
else
NEATO_PATH:=Concatenation(NEATO_PATH," ");
    MakeReadOnlyGlobal("NEATO_PATH");
fi;
fi;
#############################################################

#############################################################
if not IsBound(DISPLAY_PATH)
   then
    DISPLAY_PATH:=Filename( DirectoriesSystemPrograms( ), "display" );

if DISPLAY_PATH=fail
   or not IsExecutableFile(DISPLAY_PATH)
   then
    Info(InfoWarning,1,"display command not found. Please set DISPLAY_PATH by hand if you intend to use HAP functions that display images.");
else
DISPLAY_PATH:=Concatenation(DISPLAY_PATH," ");
    MakeReadOnlyGlobal("DISPLAY_PATH");
fi;
fi;
#############################################################

#############################################################
if not IsBound(BROWSER_PATH)
   then
    BROWSER_PATH:=Filename( DirectoriesSystemPrograms( ), "mozilla" );

if BROWSER_PATH=fail
   or not IsExecutableFile(BROWSER_PATH)
   then
    Info(InfoWarning,1,"mozilla command not found. Please set BROWSER_PATH by hand if you intend to use HAP functions that require a web browser.");
else
BROWSER_PATH:=Concatenation(BROWSER_PATH," ");
    MakeReadOnlyGlobal("BROWSER_PATH");
fi;
fi;
#############################################################

