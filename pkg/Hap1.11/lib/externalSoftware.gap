#############################################################
if not IsBound(POLYMAKE_PATH)
   then
    POLYMAKE_PATH:=Filename( DirectoriesSystemPrograms( ), "polymake" );

if POLYMAKE_PATH=fail
   or not IsExecutableFile(POLYMAKE_PATH)
   then
    Info(InfoWarning,1,"HAP warning: Set POLYMAKE_PATH manually if needed. ");
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
    Info(InfoWarning,1,"HAP warning: Set NEATO_PATH manually if needed.");
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
    Info(InfoWarning,1,"HAP warning: Set DISPLAY_PATH manually if needed.");
else
DISPLAY_PATH:=Concatenation(DISPLAY_PATH," ");
    MakeReadOnlyGlobal("DISPLAY_PATH");
fi;
fi;
#############################################################

#############################################################
if not IsBound(BROWSER_PATH)
   then
    BROWSER_PATH:=Filename( DirectoriesSystemPrograms( ), "firefox" );

if BROWSER_PATH=fail
   or not IsExecutableFile(BROWSER_PATH)
   then
    Info(InfoWarning,1,"HAP warning: Set BROWSER_PATH manually if needed.");
else
BROWSER_PATH:=Concatenation(BROWSER_PATH," ");
    MakeReadOnlyGlobal("BROWSER_PATH");
fi;
fi;
#############################################################

