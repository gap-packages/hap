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
if not IsBound(DOT_PATH)
   then
    DOT_PATH:=Filename( DirectoriesSystemPrograms( ), "dot" );

if DOT_PATH=fail
   or not IsExecutableFile(DOT_PATH)
   then
    Info(InfoWarning,1,"HAP warning: Set DOT_PATH manually if needed.");
else
DOT_PATH:=Concatenation(DOT_PATH," ");
    MakeReadOnlyGlobal("DOT_PATH");
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

#############################################################
if not IsBound(ASY_PATH)
   then
    ASY_PATH:=Filename( DirectoriesSystemPrograms( ), "asy" );

if ASY_PATH=fail
   or not IsExecutableFile(ASY_PATH)
   then
    Info(InfoWarning,1,"HAP warning: Set ASY_PATH manually if needed.");
else
ASY_PATH:=Concatenation(ASY_PATH," ");
    MakeReadOnlyGlobal("ASY_PATH");
fi;
fi;
#############################################################


