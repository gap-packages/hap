#############################################################################
##
#W  init.g                   HAPcocyclic                      Robert F. Morse  
##
##    
##  $Id: init.g,v 1.2 2008-07-11 21:09:02 unialg Exp $

#############################################################################
##
##  Declare the package
##
DeclarePackage
(  "hapcocyclic", 
   "0.1", 
   function() 

   if TestPackageAvailability( "hap", "1.8.6.1" ) = fail then
       Info
       (  InfoWarning, 
          1,
          "Loading the HAPcoyclic package: package HAP must be available" 
       );
       return fail;
   fi;

   return true; 
   end 
);

DeclarePackageDocumentation( "hapcocyclic", "doc" );

#############################################################################
##
##  Load the HAP if it is not loaded already  
##
if IsList( TestPackageAvailability( "hap", "1.8.6.1" ) ) then
    HideGlobalVariables( "BANNER" );
    BANNER := false;
    LoadPackage( "hap" );
    UnhideGlobalVariables( "BANNER" );
fi;

#############################################################################
##
##  The Banner
##
if BANNER and not QUIET then
    ReadPkg("hapcocyclic", "gap/banner.g");
fi;

#############################################################################
##
##  Read .gd files
##
ReadPkg( "HAP", "lib/HapCocyclic/gap/cocycle.gd"      );
ReadPkg( "hap", "lib/HapCocyclic/gap/ccgroup.gd"      );
ReadPkg( "hap", "lib/HapCocyclic/gap/ccelms.gd"       );


#############################################################################
##
#H
##  $Log: init.g,v $
##  Revision 1.2  2008-07-11 21:09:02  unialg
##
##  First beta release. RFM
##
##  Revision 1.1  2008-06-11 17:43:56  unialg
##  File Creation. RFM
##
