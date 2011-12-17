#############################################################################
##
##  PackageInfo.g                HAP Package                Graham Ellis 
##
#############################################################################

SetPackageInfo( rec(

  PackageName := "HAP",
  Subtitle  := "Homological Algebra Programming",
  Version := "1.2",
  Date    := "07/03/2006",
  ArchiveURL 
          := "http://hamilton.nuigalway.ie/Hap/hap1.2",
  ArchiveFormats 
          := ".tar.gz",


  Persons := [ 
    rec( 
      LastName      := "Ellis",
      FirstNames    := "Graham",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "graham.ellis@nuigalway.ie",
      WWWHome       := "http://hamilton.nuigalway.ie",
      PostalAddress := Concatenation( [
                         "Graham Ellis\n",
                         "Mathematics Department\n",
                         "NUI Galway\n",
                         "Galway\n",
                         "Ireland" ] ),
      Place         := "Galway",
      Institution   := "National University of Ireland, Galway"
    )
  ],  

  Status  := "deposited",
  #CommunicatedBy 
  #        := "none",
  #AcceptDate 
  #        := "not applicable",

  README_URL := "http://hamilton.nuigalway.ie/Hap/README.HAP",
  PackageInfoURL := "http://hamilton.nuigalway.ie/Hap/PackageInfo.g",

  AbstractHTML := 
    "This package provides some functions for group cohomology. ",

  PackageWWWHome := "http://hamilton.nuigalway.ie/Hap/www",
                  
  PackageDoc := rec(
    BookName  := "HAP",
    ArchiveURLSubset := ["doc", "www"],
    HTMLStart := "www/index.html",
    PDFFile   := "doc/manual.pdf",
    SixFile   := "doc/manual.six",
    LongTitle := "Homological Algebra Programming Package",
    Autoload := false 
  ),


  Dependencies := rec(
    GAP := ">= 4.3",
    NeededOtherPackages := [],
    SuggestedOtherPackages := [],
    ExternalConditions := [["Some functions require Polymake software",
    "http://www.math.tu-berlin.de/polymake/"]]
  ),

AvailabilityTest := ReturnTrue,

BannerString     := Concatenation( "Loading HAP ",
                            String( ~.Version ), " ...\n" ),


  Autoload := false,

#TestFile := "Hap/test/examples.test",

Keywords := [ "homology", "cohomology", "resolution", "homotopy group", 
"module of identities" ]

));

