#############################################################################
##  
##  PackageInfo.g for HAPcocyclic                             Robert F. Morse
##  

SetPackageInfo
(  rec
   (  PackageName    :=  "HAPcocyclic",
      Subtitle       :=  "Computing with group extenions",
      Version        :=  "0.1",
      Date           :=  "11/06/2008",
      ArchiveURL     :=  "", 
      ArchiveFormats :=  ".tar.gz",

      Persons := 
      [  rec
         (  LastName      :=  "Morse",
            FirstNames    :=  "Robert",
            IsAuthor      :=  true,
            IsMaintainer  :=  true,
            Email         :=  "rfmorse@evansville.edu",
            WWWHome       :=  "http://faculty.evansville.edu",
            PostalAddress :=  
              Concatenation
              (  [  "Departement of Electrical Engineering and Computer Science\n",
                    "Univeristy of Evansville\n",
                    "1800 Lincoln Avenue\n",
                    "Evansville, IN 47722 USA" 
                 ] 
              ),
            Place         :=  "Evanville",
            Institution   :=  "University of Evansville"
         ) 
      ],

      Status := "dev",
      CommunicatedBy := "",
      AcceptDate := "",

      README_URL := "",
      PackageInfoURL := "",
      AbstractHTML := "",
      PackageWWWHome := "",
               
      PackageDoc := 
        rec
        (  BookName         :=  "HAPcocyclic",
           ArchiveURLSubset :=  ["doc", "htm"],
           HTMLStart        :=  "htm/chapters.htm",
           PDFFile          :=  "doc/manual.pdf",
           SixFile          :=  "doc/manual.six",
           LongTitle        :=  "Computing with group extensions",
           Autoload         :=  false
        ),

      Dependencies := 
        rec
        (  GAP := ">=4.3",
           NeededOtherPackages    :=  [], 
           SuggestedOtherPackages :=  [],
           ExternalConditions     :=  [] 
        ),

      AvailabilityTest := ReturnTrue,
      BannerString := "",
      Autoload := false,
      Keywords := ["group extensions"]
   )
);
