@@ -0,0 +1,29 @@
# YH_Analysis
Purpose: Clean data to more concisely analyze results of a Yeast Hybrid assay.

This function pulls .TXT data from a user-specified directory and outputs .CSV files for each of the baits listed by the user.

Usage:
YH_Analysis(Directory, Bait_ID, Number_of_Replicates, Number_of_Plates_per_Replicate,Number_of_Wells)

Arguments:
Directory - The folder in which the "Identification.csv" file and all of the .TXT files are located.
Bait_ID - The list of baits used in the screen. A vector of names.
Number_of_Replicates - The number of replicates for each bait. A vector of numbers.
 - If all of the baits have the same number of replicates, then this number needs to only be listed once.
Number_of_Plates_per_Replicate - The number of plates for each bait. A vector of numbers.
 - If all of the baits use the same number of plates, then this number needs to only be listed once.
Number_of_Wells - the number of wells used by each bait. A vector of numbers.
 - If all of the baits use the same number of wells, then this number needs to only be listed once.

Within your directory, you must have the following:
Text files from a yeast hybrid screen
A two columned .CSV file labelled "Identification" whose first column is AGI and second is Family, which is specific to the screen that was run.
 - This must contain the same number of AGIs as Families.
 - The number of AGIs and Families must be the same as the number of wells.
  
  
Outputs:
A .CSV file, titled as the bait and the replicate number, with columns titled 'AGI','Family','Lum-bkg','Lum-bkg (outliers)','OD600-bkg','Growth Analysis',
   'Ratio Luc/OD600','Ratio (outliers)','Fold over Cutoff', and 'Likely True Positives'

Example:

As general as can be written,
YH_Analysis (Directory = "C:/Users/Name/Yeast-Hybrid-Folder", Bait_ID = c("BaitX", "BaitY", "BaitZ"), Number_of_Replicates = c(X1, Y1, Z1), Number_of_Plates_per_Replicate = c(X2, Y2, Z2), 
             Number_of_Wells = c(X3, Y3, Z3))
 - Note the use of a forward slash "/" not backslash "\"
 - Note the use of quotations around the Directory and Bait_ID
 - Note the use of c() indicating a list of more than one item
 - Note that BaitX corresponds to Number_of_Replicates X1, Number_of_Plates_per_Replicate X2, and Number_of_Wells X3
 - Note that X1, X2, and X3 should all be numbers when using the YH_Analysis function
 
More specific,
YH_Analysis (Directory = "C:/Users/Name/Yeast-Hybrid-Folder", Bait_ID = c("BaitX", "BaitY", "BaitZ"), Number_of_Replicates = c(A1), Number_of_Plates_per_Replicate = c(A2), 
             Number_of_Wells = c(A3))
 - Note that Number_of_Replicates, Number_of_Plates_per_Replicate, and Number_of_Wells are indicated by A1, A2, and A3
 - A1, A2, and A3 should all be numbers when using the YH_Analysis function
 - If all bait in Bait_ID use the same number of replicates, plates, and wells, then you can indicate that by using one number for each one.

Most specific,
YH_Analysis (Directory = "C:/Users/Name/Yeast-Hybrid-Folder", Bait_ID = c("BaitX", "BaitY", "BaitZ"), Number_of_Replicates = c(1), Number_of_Plates_per_Replicate = c(6), 
             Number_of_Wells = c(1999))
 - Note that this is the same function as the previous example
 - Note that this is the most commonly used numbers for number of replicates, plates, and wells.