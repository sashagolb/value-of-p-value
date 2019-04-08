# value-of-p-value
# the pvalues-2019-april.c software calculates p-values for classification and category variables
# see publication "value of p-value" in Molecular Informatics (2019)
# download the file to a working directory and create an executable in working directory:
# for linux: gcc pvalues-2019-april.c -o pvalues9.exe -lm
# for cygwin: gcc pvalues-2019-april.c -o pvalues9.exe -lm
# and use it as follows:
# USAGE: ./pvalues9 PARAMETERS (ALL IN ONE LINE)
# OUTPUTFILE
#  1 (FOR CLASSIFICATION) 2 (FOR CATEGORY)
#  NUMBER OF CLASSES OR CATEGORIES: LET IT BE N
#  FOR EACH CLASS OR CATEGORY:
#    PROBABILITIES TO ASSIGN A COMPOUND TO CLASS 1 OR CATEGORY 1
#    ...
#    PROBABILITIES TO ASSIGN A COMPOUND TO CLASS N OR CATEGORY N

