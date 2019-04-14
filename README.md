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
#  FOR EACH CLASS OF CATEGORY, SUM OF PROBABILITIES SHOULD BE 1.
#
#  USAGE EXAMPLE:
#    ./pvalues9 9.out 2 0.5 0.5 0.4 0.6 
#
#  OUTPUT IN THE FILE IS ALWAYS COPIED ON THE SCREEN
#
#  OUTPUT FORMAT
#  TABLES ARE BUILT FOR ERRORS IN EACH CLASS OR CATEGORY AND FOR THE TOTAL DISTRIBUTION OF ERRORS.
#  TABLES ARE BUILT FOR EACH CLASS AND A TOTAL DISTRIBUTION OF ERRORS IS ALSO CALCULATED.
#  TABLE FOR A CLASS CONTAINS FIVE COLUMNS.
#  COLUMN 1: ERROR.
#  COMBINATIONS: THE NUMBER OF COMBINATIONS OF COMPOUNDS GIVING AN ERROR FROM COLUMN 1. 
#    For example, FOR 5 compounds in class 1 10 combinations gives error of 2.
#  PROBABILITY: Probabilities of errors in column 1 (see also formula 3a in the main text). 
#    For example, FOR 5 compounds in class 1, error of 0 has probability of 3.125e-2 (see the main text).
#  CUMULATIVE (COLUMN 4) is sum of probabilities giving error not above the given.
#    Density
#    Cumulative density

