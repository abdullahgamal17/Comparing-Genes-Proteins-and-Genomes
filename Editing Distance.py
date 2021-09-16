# Levenshtein introduced edit distance but did not describe an algorithm for computing it, which we leave to you.

# Edit Distance Problem: Find the edit distance between two strings.

# Input: Two strings.
# Output: The edit distance between these strings.
# Code Challenge: Solve the Edit Distance Problem.

# Extra Dataset and Test Datasets

# Sample Input:

# PLEASANTLY
# MEANLY
# Sample Output:

# 5
# You have an unlimited number of attempts.
# Time limit: 5 mins

import Mymodule as mm

def Levenshteins_Edit_Distance(firstword,secondword):
    #first word is columns , second word is rows
    rows = len(secondword)+1 # add extra row 
    cols = len(firstword)+1 # add extra column
    dptable = mm.init_2d_array(rows,cols) # makes a zero 2d array in python
    for i in range(rows):
        for j in range(cols):
            char1 = firstword[j-1]
            char2 = secondword[i-1]

            if i == 0:
                dptable[i][j] = j
            
            elif j == 0:
                dptable[i][j] = i
            
            elif char1 == char2 :
                dptable[i][j] = dptable[i-1][j-1]

            else:
                up = dptable[i-1][j]
                left = dptable[i][j-1]
                corner = dptable[i-1][j-1]
                mini = min((up,left,corner))
                dptable[i][j] = 1 + mini 
            
    return dptable

str1 = "PLEASANTLY"
str2 = "MEANLY"
dptable = Levenshteins_Edit_Distance(str1, str2)
#mm.print_2d_array(dptable)
ans = dptable[-1][-1]
print(f"{ans}")