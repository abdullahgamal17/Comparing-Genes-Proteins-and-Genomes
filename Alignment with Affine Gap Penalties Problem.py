#Alignment with Affine Gap Penalties Problem
# Input: Two amino acid strings v and w (each of length at most 100).
# Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.
# Download BLOSUM62 scoring matrix

# Extra Dataset and Test Datasets

# Sample Input:

# PRTEINS
# PRTWPSEIN
# Sample Output:

# 8
# PRT---EINS
# PRTWPSEIN-
# You have an unlimited number of attempts.
# Time limit: 5 mins

#Blosum62 Matrix

#    A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
# A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
# C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
# D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
# E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
# F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
# G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
# H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
# I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
# K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
# L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
# M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
# N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
# P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
# Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
# R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
# S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
# T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
# V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
# W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
# Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7





BLOSUM62 = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0,  'W': -3, 'V': 0, 'Y': -2},

          'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},

          'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},

          'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1,  'W': -4, 'V': -3, 'Y': -3},

          'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},

          'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},

          'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},

          'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},

          'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1,  'W': -3, 'V': -2, 'Y': -2},

          'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1,  'W': -1, 'V': 1, 'Y': -1},

          'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1,  'W': -2, 'V': 1, 'Y': -1},

          'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},

          'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2,'V': -2, 'Y': -1},

          'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T':  -1, 'W': -4, 'V': -2, 'Y': -3},

          'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},

          'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1,  'W': -3, 'V': -3, 'Y': -2},

          'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5,  'W': -2, 'V': 0, 'Y': -2},

          'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},

          'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},

          'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T':  -2, 'W': 2, 'V': -1, 'Y': 7}
          }

# To Use this dictionary and find a value type Blosum62[First charachter][second character]

#Initializes 2d Array with given size of rows and columns to zero
def init_2d_array(rows,cols):
    _array_2d = [[0 for i in range(cols)] for j in range(rows)]
    return _array_2d

#prints a 2D array
def print_2d_array(_2d_array):
    rows = len(_2d_array)
    cols = len(_2d_array[0])
    for i in range(rows):
        for j in range(cols):
            print(str(_2d_array[i][j]).center(5),end=" ")
        print("")

def add_indel(string,pos):
    obj = string[:pos] + "-" + string[pos:]

    return obj

#gets a value from BLOSUM62
def getval(char1,char2):
    return BLOSUM62[char1][char2]

def showinfo(lower,middle,upper,backtrack):
    print("lower".center(50) + "\n")
    print_2d_array(lower)
    print("")
    print("middle".center(50) + "\n")
    print_2d_array(middle)
    print("")
    print("upper".center(50) + "\n")
    print("")
    print_2d_array(upper)
    print("")
    print("backtrack".center(50)+ "\n")
    print("")
    print_2d_array(backtrack)
    print("")

def solve(firststring,secondstring,sigma,epsilon,s):
    
    #On top is second word on left is first word
    
    rows = len(firststring)+1
    cols = len(secondstring)+1
    
    lower = init_2d_array(rows,cols)
    middle = init_2d_array(rows,cols)
    upper = init_2d_array(rows,cols)
    backtrack = init_2d_array(rows,cols)

    for i in range(1,rows):
        lower[i][0] = -sigma - (i-1)*epsilon
        middle[i][0] = -sigma - (i-1)*epsilon
        upper[i][0] = -5000

    for j in range(1,cols):
        lower[0][j] = -5000
        middle[0][j] = -sigma - (j-1)*epsilon
        upper[0][j] = -sigma - (j-1)*epsilon
    
    for i in range(1,rows):
        for j in range(1,cols):
            ######################## Filling lower ################################
            from_lower = lower[i-1][j] - epsilon
            from_middle_vertical = middle[i-1][j] - sigma
            lower_ij = max(from_lower,from_middle_vertical)
            lower[i][j] = lower_ij
            ######################## Filling upper ################################
            from_upper = upper[i][j-1] - epsilon
            from_middle_horizontal = middle[i][j-1] - sigma
            upper_ij = max(from_upper,from_middle_horizontal)
            upper[i][j] = upper_ij
            ######################## Filling middle ################################
            char1 = firststring[i-1]
            char2 = secondstring[j-1]
            val = getval(char1 , char2)
            from_middle_diagonal = middle[i-1][j-1] + val
            scores_middle = (lower_ij,from_middle_diagonal,upper_ij)
            middle_ij = max(scores_middle)
            middle[i][j] = middle_ij
            backtrack[i][j] = scores_middle.index(middle_ij) + 1
    

    i = len(firststring)
    j = len(secondstring)
    max_score = middle[i][j] 
        
    output1 , output2 = firststring,secondstring 
    while i*j != 0:
        if backtrack[i][j] == 1:
            i -= 1
            output2 = add_indel(output2,j)
        elif backtrack[i][j] == 3:
            j -= 1
            output1 = add_indel(output1,i)
        else:
            i -= 1
            j -= 1
    
    for repeat in range(i):
        output2 = add_indel(output2,0)
    for repeat in range(j):
        output1 = (output1,0)
    
    answer = (max_score , output1 , output2)

    if s == True:
        showinfo(lower , middle , upper , backtrack)
        print("")

    return answer

string1 = "KWRWPNPCAWDDWALWTRVHTLDFTRYRQMYKVWSSMDHATAGFGCWCWSKHCCIYYMFQKNEVYYAYTWCFHDWVEHMK"
string2 = "PNPYQLWTRVHTLGFTRSCRQMNLKVWSSAGFLCWASSQWSKHCYYMFQKNEVYYAYTWCFHYWTEHMK"
sigma = 11
epsilon = 1
show = False # make false to remove the tables , make True to show tables
answer = solve(string1,string2, sigma , epsilon , show)
print(f"{answer[0]}\n{answer[1]}\n{answer[2]}")












