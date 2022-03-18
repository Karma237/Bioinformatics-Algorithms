'''
yousra samir elamir
Karma yasser ismaiel
fatma amr biomy
omar mohamed hamed
hasnaa sayed abdo
youssef samir elamir
'''
import numpy as np

choice = input("please choose the molecule type DNA or protein (1/2): ")
if choice == "1":
    def localAligment(X, Y):
        len_x = len(X)  # M
        len_y = len(Y)  # N
        matrix = [[0 for X in range(len_x + 1)] for Y in range(len_y + 1)]
        optimal_point_row = 0
        optimal_point_column = 0
        optimal_point = 0
        for i in range(1, len_y + 1):
            for j in range(1, len_x + 1):
                if X[j - 1] == Y[i - 1]:
                    score = match
                else:
                    score = mismatch

                matrix[i][j] = max(matrix[i - 1][j] + gap,
                                   matrix[i][j - 1] + gap,
                                   matrix[i - 1][j - 1] + score,
                                   0)
                if matrix[i][j] > optimal_point:
                    optimal_point = matrix[i][j]
                    optimal_point_row = i
                    optimal_point_column = j
        print("\nTHE MATRIX :- ")
        print("-----------------")
        print(matrix)

        # TRACING

        LA_seq1 = ""
        LA_seq2 = ""
        LA_matching = ""
        i = optimal_point_row
        j = optimal_point_column

        while (True):
            if (matrix[i][j] == 0):
                break

            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            if (X[j - 1] == Y[i - 1]):
                score = match
            else:
                score = mismatch
            diagonal = matrix[i - 1][j - 1] + score

            if (matrix[i][j] == diagonal):
                LA_seq1 += X[j - 1]
                LA_seq2 += Y[i - 1]
                if score == match:
                    LA_matching += "|"
                else:
                    LA_matching += " "
                i -= 1
                j -= 1
            elif matrix[i][j] == up:
                LA_seq1 += "-"
                LA_seq2 += Y[i - 1]
                LA_matching += " "
                i -= 1
            elif matrix[i][j] == left:
                LA_seq1 += X[j - 1]
                LA_seq2 += "-"
                LA_matching += " "
                j -= 1

        # PRINTING
        LA_seq1 = LA_seq1[::-1]  # to reverse the sequence
        LA_seq2 = LA_seq2[::-1]
        LA_matching = LA_matching[::-1]
        print("\nTHE BEST LOCAL ALIGNMENT IS :")
        print("------------------------------\n")
        print(LA_seq1)
        print(LA_matching)
        print(LA_seq2)
        print("\nTHE SCORE IS :- ", optimal_point)
        print("-----------------")
        print("\nOPTIMAL SOLUTION WAS AT CELL :- ( ", optimal_point_row, " , ", optimal_point_column, " )")
        print("-------------------------------")


    # CALLING FUNCTION
    print("please enter the following info, ")

    X = input("the first DNA sequence :")
    Y = input("the second DNA sequence :")
    match = int(input("match score :"))
    mismatch = int(input("mismatch score :"))
    gap = int(input("Gap score :"))

    localAligment(X, Y)


elif choice == "2":

    blosum62 = {
        ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
        ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
        ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
        ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
        ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
        ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
        ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
        ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
        ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
        ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
        ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
        ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
        ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
        ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
        ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
        ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
        ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
        ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
        ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
        ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
        ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
        ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
        ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
        ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
        ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
        ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
        ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
        ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
        ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
        ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
        ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
        ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
        ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
        ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
        ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
        ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
        ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
        ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
        ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
        ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
        ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
        ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
        ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
        ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
        ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
        ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
        ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
        ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
        ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
        ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
        ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
        ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
        ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
        ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
        ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
        ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
        ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
        ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
        ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
        ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
        ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
        ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
        ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
        ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
        ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
        ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
        ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
        ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
        ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
    }


    def BLOSUM62(letter1, letter2):
        if letter1 == letter2:
            if (letter1, letter2) in blosum62:
                return blosum62[(letter1, letter2)]
            else:
                return blosum62[(letter2, letter1)]
        elif letter1 == '-' or letter2 == '-':
            return Gap
        else:
            if (letter1, letter2) in blosum62:
                return blosum62[(letter1, letter2)]
            else:
                return blosum62[(letter2, letter1)]


    def LocalProteinAlignment(seq1, seq2, blosumMatrix):
        optimal_point_row = 0
        optimal_point_column = 0
        optimal_point = 0
        lenSeq1 = len(seq1)
        lenSeq2 = len(seq2)
        scoreMatrix = np.zeros((lenSeq1 + 1, lenSeq2 + 1))
        # Initialization matrix
        for i in range(lenSeq1 + 1):
            scoreMatrix[i][0] = 0
        for j in range(lenSeq2 + 1):
            scoreMatrix[0][j] = 0

        for i in range(1, lenSeq1 + 1):
            for j in range(1, lenSeq2 + 1):
                diagonal = scoreMatrix[i - 1][j - 1] + int(blosumMatrix(seq1[i - 1], seq2[j - 1]))
                up = scoreMatrix[i - 1][j] + Gap
                left = scoreMatrix[i][j - 1] + Gap
                if max(diagonal, up, left) >= 0:
                    scoreMatrix[i][j] = max(diagonal, up, left)
                else:
                    scoreMatrix[i][j] = 0

                if scoreMatrix[i][j] > optimal_point:
                    optimal_point = scoreMatrix[i][j]
                    optimal_point_row = i
                    optimal_point_column = j
        print(F'score matrix = \n{scoreMatrix}')

        # TRACING

        LA_seq1 = ""
        LA_seq2 = ""
        LA_matching = ""
        i = optimal_point_row
        j = optimal_point_column

        while (True):
            if (scoreMatrix[i][j] == 0):
                break

            up = scoreMatrix[i - 1][j] + Gap
            left = scoreMatrix[i][j - 1] + Gap
            diagonal = scoreMatrix[i - 1][j - 1] + int(blosumMatrix(seq1[i - 1], seq2[j - 1]))

            if (scoreMatrix[i][j] == diagonal):
                LA_seq1 += seq1[i - 1]
                LA_seq2 += seq2[j - 1]
                if seq1[i - 1] == seq2[j - 1]:
                    LA_matching += "|"
                else:
                    LA_matching += " "
                i -= 1
                j -= 1
            elif scoreMatrix[i][j] == up:
                LA_seq1 += seq1[i - 1]
                LA_seq2 += "-"
                LA_matching += " "
                i -= 1
            elif scoreMatrix[i][j] == left:
                LA_seq1 += "-"
                LA_seq2 += seq2[j - 1]
                LA_matching += " "
                j -= 1

        # PRINTING
        LA_seq1 = LA_seq1[::-1]  # to reverse the sequence
        LA_seq2 = LA_seq2[::-1]
        LA_matching = LA_matching[::-1]
        print("\nTHE BEST LOCAL ALIGNMENT IS :")
        print("------------------------------\n")
        print(LA_seq1)
        print(LA_matching)
        print(LA_seq2)
        print("\nTHE SCORE IS :- ", optimal_point)
        print("-----------------")
        print("\nOPTIMAL SOLUTION WAS AT CELL :- ( ", optimal_point_row, " , ", optimal_point_column, " )")
        print("-------------------------------")


    seq1 = input("please enter the first protein sequence: ")
    seq2 = input("please enter the second protein sequence: ")
    Gap = -1
    # LocalProteinAlignment('GHTAYG', 'GHTAG', BLOSUM62)
    LocalProteinAlignment(seq1, seq2, BLOSUM62)

else:
    print("invalid")
