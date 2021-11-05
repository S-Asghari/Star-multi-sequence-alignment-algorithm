match = 1
mismatch = -1
openGap = -10
gapExtension = -1


def pairwise_alignment(s1, s2):

    newS1 = []
    newS2 = []

    # s1[i] is aligned to s2[j]     #a[1][j][i] is used for traceback
    a = [[[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)] for z in range(2)]

    # s1[i] is aligned to a gap         #b[1][j][i] is used for traceback
    b = [[[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)] for z in range(2)]

    # s2[j] is aligned to a gap         #c[1][j][i] is used for traceback
    c = [[[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)] for z in range(2)]

    # -----------------------------------------------------
    # initialize a, b, c

    for i in range(len(s1)+1):
        a[0][0][i] = -1 * float('inf')
        b[0][0][i] = openGap + gapExtension*(i - 1)
        c[0][0][i] = -1 * float('inf')

    for j in range(len(s2)+1):
        a[0][j][0] = -1 * float('inf')
        b[0][j][0] = -1 * float('inf')
        c[0][j][0] = openGap + gapExtension*(j - 1)

    a[0][0][0] = 0
    b[0][0][0] = openGap - gapExtension
    c[0][0][0] = openGap - gapExtension

    # -----------------------------------------------------
    # fill a, b, c matrices

    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            if s1[i-1] == s2[j-1]:
                ss = match  # ss = s(xi, yj)
            else:
                ss = mismatch

            L1 = [a[0][j-1][i-1]+ss, b[0][j-1][i-1]+ss, c[0][j-1][i-1]+ss]
            a[0][j][i] = max(L1)
            a[1][j][i] = L1.index(max(L1))

            L2 = [a[0][j][i-1]+openGap, b[0][j][i-1]+gapExtension]
            b[0][j][i] = max(L2)
            b[1][j][i] = L2.index(max(L2))

            L3 = [a[0][j-1][i]+openGap, c[0][j-1][i]+gapExtension]
            c[0][j][i] = max(L3)
            c[1][j][i] = L3.index(max(L3))

    # -----------------------------------------------------
    # traceback

    def find_traceback(x, i, j):    # x can be a or b or c

        if i == 0 and j == 0:
            return None     # stop the process

        if x == 'a':

            newS1.append(s1[i - 1])
            newS2.append(s2[j - 1])

            if a[1][j][i] == 0:
                find_traceback('a', i-1, j-1)
            elif a[1][j][i] == 1:
                find_traceback('b', i-1,  j-1)
            else:
                find_traceback('c', i-1, j-1)

        elif x == 'b':

            newS1.append(s1[i - 1])
            newS2.append('-')

            if b[1][j][i] == 0:
                find_traceback('a', i-1, j)
            else:
                find_traceback('b', i-1, j)

        else:

            newS1.append('-')
            newS2.append(s2[j - 1])

            if c[1][j][i] == 0:
                find_traceback('a', i, j-1)
            else:
                find_traceback('c', i, j-1)

    ii = len(s1)
    jj = len(s2)
    L = [a[0][jj][ii], b[0][jj][ii], c[0][jj][ii]]

    if max(L) == a[0][jj][ii]:
        find_traceback('a', ii, jj)
    elif max(L) == b[0][jj][ii]:
        find_traceback('b', ii, jj)
    else:
        find_traceback('c', ii, jj)

    return max(L), newS1, newS2


# -----------------------------------------------------
# initialization

s = input()
seq = s.split()     # seq[0], seq[1], seq[2], seq[3]

score_matrix = [[0 for x in range(4)] for y in range(4)]

# -----------------------------------------------------
# filling score_matrix

for i in range(4):
    for j in range(4):
        if i < j:
            score, newS1, newS2 = pairwise_alignment(seq[i], seq[j])
            score_matrix[j][i] = score
            score_matrix[i][j] = score

# -----------------------------------------------------
# finding Sc

sum_of_scores = [0 for x in range(4)]
for i in range(4):
    for j in range(4):
        sum_of_scores[i] = sum_of_scores[i] + score_matrix[i][j]

Sc = sum_of_scores.index(max(sum_of_scores))

# -----------------------------------------------------
# finding decreasing order of similarity to center Sc

order_of_alignments = sorted(range(len(score_matrix[Sc])), key=lambda k: score_matrix[Sc][k])

# -----------------------------------------------------
# aligning process

aligned_seqs = []

while len(order_of_alignments) > 0:

    nextt = order_of_alignments.pop()
    if nextt == Sc:
        nextt = order_of_alignments.pop()
    aligned_seqs.append(nextt)
    score, newS1, newS2 = pairwise_alignment(seq[Sc], seq[nextt])
    seq[nextt] = list(reversed(newS2))
    newS1 = list(reversed(newS1))

    i2 = 0

    for i in range(len(newS1)):
        if (i2 == len(seq[Sc])) or (newS1[i] != seq[Sc][i2]):
            for j in range(len(aligned_seqs)-1):
                seq[aligned_seqs[j]] = seq[aligned_seqs[j]][:i2] + ['-'] + seq[aligned_seqs[j]][i2:]
        else:
            i2 = i2 + 1

    seq[Sc] = newS1

# -----------------------------------------------------
# print the output

for i in range(4):
    print('new seq' + str(i+1) + ' = "' + ''.join(seq[i]) + '"')

# sample input: ATTGCCATT ATGGCCATT ATCCAATTTT ATCTTCTT
