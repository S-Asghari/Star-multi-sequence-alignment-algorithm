match = 1
mismatch = -1
openGap = -10
gapExtension = -1


def first_gap_position(s):
    for i in range(len(s)):
        if s[i] == '-':
            return i
    return len(s)


def pairwise_alignment(s1, s2):

    # s1[i] is aligned to s2[j]
    a = [[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    # s1[i] is aligned to a gap
    b = [[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    # s2[j] is aligned to a gap
    c = [[0 for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    # trace back for matrix a
    a_tr = [[[] for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    # trace back for matrix b
    b_tr = [[[] for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    # trace back for matrix c
    c_tr = [[[] for x in range(len(s1)+1)] for y in range(len(s2)+1)]

    resList1 = []
    resList2 = []

    # -----------------------------------------------------
    # initialize a, b, c

    for i in range(len(s1)+1):
        a[0][i] = -1 * float('inf')
        b[0][i] = openGap + gapExtension*(i - 1)
        c[0][i] = -1 * float('inf')

    for j in range(len(s2)+1):
        a[j][0] = -1 * float('inf')
        b[j][0] = -1 * float('inf')
        c[j][0] = openGap + gapExtension*(j - 1)

    a[0][0] = 0
    b[0][0] = openGap - gapExtension
    c[0][0] = openGap - gapExtension

    # -----------------------------------------------------
    # fill a, b, c matrices

    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            if s1[i-1] == s2[j-1]:
                ss = match  # ss = s(xi, yj)
            else:
                ss = mismatch

            L1 = [a[j - 1][i - 1] + ss, b[j - 1][i - 1] + ss, c[j - 1][i - 1] + ss]
            a[j][i] = max(L1)
            if a[j - 1][i - 1] + ss == max(L1):
                a_tr[j][i].append('a')
            if b[j - 1][i - 1] + ss == max(L1):
                a_tr[j][i].append('b')
            if c[j - 1][i - 1] + ss == max(L1):
                a_tr[j][i].append('c')

            L2 = [a[j][i - 1] + openGap, b[j][i - 1] + gapExtension]
            b[j][i] = max(L2)
            if a[j][i - 1] + openGap == max(L2):
                b_tr[j][i].append('a')
            if b[j][i - 1] + gapExtension == max(L2):
                b_tr[j][i].append('b')

            L3 = [a[j - 1][i] + openGap, c[j - 1][i] + gapExtension]
            c[j][i] = max(L3)
            if a[j - 1][i] + openGap == max(L3):
                c_tr[j][i].append('a')
            if c[j - 1][i] + gapExtension == max(L3):
                c_tr[j][i].append('c')

    # -----------------------------------------------------
    # traceback

    def find_traceback(x, i, j, newS1, newS2):    # x can be a or b or c

        if i == 0 and j == 0:
            resList1.append(newS1)
            resList2.append(newS2)
            return None     # stop the process

        if x == 'a' and i > 0 and j > 0:

            newS1.append(s1[i - 1])
            newS2.append(s2[j - 1])

            for k in range(len(a_tr[j][i])):
                find_traceback(a_tr[j][i][k], i-1, j-1, newS1.copy(), newS2.copy())

        elif x == 'b' or j == 0:

            newS1.append(s1[i - 1])
            newS2.append('-')

            for k in range(len(b_tr[j][i])):
                find_traceback(b_tr[j][i][k], i-1, j, newS1.copy(), newS2.copy())

        else:

            newS1.append('-')
            newS2.append(s2[j - 1])

            for k in range(len(c_tr[j][i])):
                find_traceback(c_tr[j][i][k], i, j-1, newS1.copy(), newS2.copy())

    ii = len(s1)
    jj = len(s2)
    L = [a[jj][ii], b[jj][ii], c[jj][ii]]

    if max(L) == a[jj][ii]:
        find_traceback('a', ii, jj, [], [])
    if max(L) == b[jj][ii]:
        find_traceback('b', ii, jj, [], [])
    if max(L) == c[jj][ii]:
        find_traceback('c', ii, jj, [], [])

    maxi = -1
    maxi_index = 0

    for n in range(len(resList1)):
        m = first_gap_position(list(reversed(resList1[n]))) + first_gap_position(list(reversed(resList2[n])))
        if m > maxi:
            maxi = m
            maxi_index = n

    return max(L), list(reversed(resList1[maxi_index])),  list(reversed(resList2[maxi_index]))


# -----------------------------------------------------
# initialization

seq = [None for x in range(4)]

for i in range(4):
    seq[i] = input()

score_matrix = [[0 for x in range(4)] for y in range(4)]

resSeq = [[None for x in range(4)] for y in range(4)]

# -----------------------------------------------------
# filling score_matrix

for i in range(4):
    for j in range(4):
        if i < j:
            score, newS1, newS2 = pairwise_alignment(seq[i], seq[j])
            score_matrix[j][i] = score
            score_matrix[i][j] = score
            resSeq[j][i] = newS1.copy()     # new seq[i] after alignment with seq[j]
            resSeq[i][j] = newS2.copy()     # new seq[j] after alignment with seq[i]

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

    i = 0
    i2 = 0

    while i < len(resSeq[nextt][Sc]):

        if (i2 == len(seq[Sc])) or (resSeq[nextt][Sc][i] != seq[Sc][i2]):

            if i2 == len(seq[Sc]):
                for j in range(len(aligned_seqs)-1):
                    resSeq[Sc][aligned_seqs[j]] = resSeq[Sc][aligned_seqs[j]][:i] + ['-'] + resSeq[Sc][aligned_seqs[j]][i:]

            elif resSeq[nextt][Sc][i] == '-':
                for j in range(len(aligned_seqs)-1):
                    resSeq[Sc][aligned_seqs[j]] = resSeq[Sc][aligned_seqs[j]][:i2] + ['-'] + resSeq[Sc][aligned_seqs[j]][i2:]

            elif seq[Sc][i2] == '-':
                resSeq[nextt][Sc] = resSeq[nextt][Sc][:i] + ['-'] + resSeq[nextt][Sc][i:]
                resSeq[Sc][nextt] = resSeq[Sc][nextt][:i] + ['-'] + resSeq[Sc][nextt][i:]
                i2 = i2 + 1

        else:
            i2 = i2 + 1

        if i == len(resSeq[nextt][Sc])-1 and i2 != len(seq[Sc]):
            resSeq[nextt][Sc] = resSeq[nextt][Sc][:i+1] + ['-'] + resSeq[nextt][Sc][i+1:]
            resSeq[Sc][nextt] = resSeq[Sc][nextt][:i+1] + ['-'] + resSeq[Sc][nextt][i+1:]
            i = i + 1
            i2 = i2 + 1

        i = i + 1

    seq[Sc] = resSeq[nextt][Sc].copy()

# -----------------------------------------------------
# print the output

for i in range(4):
    if i==Sc:
        print(''.join(seq[Sc]))
    else:
        print(''.join(resSeq[Sc][i]))

