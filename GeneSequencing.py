#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []

            for j in range(len(sequences)):

                if (j < i):
                    s = {}
                else:
                    # Compute the alignment score and actual character-by-character alignment
                    if banded:
                        score, alignment1, alignment2 = self.banded_needleman_wunsch(sequences[i][:align_length],
                                                                                     sequences[j][:align_length],
                                                                                     MAXINDELS)
                    else:
                        score, alignment1, alignment2 = self.unrestricted_needleman_wunsch(sequences[i][:align_length],
                                                                                           sequences[j][:align_length])

                    s = {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
                    table.item(i, j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.update()

                jresults.append(s)

            results.append(jresults)

        return results

    # Implement Needleman-Wunsch algorithm for unrestricted version
    # Time complexity  = O(nm) where n is length of first sequence and m is length of second sequence
    # Space complexity = O(nm)
    def unrestricted_needleman_wunsch(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)

        # Initialize DP table - O(nm)
        dp = [[0] * (m + 1) for i in range(n + 1)]
        for i in range(n + 1):
            dp[i][0] = i * INDEL

        for j in range(m + 1):
            dp[0][j] = j * INDEL

        # Fill in DP table - O(nm)
        for i in range(1, n + 1):

            for j in range(1, m + 1):
                match = dp[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB)
                delete = dp[i - 1][j] + INDEL
                insert = dp[i][j - 1] + INDEL
                dp[i][j] = min(match, delete, insert)

        # Compute actual alignment O(n + m)
        alignment1 = ""
        alignment2 = ""
        i = n
        j = m
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB):
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = seq2[j - 1] + alignment2
                i -= 1
                j -= 1
            elif i > 0 and dp[i][j]:
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = '-' + alignment2
                i -= 1
            else:
                alignment1 = '-' + alignment1
                alignment2 = seq2[j - 1] + alignment2
                j -= 1
        return dp[n][m], alignment1, alignment2

    # Time complexity = O(kn) where k is the max indel, n is the length of first sequence, m is length of second.
    # Space complexity = O(kn)
    def banded_needleman_wunsch(self, seq1, seq2, max_indels):
        n = len(seq1)
        m = len(seq2)

        # Compute bandwidth for banded version - O(1)
        k = 2 * max_indels + 1

        # Initialize DP table - O(kn)
        dp = [[math.inf] * (m + 1 + k) for i in range(n + 1)]

        # Initialize the first row and column - O(1)
        dp[0][0] = 0
        for i in range(1, n + 1):
            dp[i][0] = i * INDEL
        for j in range(1, min(m + 1, k + 1)):
            dp[0][j] = j * INDEL

        # Fill in DP table - O(kn)
        for i in range(1, n + 1):
            for j in range(max(1, i - k), min(m + 1, i + k + 1)):
                match = dp[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB)
                delete = dp[i - 1][j] + INDEL
                insert = dp[i][j - 1] + INDEL
                dp[i][j] = min(match, delete, insert)

        # Compute actual alignment - O(n + m)
        alignment1 = ""
        alignment2 = ""
        i = n
        j = m
        while i > 0 or j > 0:
            if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB):
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = seq2[j - 1] + alignment2
                i -= 1
                j -= 1
            elif i > 0 and dp[i][j] == dp[i - 1][j] + INDEL:
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = '-' + alignment2
                i -= 1
            else:
                alignment1 = '-' + alignment1
                alignment2 = seq2[j - 1] + alignment2
                j -= 1
        return dp[n][m], alignment1, alignment2
