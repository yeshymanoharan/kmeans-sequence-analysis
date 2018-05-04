'''
A Python implementation of the Smith-Waterman algorithm for local alignment
of nucleotide sequences.
'''

import argparse
import os
import re
import sys
import unittest

class SmithWaterman():
    # These scores are taken from Wikipedia.
    # en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    match    = 10
    mismatch = -5
    gap      = -10
    seq1     = None
    seq2     = None
    debug = True

    def main(sequence1, sequence2, debug=debug):
        SmithWaterman.debug = debug
        # Add static sequence variables.
        SmithWaterman.seq1 = sequence1
        SmithWaterman.seq2 = sequence2

        # The scoring matrix contains an extra row and column for the gap (-), hence
        # the +1 here.
        rows = len(SmithWaterman.seq1) + 1
        cols = len(SmithWaterman.seq2) + 1

        # Initialize the scoring matrix.
        score_matrix, start_pos, no_alignment = SmithWaterman.create_score_matrix(rows, cols)
        if no_alignment == 1:
            return 0, (0,0)

        # Traceback. Find the optimal path through the scoring matrix. This path
        # corresponds to the optimal local sequence alignment.
        seq1_aligned, seq2_aligned = SmithWaterman.traceback(score_matrix, start_pos)
        #print(seq1_aligned, seq2_aligned)
        assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'
        #print(seq1_aligned, seq2_aligned)
        # Pretty print the results. The printing follows the format of BLAST results
        # as closely as possible.
        alignment_str, idents, gaps, mismatches = SmithWaterman.alignment_string(seq1_aligned, seq2_aligned)
        alength = max(len(sequence1), len(sequence2))
        if (SmithWaterman.debug):
            print()
            print('sequence: ' + str(''.join(sequence1)))
            print('centroid: ' + str(''.join(sequence2)))
            print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
                alength, idents / alength, gaps, alength, gaps / alength))
            print()
            for i in range(0, alength, 100):
                seq1_slice = seq1_aligned[i:i+100]
                print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
                print('             {0}'.format(alignment_str[i:i+100]))
                seq2_slice = seq2_aligned[i:i+100]
                print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
                print()
        string = seq2_aligned.replace('-','')
        start_position=start_pos[1] - len(string)
        end_position=start_pos[1]
        return idents / alength * 1.0, (start_position,end_position)


    def create_score_matrix(rows, cols):
        '''Create a matrix of scores representing trial alignments of the two sequences.
        Sequence alignment can be treated as a graph search problem. This function
        creates a graph (2D matrix) of scores, which are based on trial alignments
        of different base pairs. The path with the highest cummulative score is the
        best alignment.
        '''
        score_matrix = [[0 for col in range(cols)] for row in range(rows)]

        # Fill the scoring matrix.
        max_score = 0
        max_pos   = None    # The row and columbn of the highest score in matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                score = SmithWaterman.calc_score(score_matrix, i, j)
                if score > max_score:
                    max_score = score
                    max_pos   = (i, j)

                score_matrix[i][j] = score

        if max_pos is None:
            #'the x, y position with the highest score was not found'
            return 0, 0, 1

        return score_matrix, max_pos, 0


    def calc_score(matrix, x, y):
        '''Calculate score for a given x, y position in the scoring matrix.
        The score is based on the up, left, and upper-left neighbors.
        '''
        similarity = SmithWaterman.match if SmithWaterman.seq1[x - 1] == SmithWaterman.seq2[y - 1] else SmithWaterman.mismatch

        diag_score = matrix[x - 1][y - 1] + similarity
        up_score   = matrix[x - 1][y] + SmithWaterman.gap
        left_score = matrix[x][y - 1] + SmithWaterman.gap

        return max(0, diag_score, up_score, left_score)


    def traceback(score_matrix, start_pos):
        '''Find the optimal path through the matrix.
        This function traces a path from the bottom-right to the top-left corner of
        the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
        or both of the sequences being aligned. Moves are determined by the score of
        three adjacent squares: the upper square, the left square, and the diagonal
        upper-left square.
        WHAT EACH MOVE REPRESENTS
            diagonal: match/mismatch
            up:       gap in sequence 1
            left:     gap in sequence 2
        '''

        END, DIAG, UP, LEFT = range(4)
        arr = ['end', 'diag', 'up', 'left']
        aligned_seq1 = []
        aligned_seq2 = []
        x, y         = start_pos
        move         = SmithWaterman.next_move(score_matrix, x, y)
        #print(start_pos)
        #print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in score_matrix]))
        while move != END:
            #print(arr[move])
            if move == DIAG:
                aligned_seq1.append(SmithWaterman.seq1[x - 1])
                aligned_seq2.append(SmithWaterman.seq2[y - 1])
                x -= 1
                y -= 1
            elif move == UP:
                aligned_seq1.append(SmithWaterman.seq1[x - 1])
                aligned_seq2.append('-')
                x -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(SmithWaterman.seq2[y - 1])
                y -= 1

            move = SmithWaterman.next_move(score_matrix, x, y)

        aligned_seq1.append(SmithWaterman.seq1[x - 1])
        aligned_seq2.append(SmithWaterman.seq2[y - 1])


        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


    def next_move(score_matrix, x, y):
        diag = score_matrix[x - 1][y - 1]
        up   = score_matrix[x - 1][y]
        left = score_matrix[x][y - 1]
        if diag >= up and diag >= left:     # Tie goes to the DIAG move.
            #print(diag, up, left)
            return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
        elif up > diag and up >= left:      # Tie goes to UP move.
            #print(diag, up, left)
            return 2 if up != 0 else 0      # UP move or end.
        elif left > diag and left > up:
            #print(diag, up, left)
            return 3 if left != 0 else 0    # LEFT move or end.
        else:
            # Execution should not reach here.
            raise ValueError('invalid move during traceback')


    def alignment_string(aligned_seq1, aligned_seq2):
        '''Construct a special string showing identities, gaps, and mismatches.
        This string is printed between the two aligned sequences and shows the
        identities (|), gaps (-), and mismatches (:). As the string is constructed,
        it also counts number of identities, gaps, and mismatches and returns the
        counts along with the alignment string.
        AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
        ::||::::::||:|::::::: |:  :||:|   <-- alignment string
        CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
        '''
        # Build the string as a list of characters to avoid costly string
        # concatenation.
        idents, gaps, mismatches = 0, 0, 0
        alignment_string = []
        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == base2:
                alignment_string.append('|')
                idents += 1
            elif '-' in (base1, base2):
                alignment_string.append(' ')
                gaps += 1
            else:
                alignment_string.append(':')
                mismatches += 1

        return ''.join(alignment_string), idents, gaps, mismatches


    def print_matrix(matrix):
        '''Print the scoring matrix.
        ex:
        0   0   0   0   0   0
        0   2   1   2   1   2
        0   1   1   1   1   1
        0   0   3   2   3   2
        0   2   2   5   4   5
        0   1   4   4   7   6
        '''
        for row in matrix:
            for col in row:
                print('{0:>4}'.format(col))
            print()


class ScoreMatrixTest(unittest.TestCase):
    '''Compare the matrix produced by create_score_matrix() with a known matrix.'''
    def test_matrix(self):
        # From Wikipedia (en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
        #                -   A   C   A   C   A   C   T   A
        known_matrix = [[0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
                        [0,  2,  1,  2,  1,  2,  1,  0,  2],  # A
                        [0,  1,  1,  1,  1,  1,  1,  0,  1],  # G
                        [0,  0,  3,  2,  3,  2,  3,  2,  1],  # C
                        [0,  2,  2,  5,  4,  5,  4,  3,  4],  # A
                        [0,  1,  4,  4,  7,  6,  7,  6,  5],  # C
                        [0,  2,  3,  6,  6,  9,  8,  7,  8],  # A
                        [0,  1,  4,  5,  8,  8, 11, 10,  9],  # C
                        [0,  2,  3,  6,  7, 10, 10, 10, 12]]  # A

        global seq1, seq2
        seq1 = 'AGCACACA'
        seq2 = 'ACACACTA'
        rows = len(seq1) + 1
        cols = len(seq2) + 1

        matrix_to_test, max_pos = SmithWaterman.create_score_matrix(rows, cols)
        self.assertEqual(known_matrix, matrix_to_test)

class SmithWatermanBatch():

    def main(sequences):
        similarity_array = []
        num_sequences = len(sequences)
        for i in range(0, num_sequences):
            print('\r', round(i / num_sequences * 1.0, 2))
            similarity_row = []
            for j in range(i, num_sequences):
                SmithWaterman.debug = False
                similarity_row.append(round(SmithWaterman.main(sequences[i], sequences[j]),3))
            similarity_array.append(similarity_row)
        return similarity_array

if __name__ == '__main__':
    print(SmithWaterman.main('CTTCTACAGTAACTTTATGGGTTAGAATTATATCTTGTGA', 'AAGGATATTTCATAGTAAGTGAATTAAGTGTACATTATTA'))
