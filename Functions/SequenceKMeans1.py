import numpy as np
import time
from random import *
from .SequenceFunctions import SequenceFunctions
from .SmithWaterman import SmithWaterman
import matplotlib.pyplot as plt
#adding intracohesion broke code
class SequenceKMeans():
    '''
    Initialize centroids
    sequences_input: List of sequences (eg. [['A','T','G'],['T','G','C']])
    number_of_means: Number of centroid sequences (K) (eg. 2)
    '''
    def initialize(sequences_input, number_of_means):
        global sequences, centroids, intracohesion, aligned_centroids
        sequences = sequences_input
        shuffle(SequenceFunctions.randomize_sequence_batch(sequences))
        centroids = {}
        intracohesion = {}
        aligned_centroids = {}
        for k in range(number_of_means):
            key = ''.join(sequences.pop(k))
            centroids[key] = []
            intracohesion[key] = []
            aligned_centroids[key] = []


    def group():
        # Put sequences to the closest centroid.
        global sequences, centroids, intracohesion, aligned_centroids
        for sequence in sequences:
            best_centroid = None
            best_similarity_value = 0
            for centroid in centroids.keys():
                #print(sequence, centroid)
                similarity, aligned_seq1, aligned_seq2 = SmithWaterman.main(sequence, [x for x in centroid])
                #print(sequence, centroid)
                #print(similarity)
                if similarity > best_similarity_value:
                    best_centroid = centroid
                    best_similarity_value = similarity
            #print(centroids)
            centroids[best_centroid].append(sequence)
            #print(intracohesion)
            intracohesion[best_centroid].append(best_similarity_value)
            aligned_centroids[best_centroid].append(aligned_seq2)

    def calculate_new_centroid():
        global centroids, intracohesion, aligned_centroids

        keys = list(centroids.keys())
        #print(keys)
        for centroid in keys:
            #print(centroid)
            sequences = centroids.pop(centroid)
            intracohesion.pop(centroid)
            aligned_sequences = aligned_centroids.pop(centroid)
            # Determine length of smallest sequence
            if len([len(sequence) for sequence in sequences]) != 0:
                minimum_length = max([len(sequence) for sequence in sequences])
                print(minimum_length)
                new_sequence_letter_count = [{'A':0, 'G':0, 'T':0, 'C':0, 'N':0, 'a':0, 'b':0, 'c':0, 'd':0, 'e':0, 'f':0, 'g':0, 'h':0, 'i':0, 'j':0, 'k':0, '-':0 } for i in range(minimum_length)]
                for sequence in sequences:
                    for index in range(minimum_length):
                        try:
                            letter = sequence[index]
                            new_sequence_letter_count[index][letter] += 1
                        except:
                            p =0
                new_sequence = []
                for index in range(minimum_length):
                    best_key = None
                    best_value = 0
                    for key, value in new_sequence_letter_count[index].items():
                        if value > best_value:
                            if key != '-':
                                best_value = value
                                best_key = key
                        #print(best_value)

                    if best_value != 0: #add condition for if best_value is equal to zero, bc we dont always want A.
                        new_sequence.append(best_key)
                print(new_sequence)
                centroids[''.join(new_sequence)] = []
                intracohesion[''.join(new_sequence)] = []
                aligned_centroids[''.join(new_sequence)] = []
            else:
                #print("What")
                time.sleep(1)
                centroids[centroid] = []
                intracohesion[centroid] = []
                aligned_centroids[centroid] = []
        #print(centroids.keys())

    def print_stats(print_keys = False):
        global centroids, intracohesion
        cohesion_list = []
        for centroid in centroids.keys():
            print('Representative Centroid Sequence: \n', centroid, '\n')
            print('Number of sequences in group: ', str(len(centroids[centroid])), '\n')
            print('Cohesion: ', str(np.mean(intracohesion[centroid])))
            if np.mean(intracohesion[centroid]) >= 0 :
                cohesion_list.append(np.mean(intracohesion[centroid]))
        print('\n')
        return np.mean(cohesion_list)

    def main(sequences_input, num_means, num_iterations):
        SequenceKMeans.initialize(sequences_input, num_means)
        cohesion_overall_list = []
        for i in range(num_iterations):
            print('Iteration: ' + str(i))
            SequenceKMeans.group()
            cohesion_overall_mean = SequenceKMeans.print_stats()
            cohesion_overall_list.append(cohesion_overall_mean)
            SequenceKMeans.calculate_new_centroid()

        plt.figure(1)
        plt.title('Momentum Strategy')
        ax1 = plt.subplot(111)
        ax1.set_title('PORTFOLIO VALUE')
        print(cohesion_overall_list)
        ax1.plot(range(len(cohesion_overall_list)),cohesion_overall_list, color='green', marker='o')
        plt.show()
