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
        global sequences, centroids, intracohesion, intercohesion, aligned_centroids
        sequences = sequences_input
        #s might have repeats
        s = SequenceFunctions.randomize_sequence_batch(sequences)
        shuffle(s)
        print(sequences)
        centroids = {}
        intracohesion = {}
        aligned_centroids = {}
        intercohesion = {}
        for k in range(number_of_means):
            key = ''.join(s.pop(k))
            while key in centroids.keys():
                y = SequenceFunctions.randomize_sequence_batch(sequences)
                shuffle(y)
                key = ''.join(y.pop(0))
            centroids[key] = []
            intracohesion[key] = []
            intercohesion[key] = []
            aligned_centroids[key] = []
        centroids['misc'] = []
        intracohesion['misc'] = []
        aligned_centroids['misc'] = []
        intercohesion['misc'] = []
        print('init', str(len(sequences)))


    def group():
        # Put sequences to the closest centroid.
        global sequences, centroids, intracohesion, aligned_centroids, intercohesion
        #print("Num Sequences: " + str(len(sequences)))
        for sequence in sequences:
            best_centroid = 'misc'
            first_best_similarity_value = 0
            second_best_similarity_value = 0
            for centroid in centroids.keys():
                #print(sequence, centroid)
                similarity, consensus_position = SmithWaterman.main(sequence, [x for x in centroid])
                centroids[centroid].append((sequence, similarity))
                #print(sequence, centroid)
                #print(similarity)

                if similarity > first_best_similarity_value:
                    if first_best_similarity_value != 0:
                        second_best_similarity_value = first_best_similarity_value
                    best_centroid = centroid
                    first_best_similarity_value = similarity
                elif similarity > second_best_similarity_value:
                    second_best_similarity_value = similarity

            '''
            mean = np.mean([item[0] for item in similarity_arr])
            std = np.std([item[0] for item in similarity_arr])
            norm_similarity_arr = [((item[0]-mean)/std, item[1]) for item in similarity_arr]
            norm_min = min([item[0] for item in norm_similarity_arr])
            norm_max = max([item[0] for item in norm_similarity_arr])
            scaled_similarity_arr = [((float(item[0]) - norm_min)/(norm_max - norm_min), item[1]) for item in norm_similarity_arr]
            sorted_centroid_weights = sor
            '''

            #print(centroids)
            #centroids[best_centroid].append(sequence)
            #print(intracohesion)
            intracohesion[best_centroid].append(first_best_similarity_value)
            intercohesion[best_centroid].append(second_best_similarity_value)
            aligned_centroids[best_centroid].append(consensus_position)

    def calculate_new_centroid():
        global centroids, intracohesion, aligned_centroids, intercohesion, sequences

        keys = list(centroids.keys())
        #print(keys)
        for centroid in keys:
            #print(centroid)
            consensus_sequences = centroids.pop(centroid)
            intracohesion.pop(centroid)
            consensus_positions = aligned_centroids.pop(centroid)
            intercohesion.pop(centroid)
            #print(consensus_positions)
            # Determine length of smallest sequence
            if len(consensus_sequences) != 0:
                minimum_length = min([len(sequence) for sequence in consensus_sequences])
                #print(minimum_length)
                new_sequence_letter_count = [{'A':0, 'G':0, 'T':0, 'C':0, 'N':0, 'a':0, 'b':0, 'c':0, 'd':0, 'e':0, 'f':0, 'g':0, 'h':0, 'i':0, 'j':0, 'k':0, '-':0 } for i in range(minimum_length)]
                consensus_list = []
                for list_index in range(len(consensus_sequences)):
                    sequence = consensus_sequences[list_index]
                    consensus_start_position, consensus_end_position = consensus_positions[list_index]
                    consensus_list.append(consensus_end_position - consensus_start_position)
                    for index in range(minimum_length):
                        try:
                            #if index >= consensus_start_position and index < consensus_end_position:
                            letter = sequence[index]
                            new_sequence_letter_count[index][letter] += 1
                            #else:
                            #new_sequence_letter_count[index]['-'] -= 1
                        except:
                            pass
                #print(np.mean(consensus_list))
                new_sequence = []
                for index in range(minimum_length):
                    best_key = None
                    best_value = 0
                    for key, value in new_sequence_letter_count[index].items():
                        if value > best_value:
                            best_value = value
                            best_key = key
                        #print(best_value)

                    if best_value != 0 and best_key != '-': #add condition for if best_value is equal to zero, bc we dont always want A.
                        new_sequence.append(best_key)
                #print(''.join(new_sequence))
                new_centroid = ''.join(new_sequence)
                if 'misc' == centroid:
                    new_centroid = 'misc'
                while new_centroid in centroids.keys():
                    s = SequenceFunctions.randomize_sequence_batch(sequences)
                    shuffle(s)
                    new_centroid = ''.join(s.pop(0))
                centroids[new_centroid] = []
                intracohesion[new_centroid] = []
                aligned_centroids[new_centroid] = []
                intercohesion[new_centroid] = []
            else:
                #print("What")
                #time.sleep(1)
                s = SequenceFunctions.randomize_sequence_batch(sequences)
                shuffle(s)
                new_centroid = ''.join(s.pop(0))
                if 'misc' == centroid:
                    new_centroid = 'misc'
                while new_centroid in centroids.keys():
                    s = SequenceFunctions.randomize_sequence_batch(sequences)
                    shuffle(s)
                    new_centroid = ''.join(s.pop(0))
                centroids[new_centroid] = []
                intracohesion[new_centroid] = []
                aligned_centroids[new_centroid] = []
                intercohesion[new_centroid] = []

        centroids['misc'] = []
        intracohesion['misc'] = []
        aligned_centroids['misc'] = []
        intercohesion['misc'] = []

        #print(centroids.keys())

    def print_stats(print_keys = False):
        global sequences, centroids, intracohesion, aligned_centroids, intercohesion
        cohesion_list = []
        intercohesion_list = []
        print('________________________________________________________')
        print("Num Groups: " + str(len(centroids.keys())))
        for centroid in centroids.keys():
            #print('Representative Centroid Sequence: \n', centroid, '\n')
            #print('Number of sequences in group: ', str(len(centroids[centroid])), '\n')

            #print('Cohesion: ', str(np.mean(intracohesion[centroid])))
            if np.mean(intracohesion[centroid]) >= 0 :
                cohesion_list.append(np.mean(intracohesion[centroid]) * len(intracohesion[centroid])/len(sequences))
            #print('Intercohesion: ', str(np.mean(intercohesion[centroid])))
            if np.mean(intercohesion[centroid]) >= 0:
                intercohesion_list.append(np.mean(intercohesion[centroid]) * len(intercohesion[centroid])/len(sequences))

            if len(aligned_centroids[centroid]):
                start = [position[0] for position in aligned_centroids[centroid]]
                end = [position[1] for position in aligned_centroids[centroid]]
                #print(aligned_centroids[centroid])
                #print('Start:' +str(np.mean(start)) +"-"+ str(np.std(start)))
                #print('End:' +str(np.mean(end)) +"+"+ str(np.std(end)))
        print('\n')
        return np.sum(cohesion_list), np.sum(intercohesion_list), centroids.keys()

    def main(sequences_input, num_means, num_iterations):
        SequenceKMeans.initialize(sequences_input, num_means)
        cohesion_overall_list = []
        intercohesion_overall_list = []
        keys = []
        SequenceKMeans.group()
        cohesion_overall_mean, intercohesion_overall_mean, keys = SequenceKMeans.print_stats()
        cohesion_overall_list.append(cohesion_overall_mean)
        intercohesion_overall_list.append(intercohesion_overall_mean)
        for i in range(num_iterations):
            print('Iteration: ' + str(i))
            SequenceKMeans.calculate_new_centroid()
            SequenceKMeans.group()
            cohesion_overall_mean, intercohesion_overall_mean, keys = SequenceKMeans.print_stats()
            cohesion_overall_list.append(cohesion_overall_mean)
            intercohesion_overall_list.append(intercohesion_overall_mean)

        '''
        plt.figure(1)
        plt.title('Inter/Intra-Cohesion Values')
        ax1 = plt.subplot(111)
        ax1.set_title('Group Cohesion over Iterations')
        print(cohesion_overall_list)
        print(intercohesion_overall_list)
        ax1.plot(range(len(cohesion_overall_list)),cohesion_overall_list, color='green', marker='o')
        ax1.plot(range(len(intercohesion_overall_list)),intercohesion_overall_list, color='red', marker='o')
        #plt.show()
        '''
        return cohesion_overall_list, intercohesion_overall_list, keys
