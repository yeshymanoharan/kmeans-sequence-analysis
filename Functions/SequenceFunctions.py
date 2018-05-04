import random
import itertools

'''
Provides utility functions for sequence manipulation
'''
class SequenceFunctions():

    '''
    Splits nucleotide sequence based on recurring pattern found in feature

    ex.
    input = AAGGTC,
            SSSSBB
    > AAGG, TC

    '''
    def separate_sequence():
        print("Not Implemented Yet")

    '''
    Function creates random sequence based on given sequence.
    "combine_repetitions" == True: Combines repetitions
    Ex. AATGTC -> AA T G T C -> T T AA C G -> TTAACG
    '''
    def randomize_sequence(sequence, combine_repetitions=True):
        new_sequence = []
        previous_char_array = []
        previous_char_array.append(sequence[0])
        for index in range(1, len(sequence)):
            if combine_repetitions == True:
                if sequence[index] == previous_char_array[0]:
                    previous_char_array.append(sequence[index])
                else:
                    new_sequence.append(previous_char_array)
                    previous_char_array = []
                    previous_char_array.append(sequence[index])
            else:
                new_sequence.append(previous_char_array)
                previous_char_array = []
                previous_char_array.append(sequence[index])
        new_sequence.append(previous_char_array)
        random.shuffle(new_sequence)
        return list(itertools.chain.from_iterable(new_sequence))

    def randomize_sequence_batch(sequences, combine_repetitions=True):
        sequences_array = []
        for sequence in sequences:
            sequences_array.append(''.join(SequenceFunctions.randomize_sequence(sequence, combine_repetitions)))
        return sequences_array

    def strip_head_tail(sequences, head=None, tail=None):
        stripped_sequences = []
        for sequence in sequences:
            #sequence = list(sequence[0])
            if head != None:
                head = list(head)
                #print("Removing Head.")
                #print(len(sequence))
                #print(sequence)
                index=0
                while len(head) > index and head[index] == sequence[0]:
                    sequence.pop(0)
                    index += 1
                #print(sequence)


            if tail != None:
                tail = list(tail)
                index=len(tail) - 1
                while index >= 0 and tail[index] == sequence[-1]:
                    sequence.pop(-1)
                    index -= 1
                #print(sequence)
                #print("Removing Tail.")
            #print(sequence)
            stripped_sequences.append(sequence)
        return stripped_sequences

    def generate_random_sequences(patterns, weights, length_of_sequence, number_of_sequences):

        def random_character():
            num = random.randint(0,3)
            if num == 0:
                return 'A'
            elif num == 1:
                return 'G'
            elif num == 2:
                return 'T'
            else:
                return 'C'

        generated_sequences = []
        generated_names = []
        for i in range(number_of_sequences):
            sequence = ""
            num = random.random()
            weight_sum = weights[0]
            pattern_index = 0
            while num > weight_sum:
                pattern_index += 1
                weight_sum += weights[pattern_index]
            pattern_starting_index_limit = length_of_sequence - len(patterns[pattern_index])
            num = int(random.random() * pattern_starting_index_limit)
            for i in range(num):
                sequence+= random_character()
            sequence+= patterns[pattern_index]
            for i in range(pattern_starting_index_limit-num):
                sequence+= random_character()
            generated_sequences.append(sequence)
            generated_names.append(">testseq" + str(i))

        return generated_sequences, generated_names


if __name__ == "__main__":
    sequence = SequenceFunctions.randomize_sequence("AATGTC")
    print(sequence)
