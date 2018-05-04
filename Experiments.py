from Functions.DataReader import DataReader
from Functions.SmithWaterman import SmithWatermanBatch
from Functions.SequenceFunctions import SequenceFunctions
from Functions.SequenceKMeans import SequenceKMeans
from Functions.Utils import mkdir_p, safe_open, getopts
from sys import argv
import matplotlib.pyplot as plt
import random
import numpy as np
import csv
import math
import json

'''
Main code that contains experiments and helper functions to carry out those experiments.
'''

'''
Function that calculates similarity between sequences with respect to primary and secondary structure.
'''
def primary_secondary_similarity():
    print("Reading...")
    ## Read sequences from file
    sequences = DataReader.convert_file_to_array("Data/MKP3PositiveSequences.csv")

    ## Index 2 corresponds to primary structure.
    primary_sequences = DataReader.strip_array_by_indices(sequences,9,[2])

    ## Index 7 corresponds to secondary structure. Needs to be formatted to one character representations.
    secondary_sequences = DataReader.strip_array_by_indices(sequences,9,[7])
    secondary_input_dictionary = {'S':'a','B':'b', 'D':'c', 'DH':'d', 'BH':'e', 'IS':'f', 'IA':'g','ISH':'h', 'IAH':'i', 'M':'j', 'MH':'k'}
    secondary_sequences = DataReader.replace_items_in_array(secondary_sequences, secondary_input_dictionary)

    ## Randomize Sequences
    randomized_primary_sequences = SequenceFunctions.randomize_sequence_batch(primary_sequences)
    print(randomized_primary_sequences)
    randomized_secondary_sequences = SequenceFunctions.randomize_sequence_batch(secondary_sequences)

    print("Calculating...")
    primary_sequences_similarity = SmithWatermanBatch.main(primary_sequences)
    secondary_sequences_similarity = SmithWatermanBatch.main(secondary_sequences)
    randomized_primary_sequences_similarity = SmithWatermanBatch.main(randomized_primary_sequences)
    randomized_secondary_sequences_similarity = SmithWatermanBatch.main(randomized_secondary_sequences)
    combined_array = DataReader.combine_arrays(primary_sequences_similarity, secondary_sequences_similarity, "Primary Similarity", "Secondary Similarity")
    randomized_combined_array = DataReader.combine_arrays(randomized_primary_sequences_similarity, randomized_secondary_sequences_similarity, "Randomized Primary Similarity", "Randomized Secondary Similarity")

    ## Plot
    print("Plotting...")
    X = [float(point[0]) for point in combined_array[1:]]
    Y = [float(point[1]) for point in combined_array[1:]]
    plt.subplot(2,2,1)
    plt.scatter(X,Y, s=1)
    plt.xlabel(combined_array[0][0])
    plt.ylabel(combined_array[0][1])

    X = [float(point[0]) for point in randomized_combined_array[1:]]
    Y = [float(point[1]) for point in randomized_combined_array[1:]]
    plt.subplot(2,2,2)
    plt.scatter(X,Y, s=1)
    plt.xlabel(randomized_combined_array[0][0])
    plt.ylabel(randomized_combined_array[0][1])
    plt.show()

'''
Main K-Means code
'''
def k_means(**kwargs):
    names, sequences = load_sequences(
        filename=kwargs['filename'],
        sequence_type=kwargs['sequence_type'],
        head=kwargs['head'],
        tail=kwargs['tail'],
        num_index_per_group=kwargs['num_index_per_group_in_file'],
        num_index_primary=kwargs['num_index_primary_in_file'],
        num_index_secondary=kwargs['num_index_secondary_in_file'],
        select_top_sequence=kwargs['select_top_sequence']
    )
    if kwargs['run_k_means']:
        k_means_run(
            filename=kwargs['filename'],
            sequence_type=kwargs['sequence_type'],
            generated_sequences=sequences,
            lower_mean=kwargs['lower_mean'],
            upper_mean=kwargs['upper_mean'],
            num_algo_restarts=kwargs['num_algo_restarts'],
            num_algo_iterations=kwargs['num_algo_iterations'],
            plot=kwargs['plot']
        )
    else:
        result = analyse_k_means(
            filename=kwargs['filename'],
            sequence_type=kwargs['sequence_type'],
            sequences=sequences,
            names=names,
            plot=kwargs['plot'],
            upper_mean=kwargs['upper_mean'],
            lower_mean=kwargs['lower_mean']
        )
        return result
        #k_means_test(filename=filename, sequence_type=sequence_type, sequences=sequences,names=names, centroids=centroids)

def analyse_k_means(filename, sequence_type, sequences, names, plot, upper_mean, lower_mean):
    k_means_data = {
        "intracohesion_min":[],
        "intracohesion_max":[],
        "intercohesion_min":[],
        "intercohesion_max":[],
        "centroids_min":[],
        "centroids_max":[]
    }
    try:
        for j in range(lower_mean,upper_mean):
            with safe_open("Data/k_means_data/" + filename.split('/')[-1][:-4] + "/" + sequence_type + "/means" + str(j) + ".csv", "r") as f:
                reader = csv.reader(f)
                key = None
                for row in reader:
                    #print(key)

                    if key == None:
                        key = row[0]
                    else:
                        #print(k_means_data[key])
                        if key == "centroids_min" or key == "centroids_max":
                            k_means_data[key].append(row)
                        else:
                            k_means_data[key].append(row[-1])
                        key = None

    except Exception as e:
        print(e)
        pass
    print(k_means_data)


    def rotate(key):
        y_arr = k_means_data[key]
        x_arr = range(2, len(y_arr) + 2)
        last_point_y = float(y_arr[-1]) - float(y_arr[0])
        last_point_x = x_arr[-1] - 2
        angle_radians = -np.arctan2(last_point_y, last_point_x)
        #print(last_point_y)
        #print(last_point_x)
        #print(math.degrees(angle_radians))
        transformed_y_arr = []
        regular_y_arr = []
        for index in range(len(y_arr)):
            #y' = y cos f + x sin f
            transformed_y_arr.append((float(y_arr[index])-float(y_arr[0]))  * math.cos(angle_radians) + (x_arr[index] - 2) * math.sin(angle_radians))
            regular_y_arr.append(float(y_arr[index]))
        #print(transformed_y_arr)
        return regular_y_arr, transformed_y_arr

    intracohesion_max, transformed_intracohesion_max = rotate("intracohesion_max")
    intercohesion_max, transformed_intercohesion_max = rotate("intercohesion_max")
    #print(transformed_intracohesion_max)
    #print(k_means_data["intracohesion_max"])

    plt.figure(1)
    plt.title('Inter/Intra-Cohesion Values')
    ax1 = plt.subplot(211)
    ax1.set_title('Intracohesion v Num. Means')
    ax1.plot(range(2, len(intracohesion_max) + 2), intracohesion_max, color='green', marker='v', label='intragroup cohesion minimum')
    ax2 = plt.subplot(212)
    ax2.set_title('Transformed Intracohesion v Num. Means')
    ax2.plot(range(2, len(intracohesion_max) + 2), transformed_intracohesion_max, color='red', marker='v', label='intragroup cohesion minimum')

    plt.figure(2)
    plt.title('Inter-Cohesion Values')
    ax1 = plt.subplot(211)
    ax1.set_title('Intercohesion v Num. Means')
    ax1.plot(range(2, len(intercohesion_max) + 2), intercohesion_max, color='green', marker='v', label='intragroup cohesion minimum')
    ax2 = plt.subplot(212)
    ax2.set_title('Transformed Intercohesion v Num. Means')
    ax2.plot(range(2, len(intercohesion_max) + 2), transformed_intercohesion_max, color='red', marker='v', label='intragroup cohesion minimum')

    ideal_means = np.argmax(transformed_intracohesion_max)
    print("IDEAL AMOUNT OF MEANS: " + str(ideal_means+2))
    centroids = SequenceKMeans.test(sequences_input=sequences, names_input=names, initial_centroids=k_means_data["centroids_max"][ideal_means], debug=True)
    print(centroids)

    if plot:
        plt.show()
    return centroids


def load_sequences(filename=None, sequence_type="primary", head=None, tail=None, num_index_per_group=9, num_index_primary=2, num_index_secondary=7, select_top_sequence=False):
    if filename != None:
        print("Reading...")
        ## Read sequences from file
        sequences = DataReader.convert_file_to_array(filename)

        ## Index 2 corresponds to primary structure.
        primary_sequences = DataReader.strip_array_by_indices(sequences,num_index_per_group,[num_index_primary])

        ## Index 7 corresponds to secondary structure. The sequence data is formatted to one character representations.
        secondary_sequences = DataReader.strip_array_by_indices(sequences,num_index_per_group,[num_index_secondary])
        secondary_input_dictionary = {'S':'a','B':'b', 'D':'c', 'DH':'d', 'BH':'e', 'IS':'f', 'IA':'g','ISH':'h', 'IAH':'i', 'M':'j', 'MH':'k'}
        secondary_sequences = DataReader.replace_items_in_array(secondary_sequences, secondary_input_dictionary)

        ## Index 0 corresponds to the name of the sequence.
        generated_names = DataReader.strip_array_by_indices(sequences,num_index_per_group,[0])

        ## Removes the head/tail primer regions for primary sequences
        if sequence_type == "primary":
            generated_sequences = SequenceFunctions.strip_head_tail(sequences=primary_sequences, head=head, tail=tail)
        else:
            generated_sequences = SequenceFunctions.strip_head_tail(sequences=secondary_sequences, head=None, tail=None)

        ## This is to remove any sequences that doubly occur. (Had two variations of structures in the dataset)
        if select_top_sequence:
            indices = []
            for i in range(len(generated_names)):
                if generated_names[i][0][generated_names[i][0].index("-") + 1] == '1':
                    indices.append(i)
            generated_sequences = DataReader.strip_array_by_indices(generated_sequences,len(sequences),indices)
            generated_names = DataReader.strip_array_by_indices(generated_names,len(sequences),indices)


    else:

        patterns = ["GGGGG", "AGAGA", "AAAAA", "CCCCC", "CTCTC", "TTTTT"]
        weights = [float(1)/len(patterns) for p in patterns]
        length_of_sequence = 7
        number_of_sequences = 50

        generated_sequences, generated_names = SequenceFunctions.generate_random_sequences(
            patterns=patterns,
            weights=weights,
            length_of_sequence=length_of_sequence,
            number_of_sequences=number_of_sequences
        )

    return generated_names, generated_sequences

def k_means_run(filename, sequence_type, generated_sequences, lower_mean, upper_mean, num_algo_restarts, num_algo_iterations, plot=False):
    print(len(generated_sequences))
    centroids_overall = []
    cohesion_overall = []
    intercohesion_overall = []
    means = range(lower_mean,upper_mean)
    for j in means:
        cohesion_results = []
        intercohesion_results = []
        centroids_results = []
        for i in range(num_algo_restarts):
            cohesion_overall_list, intercohesion_overall_list, centroids = SequenceKMeans.main(generated_sequences, num_means=j, num_iterations=num_algo_iterations)
            cohesion_results.append(cohesion_overall_list)
            intercohesion_results.append(intercohesion_overall_list)
            centroids_results.append(centroids)
        index_min = np.argmin([item[-1] for item in intercohesion_results])
        index_max = np.argmax([item[-1] for item in cohesion_results])
        centroids_overall.append((centroids_results[index_min], centroids_results[index_max]))
        cohesion_overall.append((cohesion_results[index_min], cohesion_results[index_max]))
        intercohesion_overall.append((intercohesion_results[index_min], intercohesion_results[index_max]))
        with safe_open("Data/k_means_data/" + filename.split('/')[-1][:-4] + "/" + sequence_type + "/means" + str(j) + ".csv", "w") as f:
            #writer = csv.DictWriter(f, fieldnames=["cohesion_min", "cohesion_max", "intercohesion_min", "intercohesion_max", "centroids_min", "centroids_max"])
            writer = csv.writer(f)
            writer.writerows([["intracohesion_min"]])
            writer.writerows([cohesion_results[index_min]])
            writer.writerows([["intracohesion_max"]])
            writer.writerows([cohesion_results[index_max]])
            writer.writerows([["intercohesion_min"]])
            writer.writerows([intercohesion_results[index_min]])
            writer.writerows([["intercohesion_max"]])
            writer.writerows([intercohesion_results[index_max]])
            writer.writerows([["centroids_min"]])
            writer.writerows([centroids_results[index_min]])
            writer.writerows([["centroids_max"]])
            writer.writerows([centroids_results[index_max]])
        for i in range(num_algo_restarts):
            plt.figure(j)
            plt.title('Inter/Intra-Cohesion Values')
            ax1 = plt.subplot(111)
            ax1.set_title('Group Cohesion v Iterations (Num. Means: ' + str(j) + ')')
            ax1.plot(range(len(cohesion_results[i])),cohesion_results[i], color='black', marker='o')
            ax1.plot(range(len(intercohesion_results[i])),intercohesion_results[i], color='black', marker='o')

        plt.figure(j)
        plt.title('Inter/Intra-Cohesion Values')
        ax1 = plt.subplot(111)
        ax1.set_title('Group Cohesion v Iterations (Num. Means: ' + str(j) + ')')
        ax1.plot(range(len(cohesion_results[index_min])),cohesion_results[index_min], color='green', marker='o')
        ax1.plot(range(len(intercohesion_results[index_min])),intercohesion_results[index_min], color='red', marker='o')
        plt.xlabel('Number of Iterations')
        plt.ylabel('Cohesion Value')

    plt.figure(100)
    plt.title('Inter/Intra-Cohesion Values')
    ax1 = plt.subplot(111)
    ax1.set_title('Cohesion v Num. Means')
    min_distance_index = np.argmax([cohesion_overall[i][0][-1] - intercohesion_overall[i][0][-1] for i in range(len(cohesion_overall))])
    max_distance_index = np.argmax([cohesion_overall[i][1][-1] - intercohesion_overall[i][1][-1] for i in range(len(cohesion_overall))])
    ax1.axvline(means[min_distance_index],color='blue', label='k value for greatest distance between maximum cohesion values')
    ax1.axvline(means[min_distance_index],color='#95d0fc', linestyle='--', label='k value for greatest distance between maximum cohesion values')
    ax1.plot(means,[item[0][-1] for item in cohesion_overall], color='green', marker='v', label='intragroup cohesion minimum')
    ax1.plot(means,[item[0][-1] for item in intercohesion_overall], color='red', marker='v', label='intergroup cohesion minimum')
    ax1.plot(means,[item[1][-1] for item in cohesion_overall], color='#96f97b', marker='^', label='intragroup cohesion maximum')
    ax1.plot(means,[item[1][-1] for item in intercohesion_overall], color='#cf6275', marker='^', label='intergroup cohesion maximum')
    ax1.legend(loc='best', fontsize='xx-small')
    plt.xlabel('Number of Means')
    plt.ylabel('Cohesion Value')
    index = np.argmax([item[-1] for item in cohesion_overall])
    #print(centroids_overall[index])
    if plot:
        plt.show()

'''
Generates JSON file for D3 visualizer to use.
'''
def generate_d3_json_file(filename, primary_centroids, secondary_centroids):
    print('------')
    print(primary_centroids)
    print(secondary_centroids)
    print('_______')

    ## Gotta reverse map the key, because it's in the single character notation.
    def reversemapsecondarysequence(sequence):
        if sequence != "misc":
            secondary_input_dict = {'S':'a','B':'b', 'D':'c', 'DH':'d', 'BH':'e', 'IS':'f', 'IA':'g','ISH':'h', 'IAH':'i', 'M':'j', 'MH':'k'}
            secondary_input_dict_inverse = {v: k for k, v in secondary_input_dict.items()}
            new_sequence = ''
            for c in sequence:
                new_sequence += secondary_input_dict_inverse[c] + ' '
            return new_sequence
        else:
            return sequence
    with safe_open("Data/k_means_data/" + filename.split('/')[-1][:-4] + "_graphs.json", "w") as f:
        nodes = []
        links = []
        name = "group_primary_"
        count = 1
        with safe_open("Data/k_means_data/" + filename.split('/')[-1][:-4] + "_centroids_info.json", "w") as m:
            for key, values in primary_centroids.items():
                group_name = name +str(count)
                m.write(group_name + ': '+ key + '\n')
                for value in values:
                    print(value)
                    m.write("{0:.2f}".format(round(value['similarity'],2)) + ' || ' + value['name'] + ': ' + value['sequence'] + '\n')
                    nodes.append({"id": value['name'], "group": 1})
                    links.append({"source":value['name'],"target":group_name,"value":20 * value['similarity']})
                m.write('\n')
                nodes.append({"id": group_name, "group": 2})
                count += 1
            m.write('\n')
            name = "group_secondary_"
            count = 1
            for key, values in secondary_centroids.items():
                group_name = name +str(count)
                new_key = reversemapsecondarysequence(key)
                m.write(group_name + ': '+ new_key + '\n')
                for value in values:
                    m.write("{0:.2f}".format(round(value['similarity'],2)) + ' || ' + value['name'] + ': ' + reversemapsecondarysequence(value['sequence']) + '\n')
                    links.append({"source":value['name'],"target":group_name,"value":20 * value['similarity']})
                nodes.append({"id": group_name, "group": 3})
                count += 1
            json_dict = {"nodes":nodes, "links":links}
            json_obj = json.dump(json_dict, f)

'''
Running K-means on primary and secondary sequences of AKT2 oxidized 1 dataset.
'''
def Experiment1(run_k_means):
    arguments = {
    'run_k_means': True,
    'filename':'Data/Sequences/Processed/AKT2/Oxidized/akt2oxidized1.csv',
    'sequence_type': 'primary',
    'head': 'GGGACAGGGCTAGC',
    'tail': 'GAGGCAAAGCTTCCG',
    'num_index_per_group_in_file': 9,
    'num_index_primary_in_file': 2,
    'num_index_secondary_in_file': 7,
    'select_top_sequence': True,
    'lower_mean': 12,
    'upper_mean': 13,
    'num_algo_restarts': 15,
    'num_algo_iterations': 10,
    'plot': True
    }

    if run_k_means:
        k_means(**arguments)
    arguments['run_k_means'] = False
    primary_centroids = k_means(**arguments)

    arguments['sequence_type'] = 'secondary'
    if run_k_means:
        arguments['run_k_means']= True
        k_means(**arguments)
        arguments['run_k_means'] = False
    secondary_centroids = k_means(**arguments)

    generate_d3_json_file(filename=arguments['filename'], primary_centroids=primary_centroids, secondary_centroids=secondary_centroids)

    plt.plot()

'''
Running K-means on primary and secondary sequences of AKT2 oxidized 2 dataset.
'''
def Experiment2(run_k_means):
    arguments = {
    'run_k_means': True,
    'filename':'Data/Sequences/Processed/AKT2/Oxidized/akt2oxidized2.csv',
    'sequence_type': 'primary',
    'head': 'GGGACAGGGCTAGC',
    'tail': 'GAGGCAAAGCTTCCG',
    'num_index_per_group_in_file': 9,
    'num_index_primary_in_file': 2,
    'num_index_secondary_in_file': 7,
    'select_top_sequence': True,
    'lower_mean': 2,
    'upper_mean': 15,
    'num_algo_restarts': 15,
    'num_algo_iterations': 10,
    'plot': False
    }
    if run_k_means:
        k_means(**arguments)
    arguments['run_k_means'] = False
    primary_centroids = k_means(**arguments)

    arguments['sequence_type'] = 'secondary'
    if run_k_means:
        arguments['run_k_means']= True
        k_means(**arguments)
        arguments['run_k_means'] = False
    secondary_centroids = k_means(**arguments)

    generate_d3_json_file(filename=arguments['filename'], primary_centroids=primary_centroids, secondary_centroids=secondary_centroids)

    plt.plot()

'''
Running K-means on primary and secondary sequences of AKT2 oxidized combined dataset.
'''
def Experiment3(run_k_means):
    arguments = {
    'run_k_means': True,
    'filename':'Data/Sequences/Processed/AKT2/Oxidized/akt2oxidizedcombined.csv',
    'sequence_type': 'primary',
    'head': 'GGGACAGGGCTAGC',
    'tail': 'GAGGCAAAGCTTCCG',
    'num_index_per_group_in_file': 9,
    'num_index_primary_in_file': 2,
    'num_index_secondary_in_file': 7,
    'select_top_sequence': True,
    'lower_mean': 2,
    'upper_mean': 15,
    'num_algo_restarts': 15,
    'num_algo_iterations': 10,
    'plot': False
    }
    #print(run_k_means)
    if run_k_means:
        k_means(**arguments)
    arguments['run_k_means'] = False
    primary_centroids = k_means(**arguments)

    arguments['sequence_type'] = 'secondary'
    if run_k_means:
        arguments['run_k_means']= True
        k_means(**arguments)
        arguments['run_k_means'] = False
    secondary_centroids = k_means(**arguments)

    generate_d3_json_file(filename=arguments['filename'], primary_centroids=primary_centroids, secondary_centroids=secondary_centroids)

    plt.plot()

if __name__ == "__main__":
    opts = getopts(argv)
    experiments_dict = {
        'Experiment1': Experiment1,
        'Experiment2': Experiment2,
        'Experiment3': Experiment3
    }
    try:
        func = experiments_dict.get(opts['-e'], None)
        run = opts.get('-r', False)
        if func != None:
            func(run)
    except Exception as e:
        print(e)
        print('Example: python3 Experiments.py -e Experiment1')
