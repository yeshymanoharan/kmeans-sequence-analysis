'''
Class for Excel file manipulation
'''
class DataReader():
    '''
    Converts file to array.
    Takes in parameters:
    filename - type(string)
    separator - type(string)
    '''
    def convert_file_to_array(filename, separator=None):
        csv_file = open(filename)
        line = csv_file.readline()
        csv_array = []
        while line is not '':
            line_array = line.split(separator)
            csv_array.append(line_array)
            line = csv_file.readline()
        
        return csv_array

    '''
    Converts array to file.
    Takes in parameters:
    filename, array, separator
    '''
    def convert_array_to_file(filename, array, separator=','):
        string = ""
        for row in array:
            for value in row:
                string += value + separator
            string += '\n'
        
        csv_file = open(filename, "w+")
        csv_file.write(string)
        csv_file.close()

    def strip_array_by_indices(array, numrowsforpoint, indices=[0]):
        new_array = []
        counter = 0
        for line in array:
            if (counter % numrowsforpoint) in indices:
                new_array.append(line)
            counter += 1
        return new_array
    
    def combine_arrays(array1, array2, name1="Array1", name2="Array2"):
        combined_array = [[name1,name2]]
        for i in range(0, min(len(array1), len(array2))):
            line1 = array1[i]
            line2 = array2[i]
            for j in range(0, min(len(line1), len(line2))):
                combined_array.append([line1[j], line2[j]])
        return combined_array

    def replace_items_in_array(array, replacements):
        new_array = []
        for i in range(0, len(array)):
            row = array[i]
            new_row = []
            for j in range(0, len(row)):
                if row[j] not in replacements:
                    new_row.append(row[j])
                else:
                    new_row.append(replacements.get(row[j]))
            new_array.append(new_row)
        return new_array
