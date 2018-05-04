## Takes two .csv files and combines them.
primaryMatrix = open("primaryMatrix.csv")
secondaryMatrix = open("secondaryMatrix.csv")

primaryLine = primaryMatrix.readline()
secondaryLine = secondaryMatrix.readline()

points = 'primary_structure_similarity,secondary_structure_similarity\n'
while primaryLine != '' and secondaryLine != '':
    primaryArr = primaryLine.split()
    secondaryArr = secondaryLine.split()
    for i in range(len(primaryArr)):
        points += str(primaryArr[i]) + '\n'

    primaryLine = primaryMatrix.readline()
    secondaryLine = secondaryMatrix.readline()

pointsFile = open("points.csv", "w+")
print(points)
pointsFile.write(points)
pointsFile.close()