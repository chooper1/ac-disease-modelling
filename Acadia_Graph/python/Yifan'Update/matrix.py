import csv
import os 
import numpy as np

#setting working directory:
dir_path = os.path.dirname(os.path.realpath(__file__))
print(dir_path)
os.chdir(dir_path)

#matrix for class and class tramission inside the same building:
def matrixCM(filename,C,M,classlimit):
    #find max
    max = 0
    #file is 'contact_graph.csv'
    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            if int(rows['PersonA']) > max:
                max = int(rows['PersonA'])
    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            if int(rows['PersonB']) > max:
                max = int(rows['PersonB'])

    #print("The max ID of a student:" , max)

    
    #create a table to count class size, the first two rows are left to zeroes intentionally for computing purpose
    table = np.zeros((2,(max+2)),dtype=np.int) 

    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            if rows['Type'] == "C":
                entry = int(rows['Info'])
                NotFound = True
                for i in range(int(table.size/(max+2))):
                    if table[i][0] == entry:
                        num1 = int(rows['PersonA'])
                        num2 = int(rows['PersonB'])
                        table[i][num1] += 1
                        table[i][num2] += 1
                        NotFound = False
                if(NotFound):
                    newRow = np.zeros(max+2)
                    num1 = int(rows['PersonA'])
                    num2 = int(rows['PersonB'])
                    newRow[0] = entry
                    newRow[num1] += 1
                    newRow[num2] += 1
                    table = np.vstack((table,newRow))
                      


    #calculate class size:
    for i in range(int(table.size/(max+2))):
        for j in range(1,table[i].size-1):
            if table[i][j] != 0:
                table[i][table[i].size-1] += 1

    #print(table)


    #create matrix
    matrix = np.zeros((max,max),dtype=np.float)

    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            #print(rows['Type'] == "C")
            if rows['Type'] == "C":
                # filter off class that is larger than class limit:

                #search class:
                 for i in range(0,int(table.size/(max+2))):
                     if table[i][0] == int(rows['Info']):
                         if table[i][table[i].size-1] < classlimit:
                            matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1] += C
                            matrix[int(rows['PersonB'])-1][int(rows['PersonA'])-1] += C
                            #To ensure
                            #print(matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1])
            if rows['Type'] == "M":
                matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1] += M
                matrix[int(rows['PersonB'])-1][int(rows['PersonA'])-1] += M
                #To ensure
                #print(matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1])
            
 
    #print matrix
    #print(matrix)
    #ensure not empty
    #print(np.sum(matrix))
    return matrix


#matrix for resident and social:
def matrixRS(filename,R,S):
    #find max
    max = 0
    #file is 'contact_graph.csv'
    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            if int(rows['PersonA']) > max:
                max = int(rows['PersonA'])
    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            if int(rows['PersonB']) > max:
                max = int(rows['PersonB'])

    print("The max ID of a student:" , max)

    #create matrix
    matrix = np.zeros((max,max),dtype=np.float)

    with open(filename, 'r') as file:
        csv_file = csv.DictReader(file)
        for rows in csv_file:
            #print(rows['Type'] == "C")
            if rows['Type'] == "R":
                matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1] += R
                matrix[int(rows['PersonB'])-1][int(rows['PersonA'])-1] += R
                #To ensure
                #print(matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1])
            
            if rows['Type'] == "S":
                matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1] += S
                matrix[int(rows['PersonB'])-1][int(rows['PersonA'])-1] += S
                #To ensure
                #print(matrix[int(rows['PersonA'])-1][int(rows['PersonB'])-1])
            
 
    #print matrix
    #print(matrix)
    #ensure not empty
    #print(np.sum(matrix))
    return matrix





#matrixCM("contact_graph.csv",2,3,50)

#matrixRS("contact_graph.csv",4,5)




