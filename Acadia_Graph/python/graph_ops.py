import csv
import os 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series

#matrix for class and class tramission inside the same building:
def matrixCM(filename,C,M,classlimit):
    #find max
    max = 0
    
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

    #print("The max ID of a student:" , max)

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
                
    return matrix


def generate_random_social_graph(CM_graph,RS_graph,avg_friends,dispersion,avg_contacts_per_day):

    students_CM = pd.concat([pd.Series(CM_graph.PersonA.unique()),pd.Series(CM_graph.PersonB.unique())],ignore_index = True)
    students_CM = students_CM.unique()

           
    totalNodes = len(students_CM)

    x = nx.watts_strogatz_graph(totalNodes,avg_friends,dispersion)
    nx.draw_networkx(x, node_size=10)

    matrix = nx.to_numpy_matrix(x)

    for elem in np.nditer(matrix):
        elem = int(elem)
        #print(type(elem))
        if elem == 1:
            elem = np.random.poisson(avg_contacts_per_day)
            print("hello")
            print(elem)

    print(matrix.sum(dtype='float'))

    return matrix


def partition_subgraphs(graph_data,group_type,info_array,groups):
    
    for info in info_array:
        
        data = graph_data.where((graph_data['Type']==group_type)&(graph_data['Info']== info)).dropna(axis=0, how='all')

        data.PersonA = data.PersonA.astype(int)
        data.PersonB = data.PersonB.astype(int)
        data.Info = data.Info.astype(int)

        #print(data.dtypes)

        group_members = pd.concat([pd.Series(data.PersonA.unique()),pd.Series(data.PersonB.unique())],ignore_index = True)

        #print(group_members)
        
        group_members = group_members.unique()

        print(group_members)
        
        

        #creating new rows:

        df = pd.DataFrame(columns=['PersonA','PersonB','Type','Info'])


        print(np.array_split(group_members, groups))
        for part in np.array_split(group_members, groups):
            if len(part) > 1:
                lst = list(itertools.combinations(part,2))
                #print(lst)
                for x in range(len(lst)):
                    # print(lst[x][0],lst[x][1])
                    new_row = {'PersonA':lst[x][0],'PersonB':lst[x][1],'Type':group_type,'Info':info}
                    df = df.append(new_row,ignore_index=True)
            else:
                print(part)
                print(part[0])
                #safe case:
                new_row = {'PersonA':part[0],'PersonB':part[0],'Type':group_type,'Info':info}
                df = df.append(new_row,ignore_index=True)


        print(df)
