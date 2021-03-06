import csv
import os 
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series
import itertools
import scipy.stats as st

def loadEdgeList(filename):
    data = pd.read_csv(filename)
    return data
    
def calculateClassSizes(edge_list):
    # Get number of edges for each class
    class_data = edge_list[edge_list['Type'] == 'C']
    counts = class_data.groupby('Info', as_index=False).count()
    
    # Solve quadratic equation n^2 - n - 2*edges = 0
    counts['Size'] = (np.sqrt(counts['Type']*8+1).astype(int)+1)//2
    return counts[['Info','Size']].copy()    

def edgesAsMatrices(edge_list, n=None):
    # Get count of number of repeats of personA->person
    counts = edge_list.groupby(['PersonA','PersonB','Type'], as_index=False).count()

    highest_id = 0
    if n is None:
        highest_id = max(max(edge_list['PersonA']), max(edge_list['PersonB']))
    else:
        highest_id = n

    c_matrix = np.zeros((highest_id+1,highest_id+1))
    m_matrix = c_matrix.copy()
    r_matrix = c_matrix.copy()
    s_matrix = c_matrix.copy()

    c_counts = counts[counts['Type'] == 'C']
    c_matrix[c_counts['PersonA'].to_numpy(), c_counts['PersonB'].to_numpy()] = c_counts['Info'].to_numpy()
    c_matrix = c_matrix + np.transpose(c_matrix)

    m_counts = counts[counts['Type'] == 'M']
    m_matrix[m_counts['PersonA'].to_numpy(), m_counts['PersonB'].to_numpy()] = m_counts['Info'].to_numpy()
    m_matrix = m_matrix + np.transpose(m_matrix)

    r_counts = counts[counts['Type'] == 'R']
    r_matrix[r_counts['PersonA'].to_numpy(), r_counts['PersonB'].to_numpy()] = r_counts['Info'].to_numpy()
    r_matrix = r_matrix + np.transpose(r_matrix)

    s_counts = counts[counts['Type'] == 'S']
    s_matrix[s_counts['PersonA'].to_numpy(), s_counts['PersonB'].to_numpy()] = s_counts['Info'].to_numpy()
    s_matrix = s_matrix + np.transpose(s_matrix)
    
    return c_matrix, m_matrix, r_matrix, s_matrix


def student_roommates_dist(freq=None):
    # Default dist is based on:
    # https://books.google.ca/books?id=k8IyAAAAIAAJ&pg=PA27&lpg=PA27&dq=distribution+of+number+of+roommates+student&source=bl&ots=dhntlrZkqY&sig=ACfU3U2TOHZvAnGaUpX544Wesp9LM9roRA&hl=en&sa=X&ved=2ahUKEwj64diJo6_qAhUFknIEHfSvCiAQ6AEwAHoECAsQAQ#v=onepage&q=distribution%20of%20number%20of%20roommates%20student&f=false
    if freq is None:
        freq = np.array([209,507,206,152])
    
    # this is the dist of number of students with apartments of a given size
    base_dist = (freq/np.sum(freq)) 
    
    # Readjust distribution such that the frequency is over the people in them, not
    # over the apartment sizes.
    vals = np.arange(freq.size)+1
    dist_adj = base_dist/vals
    dist_adj = dist_adj/np.sum(dist_adj)
    
    return st.rv_discrete(values=(vals, dist_adj))


def random_roommate_assignment(res_matrix,dist_fn):
    # Build off-campus matrix, all zeros
    sz = res_matrix.shape[0]
    offcampus = np.zeros((sz,sz))
    
    # Determine who does not live on campus
    res_mates = np.sum(res_matrix,axis=1)
    off_campus_people = np.where(res_mates==0)[0]
    
    # Make bins (apartments) of sizes distributed by dist_fn
    n = off_campus_people.size
    range = np.arange(n)
    group_ids = np.repeat(range, dist_fn.rvs(size=n))
    group_ids = group_ids[range]
    
    # Randomly order the off-campus people
    p2 = off_campus_people.copy()
    np.random.shuffle(p2)
    
    # For each apartment, make each pair of people in it adjacent
    groups = np.unique(group_ids)
    for id in groups:
        members = p2[group_ids==id]
        for pair in itertools.combinations(members,2):
            offcampus[pair[0]][pair[1]] = 1
    
    # Make bi-directional
    offcampus = offcampus + np.transpose(offcampus)
    
    return offcampus


#Here's a change
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


#def generate_random_social_graph(CM_graph,RS_graph,avg_friends,dispersion,avg_contacts_per_day):
def generate_random_social_graph(CM_graph,avg_friends,dispersion):

    #students_CM = pd.concat([pd.Series(CM_graph.PersonA.unique()),pd.Series(CM_graph.PersonB.unique())],ignore_index = True)
    #students_CM = students_CM.unique()
    #totalNodes = len(students_CM)

    totalNodes = CM_graph.shape[0]

    x = nx.watts_strogatz_graph(totalNodes,avg_friends,dispersion)
    # nx.draw_networkx(x, node_size=10)

    matrix = nx.to_numpy_matrix(x)

    #for elem in np.nditer(matrix):
    #    elem = int(elem)
    #    if elem == 1:
    #        elem = np.random.poisson(avg_contacts_per_day)
            

    matrix = np.array(matrix)
    return matrix


def partition_subgraphs(graph_data,group_type,info_array,groups):
    
    result = pd.DataFrame(columns=['PersonA','PersonB','Type','Info'])
    
    for info in info_array:
        
        data = graph_data.where((graph_data['Type']==group_type)&(graph_data['Info']== info)).dropna(axis=0, how='all')

        data.PersonA = data.PersonA.astype(int)
        data.PersonB = data.PersonB.astype(int)
        data.Info = data.Info.astype(int)

        #print(data.dtypes)

        group_members = pd.concat([pd.Series(data.PersonA.unique()),pd.Series(data.PersonB.unique())],ignore_index = True)

        #print(group_members)
        
        group_members = group_members.unique()

        #print(group_members)
        

        #print(np.array_split(group_members, groups))
        for part in np.array_split(group_members, groups):
            if len(part) > 1:
                lst = list(itertools.combinations(part,2))
                #print(lst)
                for x in range(len(lst)):
                    # print(lst[x][0],lst[x][1])
                    new_row = {'PersonA':lst[x][0],'PersonB':lst[x][1],'Type':group_type,'Info':info}
                    result = result.append(new_row,ignore_index=True)
            '''
            else:
                #print(part)
                #print(part[0])
                #safe case:
                new_row = {'PersonA':part[0],'PersonB':part[0],'Type':group_type,'Info':info}
                result = result.append(new_row,ignore_index=True)
            '''

    for info in info_array:
        indexNames = graph_data[ (graph_data['Type']==group_type)&(graph_data['Info']== info) ].index
        graph_data.drop(indexNames , inplace=True)
    
    
    result = pd.concat([result, graph_data], ignore_index=True)

    return result
