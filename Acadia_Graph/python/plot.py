import scipy.io
import os 
import matplotlib.pyplot as plt
import visualization as viz
import numpy as np
import json
import sys



#print(mat)






def main():
    #setting working directory:
    dir_path = os.path.dirname(os.path.realpath(__file__))
    print(dir_path)
    os.chdir(dir_path)


    mat = scipy.io.loadmat(sys.argv[1])

    #"./config.json"

    json_file = open(sys.argv[3],"r",encoding="utf-8")
    config = json.load(json_file)
    json_file.close()

    viz.gen_plots(mat,sys.argv[2],config)


if __name__== "__main__":
    main()