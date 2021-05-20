import numpy as np
import math
import matplotlib.pyplot as plt

# Function to control elevation
def ele(res):
    arr = np.empty(res)
    for i in range(res[0]):
        for j in range(res[1]):
            arr[i, j] = np.sin(j / res[1] * math.pi) * (i / res[0]) 
    return arr

#Array plot
def height_map(ele, x=[], y=[], trees = False):
    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    axs.imshow(ele, cmap='rainbow')

    #include tree diagram
    print("making tree diagram")
    if (x != []) and (y != []) and trees:
        for i in range(1, x.shape[1]-1):
            for j in range(1, y.shape[0]-1):
                if x[j, i] != 0 or y[j, i] != 0:
                    #print([i, j], [x[j, i], y[j, i]])
                    axs.plot([i, x[j, i]], [j, y[j, i]], color='black')
                
    plt.show()
    return 0
# #Creates 
# def tree_graph(x, y):


# plt.imshow(ele((num, num)), cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()
