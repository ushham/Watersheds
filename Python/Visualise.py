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

num = 100

# plt.imshow(ele((num, num)), cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.show()
