import numpy as np
import math
import rasterio as rs
import matplotlib.pyplot as plt

class PeakFinder:
    sig = 2
    res = (10, 10)
    def __init__(self, ele):
        self.ele = ele
    
    def gaussian_kernal(self, tup=res):
        def gaussian(x, y):
            pt1 = 1 / (2 * math.pi * self.sig)
            pt2 = math.exp(-(x**2 + y**2) / (2 * self.sig))
            return pt1 * pt2

        n, m = tup
        kernal = np.array([[gaussian(x, y) for y in np.linspace(-m, m, 2*m)] for x in np.linspace(-n, n, 2*n)])
        
        return kernal / kernal.sum()

    def conv(self, arr1, resolution=res):
        arr2 = self.gaussian_kernal(resolution)
        
        out_arr = np.empty_like(arr1)
        for i in range(arr1.shape[0]):
            rows = [min(max(0, x), arr1.shape[0]-1) for x in range((i - resolution[0]), (i + resolution[0]), 1)]
            for j in range(arr1.shape[1]):
                cols = [min(max(0, y), arr1.shape[1]-1) for y in range((j - resolution[1]),(j + resolution[1]), 1)]
                arr_temp = np.array([[arr1[k, l] for l in cols] for k in rows])
             
                out_arr[i, j] = np.sum(np.multiply(arr_temp, arr2))

        return out_arr

    def ele_data(self, file, input_type=False):
        if input_type:
            data = np.genfromtxt(file, delimiter=',')
        else:
            fl = rs.open(file)
            data = fl.read(1)
            data = data[7850:7950, 9650:9750]

        return data


# p = PeakFinder(1)
# arr = p.ele_data('50N030W_20101117_gmted_mea075.tif')

# num = 5
# res = (num, num)
# v = p.conv(arr, res)

# arr_out = arr - v
# print(np.average(arr_out))
# #arr_out = arr_out / np.average(arr_out)

# plt.imshow(arr_out>20, cmap='rainbow', vmin=0, vmax=1)
# plt.colorbar()
# plt.show()