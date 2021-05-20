import PeakConnect as pc
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import math

class GradDecent:
    lmda = 10
    del_h = 200

    def __init__(self, data):
        self.ele = data
        self.tree_x = np.zeros_like(self.ele)
        self.tree_y = np.zeros_like(self.ele)

        # #interpolate size
        # s = data.shape
        # X = np.linspace(0, s[1], s[1])
        # Y = np.linspace(0, s[0], s[0])

        # x, y = np.meshgrid(X, Y)
        # f = interpolate.interp2d(x, y, data, kind='linear')

        # Xnew = np.linspace(0, s[1], new_x)
        # Ynew = np.linspace(0, s[0], new_y)
        # self.ele = f(Xnew, Ynew)

    def gradient_decent(self, loc_x, loc_y):
        #use second order centre difference approximation for gradient:
        #[f(x+h)-f(x-h)]/2h
        
        #check if it is an edge point:
        if loc_x == 0 or loc_x == self.ele.shape[1]-1:
            return loc_x, loc_y

        elif loc_y == 0 or loc_y == self.ele.shape[0]-1:
            return loc_x, loc_y
        
        else:
            ini_h = self.ele[loc_y, loc_x]
            grad_x = (self.ele[loc_y, loc_x + 1] - self.ele[loc_y, loc_x - 1]) / (2 * self.del_h)
            grad_y = (self.ele[loc_y + 1, loc_x] - self.ele[loc_y - 1, loc_x]) / (2 * self.del_h)

            #find compass direction that ball will roll
           
            dirc = ((180 - (np.arctan2(grad_y , grad_x) * 180 / np.pi)) % 360)
           
            if (dirc < 67.5) or (dirc > 292.5):
                tst_x = 1
            elif (dirc > 112.5) and (dirc < 247.5):
                tst_x = -1
            else:
                tst_x = 0
            
            if (dirc > 22.5) and (dirc < 157.5):
                tst_y = -1
            elif (dirc > 202.5) and (dirc < 292.5):
                tst_y = 1
            else:
                tst_y = 0

            if (tst_x == 0) and (tst_y == 0):
                check = False
            else:
                check = True

            tst_x = loc_x + tst_x
            tst_y = loc_y + tst_y
            
            if False: #ini_h > self.data[tst_y, tst_x]:
                loc_x, loc_y = tst_x, tst_y
            else:
                tst_x, tst_y = self.gradient_decent_min(tst_x, tst_y)
                
                if (tst_x == loc_x) and (tst_y == loc_y):
                    check = False
                else:
                    loc_x, loc_y = tst_x, tst_y
        
        return loc_x, loc_y, check

    def gradient_decent_min(self, loc_x, loc_y, history, switch = True):
        #Simple check to ensure that there are no lower points in the 3x3 grid around original point.
        #If the values are the same altitude as starting point, recursively looks at those points to see if
        # they produce lower neighbour

        square = self.ele[loc_y-1:loc_y+2, loc_x-1:loc_x+2]

        min_val = np.min(square)
        i = np.where(square == min_val)
        prev_min = self.ele[loc_y, loc_x]
        tst_x, tst_y = loc_x, loc_y

        if len(i[0]) == 1:
            tst_x, tst_y = loc_x + i[1][0] - 1, loc_y + i[0][0] - 1
            history.append((tst_x, tst_y))
        else:
            for j in range(len(i[0])):
                #for each elevation of same height, find which elevation leads to lower height
                t_x, t_y = loc_x + i[1][j] - 1, loc_y + i[0][j] - 1

                chk_1 = (self.ele[t_y, t_x] <= prev_min)                            #Height is decreasing or steady
                chk_2 = ((t_x != loc_x) or (t_y != loc_y))                          #Not the same point as previous
                chk_3 = (self.tree_x[t_y, t_x] == 0 and self.tree_y[t_y, t_x] == 0) #Has not been included in previous iteration
                chk_4 = not((t_x, t_y) in history)                                  #Prevents circular searches
                chk_5 = (self.ele[t_y, t_x] >= np.min(self.ele)) and switch               
                
                switch = False if (self.ele[t_y, t_x] == np.min(self.ele)) else True

                if chk_1 and chk_2 and chk_3 and chk_4 and chk_5:
                    history.append((t_x, t_y))
                    return self.gradient_decent_min(t_x, t_y, history, switch)
                    
        
        return history

    def find_equil(self, loc_x, loc_y):
        tst_x, tst_y = 0, 0
        min_heigh = np.min(self.ele)
        check = True
        n = 0

        indi_run = [[loc_x, loc_y, self.ele[loc_y, loc_x]]]
        while check:
            points = self.gradient_decent_min(loc_x, loc_y, [])
           
            if points == []:
                check = False
            #ensure we are not repeating points
            for p in points:
                
                tst_x, tst_y = p
                if (self.tree_x[loc_y, loc_x] == 0) and (self.tree_y[loc_y, loc_x] == 0):
                    height = self.ele[tst_y, tst_x]
                    self.tree_x[loc_y, loc_x], self.tree_y[loc_y, loc_x] = tst_x, tst_y
                    indi_run.append([tst_x, tst_y, self.ele[tst_y, tst_x]])
                    loc_x, loc_y = tst_x, tst_y
    
                    if height == min_heigh:
                        check = False
                else:
                    check = False
            n += 1
        return indi_run

    def make_tree(self):
        #function to run through all points and make tree diagram
        for i in range(5, self.ele.shape[0]-5):
            for j in range(5, self.ele.shape[1]-5):
                print(self.find_equil(i, j))


p = pc.PeakFinder(1)
arr = p.ele_data('50N030W_20101117_gmted_mea075.tif')

num = 5
res = (num, num)
v = p.conv(arr, res)

arr_out = arr - v

d = GradDecent(arr)
d.make_tree()
# run = d.find_equil(208, 61)

# dat = np.array(run)

# fig, axs = plt.subplots(1, 1, figsize=(12, 8))
# axs.imshow(arr, cmap='rainbow')

# axs.plot(dat[:, 0], dat[:, 1], color='black')

# plt.show()

# plt.plot(dat[:, 2])
# plt.show()