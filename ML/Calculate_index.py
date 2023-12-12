import numpy as np
import pandas as pd
import math as mt
import os

class Indices:
    # get values from the cvs file and divide them into arrays of subbands
    
    @classmethod
    def arr_freq(self, dB, band):
        low = []
        mid = []
        high = []
        counter = 0
        for i in band:
            if i <= 10*mt.pow(10, 3):
                low.append(dB[counter])
            if i >= 5*mt.pow(10, 3) and i <= 500*mt.pow(10, 3):
                mid.append(dB[counter])
            if i >= 400*mt.pow(10, 3):
                high.append(dB[counter])
            counter += 1
        return low, mid, high

    @classmethod
    def asle(self, dB1, dB2, band):
        low_1, mid_1, high_1 = self.arr_freq(dB1, band)
        low_2, mid_2, high_2 = self.arr_freq(dB2, band)
        
        def intermediate(band_1, band_2):
            index = 0
            for i in range(len(band_1)):
                index += abs(20*mt.log10(band_1[i]) + 20*mt.log10(band_2[i]))
            return index/max(len(band_1), len(band_2))
        
        return [intermediate(low_1, low_2), intermediate(mid_1, mid_2), intermediate(high_1, high_2)]
    
    @classmethod
    def cc(self, dB1, dB2, band):
        low_1, mid_1, high_1 = self.arr_freq(dB1, band)
        low_2, mid_2, high_2 = self.arr_freq(dB2, band)

        def intermediate(band_1, band_2):
            xy = 0
            x2 = 0
            y2 = 0

            for i in range(len(band_1)):
                xy += band_1[i]*band_2[i]
                x2 += band_1[i]*band_1[i]
                y2 += band_2[i]*band_2[i]
            
            return ((xy)/(mt.sqrt(x2*y2)))
        
        return [intermediate(low_1, low_2), intermediate(mid_1, mid_2), intermediate(high_1, high_2)]
    
    @classmethod
    def dabs(self, dB1, dB2, band):
        low_1, mid_1, high_1 = self.arr_freq(dB1, band)
        low_2, mid_2, high_2 = self.arr_freq(dB2, band)

        def intermediate(band_1, band_2):
            sub = 0

            for i in range(len(band_1)):
                sub += abs(band_2[i] - band_1[i])
            
            return sub/max(len(band_1), len(band_2))
        
        return [intermediate(low_1, low_2), intermediate(mid_1, mid_2), intermediate(high_1, high_2)]
    
    @classmethod
    def rmse(self, dB1, dB2, band):
        low_1, mid_1, high_1 = self.arr_freq(dB1, band)
        low_2, mid_2, high_2 = self.arr_freq(dB2, band)

        def intermediate(band_1, band_2):
            sum = 0

            for i in range(len(band_1)):
                sum = (band_1[i] - band_2[i])*(band_1[i] - band_2[i])
            
            return mt.sqrt(sum)
        
        return [intermediate(low_1, low_2), intermediate(mid_1, mid_2), intermediate(high_1, high_2)]

    
def calc():
    directory = './FRs'
    folders = [folder for folder in os.listdir(directory) if os.path.isdir(os.path.join(directory, folder))]
    for dirs in folders:
        asle = np.empty((0,3), int)
        cc = np.empty((0,3), int)
        dabs = np.empty((0,3), int)
        rmse = np.empty((0,3), int)
        di = []
        disk = []


        path = os.path.join(directory, dirs)
        org = np.array([])
        

        if 'FR' not in dirs:
            print("No FR folders are found")
            continue
        
        new_dir = "indices_" + dirs
        os.makedirs(new_dir, exist_ok=True)
        
        for file in os.listdir(path):
            if file.endswith('.txt'):
                if "healthy" in file:
                    data = np.loadtxt(path + '/' + file)
                    org = fra = data[:, 0]
        
        for file in os.listdir(path):
            print(file)
            if "healthy" not in file:
                dis = file.split('_')[4]
                dis = float(dis.replace('.net.txt', ''))
                disK = int(file.split('_')[2])
                disk.append(disK)
                di.append(dis)
                data = np.loadtxt(path + '/' + file)
                fra = data[:, 0]
                freq = data[:, 1]
                asle_temp = Indices.asle(fra, org, freq)
                asle = np.vstack((asle, asle_temp))
                cc_temp = Indices.cc(fra, org, freq)
                cc = np.vstack((cc, cc_temp))
                dabs_temp = Indices.dabs(fra, org, freq)
                dabs = np.vstack((dabs, dabs_temp))
                rmse_temp = Indices.rmse(fra, org, freq)
                rmse = np.vstack((rmse, rmse_temp))

        ind_values = [asle, cc, dabs, rmse]
        for i in range(len(ind_values)):
            ind_values[i] = np.insert(ind_values[i], 3, di, axis=1)
            ind_values[i] = np.insert(ind_values[i], 4, disk, axis = 1)

        ind_names = ['asle', 'cc', 'dabs', 'rmse']
        for i in range(len(ind_names)):
            filename = "index_{}_{}.txt".format(ind_names[i], file)
            file_path = os.path.join(new_dir, filename)
        
            with open(file_path, "w") as nt:
                np.savetxt(nt, ind_values[i])

calc()