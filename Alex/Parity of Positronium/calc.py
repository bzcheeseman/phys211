import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


datafile = r'data/Alex_Positronium.csv'
data = np.array(pd.read_csv(datafile, usecols=(8,)))
data2 = np.array(pd.read_csv(datafile, usecols=(11,)))
parallel = [data[0],data[2],data[4],data[6],data[8],data[10],data[13],data[15]]
perpendicular = [data[1],data[3],data[5],data[7],data[9],data[11],data[12],data[14]]
parallel_err = [data2[0],data2[2],data2[4],data2[6],data2[8],data2[10],data2[13],data2[15]]
perpendicular_err = [data2[1],data2[3],data2[5],data2[7],data2[9],data2[11],data2[12],data2[14]]

para = np.sum(parallel)/16
para_err = np.sum(parallel_err)/16
perp = np.sum(perpendicular)/16
perp_err = np.sum(perpendicular_err)/16


ratio = perp/para
ratio_err = ratio * (para_err/para + perp_err/perp)

print max(parallel)
print min(perpendicular)