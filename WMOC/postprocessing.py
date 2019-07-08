# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 11:16:12 2019

@author: lx2347
"""
from __future__ import print_function
import matplotlib.pyplot as plt
def plot_head_history(pipe,H,wn,tt):
    pipeid = wn.links[pipe].id-1
    plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k') 
    plt.plot(tt,H[pipeid][0,:], 'b-',label='Start Node') 
    plt.plot(tt,H[pipeid][-1,:], 'r-',label='End Node') 
    plt.xlim([tt[0],tt[-1]])
    plt.title('Pressure Head of Pipe %s '%pipe)
    plt.xlabel("Time")
    plt.ylabel("Pressure Head (m)") 
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()
    
def plot_velocity_history(pipe,V,wn,tt):
    pipeid = wn.links[pipe].id-1
    plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(tt,V[pipeid][0,:], 'b-',label='Start Node') 
    plt.plot(tt,V[pipeid][-1,:], 'r-',label='End Node') 
    plt.xlim([tt[0],tt[-1]])
    plt.title('Velocity Head of Pipe %s ' %pipe)
    plt.xlabel("Time")
    plt.ylabel("Velocity (m/s)") 
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()




#%%
#read resultds from Hammer 
#root ='E:\\Box Sync\\utexas\\03_research_Lina\\working files\\InverseAnalysis\\HammerBC'
#filename='1-300.txt'
#df1 = pd.read_csv(os.path.join(root, filename), header=0, sep="\t",
#                 dtype={"Time (sec)": np.float64,
#                        "Base - Hydraulic Grade (m)": np.float64, 
#                        "Base - Flow (m3/s)":np.float64})
#df1= pd.DataFrame(df1)
#h300 = df1['Base - Hydraulic Grade (m)'].to_numpy()
#q300=df1['Base - Flow (m3/s)'].to_numpy()
#t = df1['Time (sec)'].to_numpy()
#
#
#filename='1-600.txt'
#df2 = pd.read_csv(os.path.join(root, filename), header=0, sep="\t",
#                 dtype={"Time (sec)": np.float64,
#                        "Base - Hydraulic Grade (m)": np.float64, 
#                        "Base - Flow (m3/s)":np.float64,
#                        "Base - Air/Vapor Volume (L)": np.float64})
#df2= pd.DataFrame(df2)
#h600 = df2['Base - Hydraulic Grade (m)'].to_numpy()
#q600 = df2['Base - Flow (m3/s)'].to_numpy()
#
#filename='1-1200.txt'
#df3 = pd.read_csv(os.path.join(root, filename), header=0, sep="\t",
#                 dtype={"Time (sec)": np.float64,
#                        "Base - Hydraulic Grade (m)": np.float64, 
#                        "Base - Flow (m3/s)":np.float64,
#                        "Base - Air/Vapor Volume (L)": np.float64})
#df3= pd.DataFrame(df3)
#h1200 = df3['Base - Hydraulic Grade (m)'].to_numpy()
#q1200 = df3['Base - Flow (m3/s)'].to_numpy()

