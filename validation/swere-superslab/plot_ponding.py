import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.font_manager as fm
from scipy.io import netcdf

# Define file directories and labels
fdir = ['out-async/','out-async/']
lgd = ['ParFlow', 'CATHY', 'HGS', 'Cast3M', 'SERGHEI', 'SERGHEI-asy']

# Define file names and read data
dataname = 'data/'
datasets = {
    'Ponding': ['PF-slab-ponding.csv', 'CATHY-slab-ponding.csv', 'HGS-slab-ponding.csv', 'cast3m-slab-ponding.csv'],
    'Outflow': ['PF-slab-outflow.csv', 'CATHY-slab-outflow.csv', 'HGS-slab-outflow.csv', 'cast3m-slab-outflow.csv']
}

data = {key: [pd.read_csv(dataname + fname, header=None) for fname in fnames] for key, fnames in datasets.items()}

# Define colors and markers
colors = [(19/255, 103/255, 158/255), (171/255, 58/255, 41/255), (208/255, 127/255, 44/255), (111/255, 109/255, 161/255)]
markers = ['o', '^', 's', '*']

# Read domain time series data
t, p, q = [], [], []
for ff in fdir:
    out = np.genfromtxt(ff + 'domainTimeSeries.out', skip_header=1)
    t.append(out[:, 0])
    p.append(out[:, 1])
    q.append(out[:, 4])

# Plotting
plt.figure(1, figsize=[9, 4])
plt.rc('font', size=10)
ss = 14

for idx, (key, ylabel, ydata) in enumerate(zip(['Ponding', 'Outflow'], 
                                               ['Ponding storage [$m^{3}$]', 'Flow Rate [$m^{3}$/h]'], 
                                               [p, q])):
    plt.subplot(1, 2, idx + 1)
    for i, (dataset, color, marker) in enumerate(zip(data[key], colors, markers)):
        plt.scatter(dataset.iloc[:, 0], dataset.iloc[:, 1], s=ss, marker=marker, facecolor='None', edgecolor=color)
    plt.plot(t[0] / 3600, ydata[0] if idx == 0 else ydata[0] * 3600, 'k-')
    plt.xlabel('Time [h]')
    plt.ylabel(ylabel)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    if idx == 0:
        plt.legend(lgd)
plt.savefig('ponding-flow.png', format='png', bbox_inches='tight', dpi=600)
plt.show()
