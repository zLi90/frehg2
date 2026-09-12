import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import netcdf
import matplotlib.font_manager as fm

# Define file directories and labels
fdir = ['out-async/','out-async/']
lgd = ['ParFlow', 'CATHY', 'HGS', 'Cast3M', 'SERGHEI', 'SERGHEI-asy']
wcs = 0.1

# Define dataset names and corresponding labels
dataset_info = {
    'pft3x0': 'PF-t3-x0.csv', 'pft3x8': 'PF-t3-x8.csv', 'pft3x40': 'PF-t3-x40.csv', 'pft6x0': 'PF-t6-x0.csv',
    'chyt3x0': 'CATHY-t3-x0.csv', 'chyt3x8': 'CATHY-t3-x8.csv', 'chyt3x40': 'CATHY-t3-x40.csv', 'chyt6x0': 'CATHY-t6-x0.csv',
    'chyt6x8': 'CATHY-t6-x8.csv', 'chyt6x40': 'CATHY-t6-x40.csv', 'HGSt3x0': 'HGS-t3-x0.csv', 'HGSt3x8': 'HGS-t3-x8.csv',
    'HGSt3x40': 'HGS-t3-x40.csv', 'HGSt6x0': 'HGS-t6-x0.csv', 'HGSt6x8': 'HGS-t6-x8.csv', 'HGSt6x40': 'HGS-t6-x40.csv',
    'cast3mt3x0': 'cast3m-t3-x0.csv', 'cast3mt3x8': 'cast3m-t3-x8.csv', 'cast3mt3x40': 'cast3m-t3-x40.csv',
    'cast3mt6x0': 'cast3m-t6-x0.csv', 'cast3mt6x8': 'cast3m-t6-x8.csv', 'cast3mt6x40': 'cast3m-t6-x40.csv'
}

# Read datasets
data = {key: pd.read_csv(f'data/{fname}', header=None) for key, fname in dataset_info.items()}

# Define colors and markers
colors = [(19/255, 103/255, 158/255), (171/255, 58/255, 41/255), (208/255, 127/255, 44/255), (111/255, 109/255, 161/255)]
markers = ['o', '^', 's', '*']

# Read domain time series data
t, v, q = [], [], []
wc3x0, wc3x8, wc3x40, wc6x0, wc6x8, wc6x40 = [], [], [], [], [], []

for ff in fdir:
    out = np.genfromtxt(ff + 'domainTimeSeries.out', skip_header=1)
    t.append(out[:, 0])
    v.append(out[:, 1])
    q.append(out[:, 4])

    fname = ff + 'output_subsurface.nc'
    with netcdf.NetCDFFile(fname, 'r') as fid:
        wc = fid.variables['wc'][:]
    print(np.shape(wc))

    wc3x0.append(wc[6, :, 0, 0] / wcs)
    wc3x8.append(wc[6, :, 0, 31] / wcs)
    wc3x40.append(wc[6, :, 0, 159] / wcs)
    wc6x0.append(wc[12, :, 0, 0] / wcs)
    wc6x8.append(wc[12, :, 0, 31] / wcs)
    wc6x40.append(wc[12, :, 0, 159] / wcs)

z = np.linspace(-4.95, -0.05, 50)

# Plotting
plt.figure(1, figsize=[16, 12])
plot_info = [
    ('pft3x0', 'chyt3x0', 'HGSt3x0', 'cast3mt3x0', wc3x0, "X=0m", "T=3h"),
    ('pft3x8', 'chyt3x8', 'HGSt3x8', 'cast3mt3x8', wc3x8, "X=8m", None),
    ('pft3x40', 'chyt3x40', 'HGSt3x40', 'cast3mt3x40', wc3x40, "X=40m", None),
    ('pft6x0', 'chyt6x0', 'HGSt6x0', 'cast3mt6x0', wc6x0, "X=0m", "T=6h"),
    ('chyt6x8', 'HGSt6x8', 'cast3mt6x8', None, wc6x8, "X=8m", None),
    ('chyt6x40', 'HGSt6x40', 'cast3mt6x40', None, wc6x40, "X=40m", None)
]

for i, (pf_key, chy_key, hgs_key, cast_key, wc, title, text) in enumerate(plot_info):
    plt.subplot(2, 3, i+1)
    if pf_key in data:
        plt.scatter(data[pf_key].iloc[:, 0], data[pf_key].iloc[:, 1] * -1, s=30, marker='o', facecolor='None', edgecolor=colors[0])
    if chy_key in data:
        plt.scatter(data[chy_key].iloc[:, 0], data[chy_key].iloc[:, 1] * -1, s=30, marker='^', facecolor='None', edgecolor=colors[1])
    if hgs_key in data:
        plt.scatter(data[hgs_key].iloc[:, 0], data[hgs_key].iloc[:, 1] * -1, s=30, marker='s', facecolor='None', edgecolor=colors[2])
    if cast_key in data:
        plt.scatter(data[cast_key].iloc[:, 0], data[cast_key].iloc[:, 1] * -1, s=30, marker='*', facecolor='None', edgecolor=colors[3])
    plt.plot(np.flipud(wc[0]), z, 'k-')
    if len(wc) > 1:
        plt.plot(np.flipud(wc[1]), z, 'r--')
    plt.title(title)
    if text:
        plt.text(0.05, 0.95, text, transform=plt.gca().transAxes, fontsize=18, verticalalignment='top')
    if i == 0:
        plt.legend(lgd)

for i in range(6):
    plt.subplot(2, 3, i+1)
    plt.ylabel('Z [m]')
    plt.xlabel('Saturation [-]')
    plt.xticks([0.2, 0.4, 0.6, 0.8, 1.0])
    plt.xlim(0, 1.0)
    plt.ylim(-5, 0)

# plt.savefig('saturation.png', format='png', bbox_inches='tight', dpi=600)
plt.show()
