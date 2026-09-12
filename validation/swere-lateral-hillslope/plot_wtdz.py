"""
	Validate Serghei against Hydrus results reported in Beegum VZJ(2018)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf

fdir = ['out-pc-cuda/','out-mp-cuda/','out-pc-cuda-vardz/']
lgd = ['SERGHEI(PC, $\Delta z=$0.25m)','SERGHEI(MP, $\Delta z=$0.25m)','SERGHEI(PC, Variable $\Delta z$)','Hydrus']
tt = 60
L = 4000
H = 15
dz = 0.25
dim = [40, 1, 60]
N = int(dim[0]*dim[1]*dim[2])

r = [1.0, 1.0, 1.05]
dz_base = [0.25, 0.25, 0.1]
dz_vec = []
	

"""
	--------------------------------------------------------------------
						  Read HYDRUS results
	--------------------------------------------------------------------
"""
hydrus1 = np.genfromtxt('HYDRUS2D.csv',delimiter=',',skip_header=1)
hydrus2 = np.genfromtxt('hydrus-time.csv',delimiter=',',skip_header=1)

"""
	--------------------------------------------------------------------
						 Read Serghei results
	--------------------------------------------------------------------
"""
h_out = []
wc_out = []
dims = []
for ff in range(len(fdir)):
	fname = fdir[ff] + 'output_subsurface.nc'
	fid = netcdf.NetCDFFile(fname,'r')
	h_out.append(fid.variables['hd'])
	fid.close()
	dims.append(np.shape(h_out[ff]))
	print(' Shape of the output fields: ', np.shape(h_out[ff]))
	# get dz 
	tmp = []
	for kk in range(dims[ff][1]):
		tmp.append(dz_base[ff] * r[ff]**kk)
	dz_vec.append(np.array(tmp))

wt_out = []
x = np.linspace(0, L, dims[0][3])
ioutput = 0

print(dims)
for ff in range(len(fdir)):
	serghei = h_out[ff][-1,:,0,:]
	wt = np.zeros((dims[0][3]))
	for ii in range(dims[0][3]):
		for kk in range(dims[ff][1]-1):
			h1 = serghei[kk,ii]
			h2 = serghei[kk+1,ii]
			z1 = H - np.nansum(dz_vec[ff][:kk]) - 0.5*dz_vec[ff][kk]
			z2 = H - np.nansum(dz_vec[ff][:kk]) - 1.5*dz_vec[ff][kk]
			#z1 = H - (kk+0.5)*dz
			#z2 = H - (kk+1.5)*dz
			if h1 < 0.0 and h2 >= 0:
				wt[ii] = z1 - h1*(z1-z2)/(h1-h2)
				break
	wt_out.append(wt)

hcentre=[]
x_time=[]

for ff in range(len(fdir)):
	hcentre.append([])
	for i in range(tt+1):
		aa=h_out[ff][i,:,0,20]
		for kk in range(dims[ff][1]-1):
			h1 = aa[kk]
			h2 = aa[kk+1]
			z1 = H - np.nansum(dz_vec[ff][:kk]) - 0.5*dz_vec[ff][kk]
			z2 = H - np.nansum(dz_vec[ff][:kk]) - 1.5*dz_vec[ff][kk]
			#z1 = H - (kk+0.5)*dz
			#z2 = H - (kk+1.5)*dz
			if h1 < 0.0 and h2 >= 0:
				bb = z1 - h1*(z1-z2)/(h1-h2)
				hcentre[ff].append(bb) 
				break
for i in range(tt):
	a=i*2592000/(3600*24)
	x_time.append(a)
x_time.append(1825)


"""
	--------------------------------------------------------------------
								Make plot
	--------------------------------------------------------------------
"""
color1=[[143/255,38/255,126/255],[141/255,159/255,209/255],[143/255,38/255,126/255]]
ls = ['-','-','--']

fs = 9
cm2inch = 1.0/2.54

idx = 1

fig = plt.figure(1, figsize=[14*cm2inch, 12*cm2inch])
ifig = 1
icurve = 0

plt.subplot(2,1, 1)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0-0.01, pos1.y0+0.05, pos1.width*1.1, pos1.height*1.05]
ax.set_position(pos2)
xticks = range(0, 5 * 365+1, 365)
xtick_labels = [str(i) for i in range(0, 5 * 365+1, 365)]
for ff in range(len(fdir)):
	plt.plot(x_time[1:], hcentre[ff][1:], color=color1[ff], linestyle=ls[ff])  
plt.scatter(hydrus2[:-1,0], hydrus2[:-1,1],  marker='^', facecolor='None', edgecolor='k')  
plt.legend(lgd, loc="lower right",fontsize=fs)
plt.xticks(xticks, xtick_labels,fontsize=fs)
plt.xlim([0,1825])
plt.xlabel('Time [day]',fontsize=fs)
plt.ylabel('Water Table [m]',fontsize=fs)
plt.text(5,10.6,'(a) Water table at domain center',fontsize=fs)


plt.subplot(2,1, 2)
ax = fig.gca()
pos1 = ax.get_position()
pos2 = [pos1.x0-0.01, pos1.y0, pos1.width*1.1, pos1.height*1.05]
ax.set_position(pos2)
for ff in range(len(fdir)):
	plt.plot(x, wt_out[ff], color=color1[ff], linestyle=ls[ff])	
plt.scatter(hydrus1[:,0], hydrus1[:,1],  marker='^', facecolor='None', edgecolor='k')
plt.xlim([0,L])
plt.xticks([0,1000,2000,3000,4000],[0,1000,2000,3000,4000],fontsize=fs)
plt.xlabel('X [m]',fontsize=fs)
plt.ylabel('Water Table [m]',fontsize=fs)
plt.text(15,1,'(b) Water table at 5 years',fontsize=fs)

#plt.savefig('t4-results.eps',format='eps')
#plt.savefig('both.png',format='png',bbox_inches='tight',dpi=600)
plt.show()

