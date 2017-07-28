#Domain [-300 300]^2 km
#Origin (41.3209371228N, 289.46309961W)
#Projection Lambert
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import h5py as hp
import scipy.ndimage.filters as filter

plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
#origin = 36.7964N,-120.822E
origin = [36.7964, -120.822]

f=hp.File('attFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 



tol =-0.25
for t in range(dim[1]):
    print t
    dx, dy = np.gradient(attdata[t,:,:],0.3,0.3)
    dxdx, dydx = np.gradient(dx,0.3,0.3)
    dxdy, dydy = np.gradient(dy,0.3,0.3) 

    dirdiv = np.empty([dim[1],dim[2]])
    mineig = np.empty([dim[1],dim[2]])

    for i in range(dim[1]):
        for j in range(dim[2]):
            eig = np.linalg.eig([[dxdx[i,j],dxdy[i,j]],[dydx[i,j],dydy[i,j]]])
            eigmin =  np.argmin(eig[0])
            dirdiv[i,j] = np.dot(eig[1][:,eigmin],[dx[i,j],dy[i,j]])
            mineig[i,j] = eig[0][eigmin]

    m = Basemap(width=2000000,height=2000000,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='h',area_thresh=10.,projection='lcc',\
        lat_1=35.,lat_0=origin[0],lon_0=origin[1])
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(20,50,2.5),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-135,-105,2.5),labels=[False,False,False,True])
    m.drawstates()
    x = np.linspace(0, m.urcrnrx, dim[1])
    y = np.linspace(0, m.urcrnry, dim[2])
    xx, yy = np.meshgrid(x, y)
    potridge = np.ma.masked_where(mineig>=tol,dirdiv)
    ridge = m.contour(xx, yy, np.transpose(potridge),levels =[0],colors='blue')
    m = 15*t
    h, minute = divmod(m,60)
    plt.title("FTLE, 6-{0}-2016 {1:02d}{2:02d} UTC".format(9+(4+h)//24, (4+h)%24, minute))
    plt.savefig('owens_attracting_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
    plt.clf()
