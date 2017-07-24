#Domain [-300 300]^2 km
#Origin (41.3209371228N, 289.46309961W)
#Projection Lambert
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py as hp

cdict = {'red':  [(0.0, 0.0000, 0.0000),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 1.0000, 1.0000)],
        'green': [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.5450, 0.5450)],
        'blue':  [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.0000, 0.0000)]}
plt.register_cmap(name='CyanOrange', data=cdict)

#origin = [41.3209371228, 289.46309961]
origin = [41.3209371228,-70.53690039]
# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS84 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(width=112000,height=96000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='f',area_thresh=0.,projection='lcc',\
            lat_1=35.,lat_0=origin[0],lon_0=origin[1])
m.drawcoastlines()
m.drawparallels(np.arange(40.6,42.2,0.2),labels=[True,False,False,False])
m.drawmeridians(np.arange(-71.1,-69.6,0.3),labels=[False,False,False,True])

f=hp.File('AttFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
#print attdata.shape
attmean = np.mean(attdata)
#print attmean
f=hp.File('RepFTLEOutput.mat','r')
repdata = f['F'][:,:,:]
f.close
#print repdata.shape
repmean = np.mean(repdata)
#print repmean

attdata = np.ma.masked_less(attdata,3*attmean)
repdata = -1*np.ma.masked_less(repdata,3*repmean)

attmax = np.max(attdata)
repmin = np.min(repdata)

x = np.linspace(0, m.urcrnrx, attdata.shape[1])
y = np.linspace(0, m.urcrnry, attdata.shape[2])

xx, yy = np.meshgrid(x, y)

attquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),shading='gouraud',vmin=repmin,vmax=attmax,cmap='CyanOrange')
repquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(repdata[0,:,:])),shading='gouraud',vmin=repmin,vmax=attmax,cmap='CyanOrange')
cbar = plt.colorbar()
cbar.set_label('Attractiveness, s^-1')
ttl = plt.title("FTLE, 7-21-17 1000 UTC/6am EDT + 0 hrs")


def init():
    attquad.set_array([])
    repquad.set_array([])
    ttl.set_text("")
    return attquad, repquad, ttl

def animate(t, *args):
    attquad.set_array(np.ravel(np.transpose(np.squeeze(args[2][t,:,:]))))
    repquad.set_array(np.ravel(np.transpose(np.squeeze(args[3][t,:,:]))))
    ampm = ['am', 'pm']
    #h = int(t)
    #minute = int(t%1*15)
    m = 15*t
    h, minute = divmod(m,60)
    ttl.set_text("FTLE, 7-{0}-17 {1:02d}{2:02d} UTC/{3:02d}:{2:02d}{4} EDT".format(21+(6+h)//24, (6+h)%24, minute, (2+h-1)%12+1, ampm[((2+h)//12)%2]))
    #plt.savefig('day3-frame-{:04d}.tif'.format(int(t*4+1)), transparent=False, bbox_inches='tight')
    return attquad, repquad, ttl
    

myargs = (xx,yy,attdata,repdata)
anim = animation.FuncAnimation(plt.gcf(),animate,fargs=myargs,frames=attdata.shape[0],repeat=False)
#plt.savefig('day3-frame-{:04d}.tif'.format(int(t*4+1)), transparent=False, bbox_inches='tight')
#plt.close()
Writer = animation.writers['ffmpeg']
writer = Writer(fps=4)
anim.save('FTLE.mp4', writer=writer)
#plt.show()

