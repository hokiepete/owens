from netCDF4 import Dataset
import numpy


ncfile = "../Downloads/nam.t06z.conusnest.hiresf00.tm00.nc"

root = Dataset(ncfile,'r')
vars = root.variables

gridspacing = 3
lat = vars["latitude"][:,:]
lon = vars["longitude"][:,:]
dim = lon.shape
print dim
latorg = 41.325 
lonorg = -70.5225+360
tol = 0.015
#Find the index of the origin
print "want lat lon", latorg, lonorg
for i in range(dim[0]):
    for j in range(dim[1]):
        #print lat[i,j], lon[i,j]
        if numpy.fabs(lat[i,j]-latorg)<tol and numpy.fabs(lon[i,j]-lonorg)<tol:
            print 'closest we can get is', lat[i,j], lon[i,j]
            print i,j
            iorg = i
    	    jorg = j

#hardcoded origin to save time
#iorg = 195
#jorg = 155
#choose 42, want 1000km domain centered at (36.75,-120.875), with 12.191km gridspacing, 42 is roughly 500km
u = 3.6*numpy.squeeze(vars["UGRD_10maboveground"][:,:,:])
v = 3.6*numpy.squeeze(vars["VGRD_10maboveground"][:,:,:])
umax = numpy.max(u)
vmax = numpy.max(v)
print "Max u ", umax
print "Max v ", vmax
space0 = gridspacing*(dim[0]-iorg)
space1 = gridspacing*(dim[1]-jorg)
print "Have " + str(space0) + " kms of space from origin N/S"
print "Have " + str(space1) + " kms of space from origin E/W"
time0 = space0/vmax
time1 = space1/umax
print "Can integrate for " + str(numpy.min([time0,time1])) + "hrs"
#print time1
'''
def stressToKmphArray(tauX,tauY,scalefactor,c0):
    dim = numpy.shape(tauX)
    for i in range(dim[0]):
        for j in range(dim[1]):
            tauMag = numpy.sqrt(tauX[i,j]*tauX[i,j] + tauY[i,j]*tauY[i,j])
            if tauMag > c0:
                mpsX[i,j] = tauX[i,j]/numpy.sqrt(scalefactor*tauMag[i,j])
                mpsY[i,j] = tauY[i,j]/numpy.sqrt(scalefactor*tauMag[i,j])
            else:
                mpsX[i,j] = c0
                mpsY[i,j] = c0
            kmphX[i,j] = 3.6*mpsX[i,j]
            kmphY[i,j] = 3.6*mpsY[i,j]

    return (kmphX, kmphY);

main()
'''
