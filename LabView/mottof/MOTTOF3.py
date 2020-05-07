fitaxis = 0 # 0: x-axis, 1: y-axis

cropx = 770     # center x position of the crop square 
cropy = 500     # center y position of the crop square
cropsize = 350  # Half of crop square dimension

rcropx = 460    # Crop for background intensity comparison
rcropy = 480
rcropsize = 100 # Half of rcrop square dimension

# IMPORT STATEMENTS
import numpy
from PIL import Image, ImageFilter
from pylab import imshow, show, gray
from matplotlib import pyplot
from matplotlib import mlab
from scipy.optimize import curve_fit
from numpy import linspace

# to be used in offsetting background from MOT image
# you can tune them if you want, but I usually leave them at 0
delx = 0
dely = 0

# Takes in MOT image and background image, and outputs the sigma value
# of the log of intensity (along 1-Dimensional cross section)
def image_to_sigma(imgback,imgfore):
    # Crop the images and take their difference
    ratio = numpy.sum(imgback[rcropy-rcropsize:rcropy+rcropsize,rcropx-rcropsize:rcropx+rcropsize])/numpy.sum(imgfore[rcropy-rcropsize+dely:rcropy+rcropsize+dely,rcropx-rcropsize+delx:rcropx+rcropsize+delx])
    imdif = numpy.array(imgback[cropy-cropsize:cropy+cropsize,cropx-cropsize:cropx+cropsize] - ratio*imgfore[cropy-cropsize+dely:cropy+cropsize+dely,cropx-cropsize+delx:cropx+cropsize+delx],'float')
    (X,Y) = numpy.shape(imdif)
    
    # find center of brightness
    # loop to find center of image
    thresh =25 
    m = numpy.zeros((X,Y))
    
    for x in range(X):
        for y in range(Y):
                m[x, y] = imdif[x, y] >= thresh
    m = m / numpy.sum(numpy.sum(m))
        
    # marginal distributions
    dx = numpy.sum(m, 1)
    dy = numpy.sum(m, 0)
    
    # expected values
    cx = numpy.sum(dx * numpy.arange(X))
    cy = numpy.sum(dy * numpy.arange(Y))
        
        
    # take horizontal and vertical cuts of I_bck - I_mot
    # and also I_bck
    imref_full = numpy.array(imgback[cropy-cropsize:cropy+cropsize,cropx-cropsize:cropx+cropsize],'float')
    if fitaxis == 0:
        imref_cut = numpy.array(imref_full[int(cx),:], 'float')
        samp = numpy.array(imdif[int(cx),:],'float')
    elif fitaxis == 1:
        imref_cut = numpy.array(imref_full[:,int(cy)], 'float')
        samp = numpy.array(imdif[:,int(cy),],'float')
    
    x = numpy.arange(len(samp))
    z = numpy.zeros(len(samp))
        
    # take log only at pixel values where log is well defined
    for i in range(len(samp)):
        if samp[i]/imref_cut[i]< 1 and imref_cut[i] > 0 and samp[i] >0:
            z[i] = -numpy.log(1.0000-samp[i]/imref_cut[i])

    # xx = numpy.mgrid[0:X+0.1:1, 0:Y+0.1:1].reshape(2,-1).T
    zz = numpy.zeros([X,Y])
    for i in range(X):
        for j in range(Y):
            if imdif[i,j]/imref_full[i,j] and imdif[i,j] > 0:
                zz[i,j] = -numpy.log(1.0000-imdif[i,j]/imref_full[i,j])
    zz1d = zz.ravel()
                 
                    
    # find average and sigma
    mean = sum(x * z) / sum(z)
    sigma = numpy.sqrt(sum(z * (x - mean)**2) / sum(z))
    ambient = 0.01
    
    # define curve to be fitted
    def Gauss2D(XX, a, x0, y0, sigma, b):
        val = a * numpy.exp(-((XX[0]-x0)**2 + (XX[1]-y0)**2) / (2*sigma**2)) + b
        return val.ravel()
          
    # create x and y indices
    xx = numpy.linspace(0, X-1, X)
    yy = numpy.linspace(0, Y-1, Y)
    xx, yy = numpy.meshgrid(xx, yy)
    popt,pcov = curve_fit(Gauss2D, (xx,yy), zz1d, p0=[numpy.amax(zz), cx, cy, sigma, ambient])
        
    print('sigma = {0:.6f},\t popt[3] = {1:.6f}'.format(sigma, popt[3]))
    return popt[3]

# IMAGE IMPORT STAGE
filelist = ['05000us_0.jpg','07000us_0.jpg','09000us_0.jpg','11000us_0.jpg','13000us_0.jpg','15000us_0.jpg','17000us_0.jpg','19000us_0.jpg','21000us_0.jpg']


# SIGMA VALUE FINDING STAGE
sig_ar = numpy.zeros(len(filelist)) # Array to store the sigma values
t_ar = numpy.zeros(len(filelist)) # time array
bgarray = numpy.array(Image.open('background.jpg').convert('L'))
i = 0
for filename in filelist:
    imgarray = numpy.array(Image.open(filename).convert('L'))
    sigma = image_to_sigma(bgarray,imgarray)
    sig_ar[i] = sigma
    t_ar[i] = float(filename[0:5])*10**-6
    i += 1
    del imgarray
    #pyplot.plot(3500 + 500*i,sigma,'ko')
# print('SIGMA VALS:')
# print(sig_ar)

# LINEAR REGRESSION/TEMPERATURE OUTPUT STAGE
# sig_ar *= 2.5*0.0053*10**-3 # convert px number to meter
sig_ar *= 13.7*10**-6 # convert px number to meter

def fit_func(x, s, k): # function to fit to sigma array
    return numpy.sqrt(s**2 + k*x**2)

popt, pcov = curve_fit(fit_func, t_ar, sig_ar) # curve fitting
print(popt[0])


# Plot fit results
v_ar = numpy.linspace(min(t_ar),max(t_ar),1000)
yar = fit_func(v_ar,popt[0],popt[1])
pyplot.plot(1000*t_ar,sig_ar,'ko',ms=5)
pyplot.plot(1000*v_ar,yar,'b-',lw=3)
pyplot.xlabel('time [ms]', fontsize=18)
pyplot.ylabel('$\sigma$', fontsize=18)
pyplot.xticks(numpy.arange(0,26,2),fontsize=18)
pyplot.yticks([])
pyplot.tight_layout()
pyplot.show()

print('Temperature: {0}'.format((popt[1]*2.206948425*10**-25)/(1.38*10**-23)))