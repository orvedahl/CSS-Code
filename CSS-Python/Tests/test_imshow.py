
import pylab
import numpy
import image_polar as myim
import curved_polar as my_im

test = 2

#a=outer(arange(0,1,0.01),ones(10)).transpose()
#a=numpy.outer(numpy.arange(0,1,0.01),numpy.arange(0,1,0.01)[::-1]).transpose()

r = numpy.arange(6.e10,7.e10,0.01e10)
theta = numpy.arange(80., 100., 0.1)

a = numpy.zeros((len(theta), len(r)))
for i in range(len(theta)):
    th = theta[i]*numpy.pi/180.
    for j in range(len(r)):
        rad = r[j]
        a[i,j] = rad*numpy.sin(th)

if test==1:
    pylab.figure()
    im = pylab.imshow(a, aspect='auto', cmap=pylab.get_cmap("jet"), 
                      origin="upper", extent=(numpy.amin(r), numpy.amax(r), 
                      numpy.amax(theta), numpy.amin(theta)))
    pylab.colorbar(im)
    pylab.ylabel("jet", rotation=90)
    pylab.show()
elif test==2:
    theta = theta*numpy.pi/180.
    #myim.image_polar(a, r, numpy.cos(theta), numpy.sin(theta), 
                            #r_bcz=[0.98], cont=True, cbar=True, add_c=False,
                            #tit='Hi', mini=6.2e10, maxi=6.8e10)
    my_im.polar_image(a, r, numpy.cos(theta), numpy.sin(theta), 
                 r_bcz=[6.75e10,6.85e10], equator=True, vmin=None, vmax=None,
                 cont=True, cmap=pylab.get_cmap("jet"), cbar=True, 
                 add_cont=False, add_data=None, cb_title="Hi")
    pylab.title("Hello")
    pylab.show()
