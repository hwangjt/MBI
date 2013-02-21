import MBI
import numpy, time


nx = 12
ny = 12
x = numpy.linspace(0,1,nx)
y = numpy.linspace(0,1,ny)
P = numpy.array(numpy.meshgrid(x,y))
P = numpy.swapaxes(P,0,1)
P = numpy.swapaxes(P,1,2)
m = MBI.MBI(P,[x,y],[5,5],[4,4])

h = 1e-5
t0 = time.time()
f0 = m.evaluate([0.25,0.25])
print time.time() - t0
print f0
print '---df/dx'
print m.evaluate([0.25,0.25],1,0)
print (m.evaluate([0.25+h,0.25]) - f0)/h
print '---df/dy'
print m.evaluate([0.25,0.25],2,0)
print (m.evaluate([0.25,0.25+h]) - f0)/h
print '---d2f/dx2'
print m.evaluate([0.25,0.25],1,1)
print (m.evaluate([0.25+h,0.25]) - 2*f0 + m.evaluate([0.25-h,0.25]))/h**2
print '---d2f/dxdy'
print m.evaluate([0.25,0.25],1,2)
print (m.evaluate([0.25+h,0.25+h]) - m.evaluate([0.25+h,0.25-h]) - m.evaluate([0.25-h,0.25+h]) + m.evaluate([0.25-h,0.25-h]))/4.0/h**2
print '---d2f/dydx'
print m.evaluate([0.25,0.25],2,1)
print (m.evaluate([0.25+h,0.25+h]) - m.evaluate([0.25+h,0.25-h]) - m.evaluate([0.25-h,0.25+h]) + m.evaluate([0.25-h,0.25-h]))/4.0/h**2
print '---d2f/dy2'
print m.evaluate([0.25,0.25],2,2)
print (m.evaluate([0.25,0.25+h]) - 2*f0 + m.evaluate([0.25,0.25-h]))/h**2
