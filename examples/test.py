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

P = numpy.array([[0.25,0.25],[0.25,0.25],[0.25,0.25]])
e1 = numpy.zeros((3,2))
e2 = numpy.zeros((3,2))
e1[:,0] = 1.0
e2[:,1] = 1.0

h = 1e-5
t0 = time.time()
f0 = m.evaluate(P)
print time.time() - t0
print '---f'
print f0
print '---df/dx'
print m.evaluate(P,1,0)
print (m.evaluate(P+h*e1) - f0)/h
print '---df/dy'
print m.evaluate(P,2,0)
print (m.evaluate(P+h*e2) - f0)/h
print '---d2f/dx2'
print m.evaluate(P,1,1)
print (m.evaluate(P+h*e1) - 2*f0 + m.evaluate(P-h*e1))/h**2
print '---d2f/dxdy'
print m.evaluate(P,1,2)
print (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2
print '---d2f/dydx'
print m.evaluate(P,2,1)
print (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2
print '---d2f/dy2'
print m.evaluate(P,2,2)
print (m.evaluate(P+h*e2) - 2*f0 + m.evaluate(P-h*e2))/h**2
