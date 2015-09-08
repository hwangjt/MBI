import unittest

import MBI
import numpy
from numpy.testing import assert_allclose
import time

class MBITestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_basic(self):
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

        #print time.time() - t0
        #print '---f'
        #print f0
        
        # print '---df/dx'
        assert_allclose(m.evaluate(P,1,0), (m.evaluate(P+h*e1) - f0)/h, rtol=1.e-4)
        # print '---df/dy'
        assert_allclose(m.evaluate(P,2,0), (m.evaluate(P+h*e2) - f0)/h, rtol=1.e-4)
        # print '---d2f/dx2'
        assert_allclose(m.evaluate(P,1,1), 
                       (m.evaluate(P+h*e1) - 2*f0 + m.evaluate(P-h*e1))/h**2, rtol=1.e-3)
        # print '---d2f/dxdy'
        assert_allclose(m.evaluate(P,1,2), 
                       (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - 
                        m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2, rtol=1.e-4)
        # print '---d2f/dydx'
        assert_allclose(m.evaluate(P,2,1), 
                       (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - 
                        m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2, rtol=1.e-4)
        # print '---d2f/dy2'
        assert_allclose(m.evaluate(P,2,2), 
                       (m.evaluate(P+h*e2) - 2*f0 + m.evaluate(P-h*e2))/h**2, rtol=1.e-3)

        # print '---df/dx'
        # print m.evaluate(P,1,0)
        # print (m.evaluate(P+h*e1) - f0)/h
        # print '---df/dy'
        # print m.evaluate(P,2,0)
        # print (m.evaluate(P+h*e2) - f0)/h
        # print '---d2f/dx2'
        # print m.evaluate(P,1,1)
        # print (m.evaluate(P+h*e1) - 2*f0 + m.evaluate(P-h*e1))/h**2
        # print '---d2f/dxdy'
        # print m.evaluate(P,1,2)
        # print (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2
        # print '---d2f/dydx'
        # print m.evaluate(P,2,1)
        # print (m.evaluate(P+h*e1+h*e2) - m.evaluate(P+h*e1-h*e2) - m.evaluate(P-h*e1+h*e2) + m.evaluate(P-h*e1-h*e2))/4.0/h**2
        # print '---d2f/dy2'
        # print m.evaluate(P,2,2)
        # print (m.evaluate(P+h*e2) - 2*f0 + m.evaluate(P-h*e2))/h**2

if __name__ == '__main__':
    unittest.main()

