from __future__ import division
import numpy, scipy.sparse, scipy.sparse.linalg, time
import MBIlib


class MBI(object):

    def __init__(self, P, xs, ms0=None, ks0=None):
        nx = len(xs)
        nf = 1 if len(P.shape) is nx else P.shape[nx]

        ns = numpy.array(P.shape[:nx],order='F')
        ms = numpy.array([int(ns[i]/3) for i in range(nx)])
        ks = 4*numpy.ones(ns.shape[0],int)
        ms[:] = ms[:] if ms0 is None else ms0
        ks[:] = ks[:] if ks0 is None else ks0
        for i in range(nx):
            ks[i] = min(ks[i], ms[i])
        nP = numpy.prod(ns)

        self.xs = xs
        self.ns, self.ms, self.ks = ns, ms, ks
        self.nx, self.nf, self.nP = nx, nf, nP

        ts = []
        for i in range(nx):
            k, m, n = ks[i], ms[i], ns[i]
            d = MBIlib.knotopen(k, k+m)
            ts.append(MBIlib.paramuni(k+m, m, n, d))

        t = numpy.zeros(P.shape[:nx]+(nx,),order='F')
        for ind,x in numpy.ndenumerate(t):
            t[ind] = ts[ind[-1]][ind[ind[-1]]]
        self.t = t.reshape((nP,nx),order='F')
        P = numpy.reshape(P,(nP,nf),order='F')

        B = self.assembleJacobian(0, 0)
        BT = B.transpose()
        BTB = BT.dot(B)
        BTP = BT.dot(P)

        nC = numpy.prod(ms)
        C = numpy.zeros((nC,nf),order='F')
        for i in range(nf):
            C[:,i] = scipy.sparse.linalg.cg(BTB,BTP[:,i])[0]

        Cx = []
        for i in range(nx):
            k, m, n = ks[i], ms[i], ns[i]
            nB = k*n
            Ba, Bi, Bj = MBIlib.computejacobian(0, 0, 1, n, nB, k, m, ts[i])
            B = scipy.sparse.csc_matrix((Ba,(Bi,Bj)))
            BT = B.transpose()
            BTB = BT.dot(B)
            BTP = BT.dot(xs[i])
            Cx.append(scipy.sparse.linalg.cg(BTB,BTP)[0])
            Cx[-1][0] = xs[i][0]
            Cx[-1][-1] = xs[i][-1]

        self.C, self.Cx = C, Cx

    def assembleJacobian(self, d1, d2):
        nx, nf, nP = self.nx, self.nf, self.nP
        ns, ms, ks = self.ns, self.ms, self.ks
        nB = nP*numpy.prod(ks)
        Ba, Bi, Bj = MBIlib.computejacobian(d1, d2, nx, nP, nB, ks, ms, self.t)
        return scipy.sparse.csc_matrix((Ba,(Bi,Bj)))

    def evaluate(self, x, d1=0, d2=0):
        getf = lambda d1, d2: MBIlib.evaluatept(d1, d2, nx, nf, self.C.shape[0], ks, ms, t, self.C)
        getx = lambda i, d1, d2: MBIlib.evaluatept(d1, d2, 1, 1, self.Cx[i].shape[0], ks[i], ms[i], t[i], self.Cx[i])

        nx, nf, nP = self.nx, self.nf, self.nP
        ns, ms, ks = self.ns, self.ms, self.ks

        t = numpy.zeros(nx)
        for i in range(nx):
            k, m, n = ks[i], ms[i], ns[i]
            t[i] = MBIlib.inversemap(k, m, x[i], self.Cx[i])

        if d1 is 0 and d2 is 0:
            return getf(0,0)
        elif d1 is not 0 and d2 is 0:
            return getf(d1,0)/getx(d1-1,1,0)
        elif d1 is not 0 and d1 is d2:
            return (getf(d1,d1)*getx(d1-1,1,0) - getf(d1,0)*getx(d1-1,1,1))/getx(d1-1,1,0)**3
        elif d1 is not 0 and d2 is not 0 and d1 is not d2:
            return getf(d1,d2)/getx(d1-1,1,0)/getx(d2-1,1,0)
