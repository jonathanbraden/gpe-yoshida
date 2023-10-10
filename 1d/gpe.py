class GPEData(object):
    def __init__(self, dataFile, logFile, trapFile, nLat):
        self.grid = Grid_1D(logFile, trapFile)
        self.nt = self.grid.nt
        self.nx = self.grid.nx
        self.dat = np.fromfile(dataFile).reshape((nt,2,2,nx))

        self.psi = self.dat[:,:,0] + 1j*self.dat[:,:,1]
        self.rho = np.sum(self.dat**2,axis=-2)
        
        return

    def saveHDF5(self):
        return

class Grid_1D(object):
    def __init__(self, logFile, trapFile):
        d = np.loadtxt(trapFile, usecols=[0,2])
        xVals = d[:,0]
        weights = d[:,1]

        tVals = np.loadtxt(logFile)[:,0]

        nt = tVals.shape[0]
        nx = xVals.shape[0]
        return
