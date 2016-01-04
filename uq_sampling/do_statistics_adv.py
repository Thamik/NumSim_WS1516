import numpy as np

from uq_sampling import Sample
from uq_statistics_adv import Statistics, SimData

if __name__ == '__main__':
        # only run when executed as a script
        stats = Statistics()
        stats.loadSimData('VTK')
        stats.setTimes(np.linspace(0,30,100))
	stats.setGrid(np.linspace(1000,2000,200))
        stats.compute()
	stats.computeTrap()
	stats.computeConvergence()
        stats.plot()
