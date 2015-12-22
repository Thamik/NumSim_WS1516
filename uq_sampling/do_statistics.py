import numpy as np

from uq_sampling import Sample
from uq_statistics import Statistics, SimData

if __name__ == '__main__':
        # only run when executed as a script
        stats = Statistics()
        stats.loadSimData('VTK')
        stats.setTimes(np.linspace(0,30,100))
        stats.compute()
        stats.plot()
