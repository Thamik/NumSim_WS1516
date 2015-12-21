import subprocess
import random

class SimFrame:
	
        def __init__(self, no_samples, no_processes):
                self.__no_samples = no_samples
                self.__no_processes = no_processes

                if self.__no_samples <= 0:
                        print 'Warning: number of samples smaller or equal to zero. Defaulting to 1.'
                        self.__no_samples = 1

                if self.__no_processes < self.__no_samples:
                        print 'Warning: number of processes smaller than number of samples.'
                        self.__no_processes = self.__no_samples

        def run():
                self.__procs = []


class Sampler:

        def __init__(self, no_samples):
                self.__no_samples = no_samples

        def setNormalDistribution(mu, sigma):
                self.__mu = mu
                self.__sigma = sigma
                self.__strategy = 'normaldist'

        def setTrapezoidalRule(mu, sigma):
                self.__mu = mu
                self.__sigma = sigma
                self.__strategy = 'trapezoidal'

        def getSamples():
                if self.__strategy == 'normaldist':
                        random.seed()
                samples = []
                for ii in range(self.__no_samples):
                        if self.__strategy == 'normaldist':
                                re = random.normalvariate(self.__mu, self.__sigma)
                                weight = 1.0
                        elif self.__strategy == 'trapezoidal':
                                re = self.__mu-3.0*self.__sigma + 6.0*self.__sigma * float(ii)/float(self.__no_samples-1)
                                weight = 1.0 # not sure if this is right
                        sample_id = random.randint(1000000,9999999)
                        samples.append(Sample(sample_id, re, weight))
                return samples


class Sample:

        def __init__(self, sample_id, reynolds_number, weight=1.0):
                self.sid = sample_id
                self.re = reynolds_number
                self.w = weight

