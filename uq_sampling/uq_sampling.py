import subprocess
import random
import time

class SimFrame:
	
        def __init__(self, no_samples, no_processes):
                self.__no_samples = no_samples
                self.__no_processes = no_processes

                if self.__no_samples <= 0:
                        print 'Warning: number of samples smaller or equal to zero. Defaulting to 1.'
                        self.__no_samples = 1

                if self.__no_processes <= 0:
                        print 'Warning: number of processes smaller or equal to zero. Defaulting to 1.'
                        self.__no_processes = 1

                if self.__no_samples < self.__no_processes:
                        print 'Warning: number of samples smaller than number of processes. Defaulting.'
                        self.__no_processes = self.__no_samples

        def generateSamples(self):
                self.__sampler = Sampler(self.__no_samples)
                self.__sampler.setNormalDistribution(2000.0, 1000.0/6.0)
                #self.__sampler.setTrapezoidalRule(2000.0, 1000.0/6.0)
                self.__sampler.computeSamples()

        def run(self):
                pool = Pool(self.__no_processes)
                
                while self.__sampler.hasMoreSamples():
                        s = self.__sampler.nextSample()
                        while not pool.hasFreeSlot():
                                time.sleep(0.1) # wait a bit
                        cmds = [ (['./geomgen_release', '16384', '0', str(s.re), str(s.sid)], 'geom_gen'),
                                (['./numsim', 'uq_data/uq_parameter_' + str(s.sid) + '.params', 'uq_data/uq_geometry.geom', str(s.sid)], '.') ]
                        #print "Run process: " + str(cmds)
                        print "Running sample %s/%s" % (str(s.n+1), self.__no_samples)
                        pool.runCommands(cmds)
                
                # wait until all processes returned
                while not pool.allFinished():
                        time.sleep(0.1)
                
                # everything is done
                print 'Running all samples finished.'

class Pool:

        def __init__(self, no_procs):
                self.__no_procs = no_procs
                self.__procs = []

        def runCommand(self, command, cwd=None):
                self.__procs.append(subprocess.Popen(command, cwd=cwd))

        def runCommands(self, cmds):
                for ii, data in enumerate(cmds):
                        cmd, cwd = data
                        p = subprocess.Popen(cmd, cwd=cwd)
                        if not ii == len(cmds)-1:
                                while p.poll() is None:
                                        time.sleep(0.1)
                # dont wait for the last command, just add it to the list
                self.__procs.append(p)

        def __numberOfFinishedProcesses(self):
                return sum([ 1 if p.poll() is not None else 0 for p in self.__procs ])

        def __removeFinishedProcesses(self):
                while self.__numberOfFinishedProcesses() > 0:
                        for p in self.__procs:
                                if p.poll() is not None:
                                        self.__procs.remove(p)
                                        break

        def hasFreeSlot(self):
                self.__removeFinishedProcesses()
                return len(self.__procs) < self.__no_procs

        def allFinished(self):
                self.__removeFinishedProcesses()
                return len(self.__procs) == 0


class Sampler:

        def __init__(self, no_samples):
                self.__no_samples = no_samples
		self.__samples = []

        def setNormalDistribution(self, mu, sigma):
                self.__mu = mu
                self.__sigma = sigma
                self.__strategy = 'normaldist'

        def setTrapezoidalRule(self, mu, sigma):
                self.__mu = mu
                self.__sigma = sigma
                self.__strategy = 'trapezoidal'

        def computeSamples(self):
                if self.__strategy == 'normaldist':
                        random.seed()
                #self.__samples = []
                for ii in range(self.__no_samples):
                        if self.__strategy == 'normaldist':
                                re = random.normalvariate(self.__mu, self.__sigma)
                                weight = 1.0
                        elif self.__strategy == 'trapezoidal':
                                re = self.__mu-3.0*self.__sigma + 6.0*self.__sigma * float(ii)/float(self.__no_samples-1)
                                weight = 1.0 # not sure if this is right
                        sample_id = random.randint(100000000,999999999)
                        self.__samples.append(Sample(ii, sample_id, re, weight))

        def hasMoreSamples(self):
                return len(self.__samples)>0

        def nextSample(self):
                return self.__samples.pop(0)

        def numberOfRemainingSamples(self):
                return len(self.__samples)

        def __str__(self):
                return self.__samples.__str__()


class Sample:

        def __init__(self, number, sample_id, reynolds_number, weight=1.0):
                self.n = number
                self.sid = sample_id
                self.re = reynolds_number
                self.w = weight

	def __str__(self):
		return (self.sid, self.re, self.w).__str__()

