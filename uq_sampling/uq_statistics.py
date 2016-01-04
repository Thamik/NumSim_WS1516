import os
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

from uq_sampling import Sample

class Statistics:

        def __init__(self):
                self.__data = []

        def loadSimData(self, datapath='VTK'):
                print 'Loading simulation data...'
                simFiles = [ (datapath+'/'+path) for path in os.listdir(datapath) if path.startswith('uq_data') ]
                for filename in simFiles:
                        self.__loadFile(filename)
                print 'Loading simulation data finished. Total number of samples: ' + str(len(self.__data))

        def __loadFile(self, filename):
                f = open(filename, 'r')
                sid = int(f.readline())
                re = float(f.readline())
                s = SimData(sid, re)
                
                line = f.readline()
                while not line == '':
                        t, u1, u2, u3 = map(float, line.split())
                        s.data[t] = (u1,u2,u3)
                        line = f.readline()

                self.__data.append(s)

        def setTimes(self, times):
                self.__times = sorted(times)
                # check on the given times
                if any([ t<0 or t>30 for t in times]):
                        # unvalid times
                        raise AttributeError()
                # interpolate
                print 'Interpolating data to match the given times...'
                self.__interpolate_to_times()

        def compute(self):
                print 'Computing mean values and standard deviations over time...'
                self.__means = []
                self.__stds = []
                for ii, t in enumerate(self.__times):
                        mean = [0.0] * len(self.__data_interp[ii][0])
                        for jj in range(len(self.__data_interp[ii])):
                                for kk in range(len(self.__data_interp[ii][0])):
                                        mean[kk] += self.__data_interp[ii][jj][kk]
                        for kk in range(len(self.__data_interp[ii][0])):
                                mean[kk] /= len(self.__data_interp[ii])

                        std = [0.0] * len(self.__data_interp[ii][0])
                        for jj in range(len(self.__data_interp[ii])):
                                for kk in range(len(self.__data_interp[ii][0])):
                                        std[kk] += (self.__data_interp[ii][jj][kk] - mean[kk]) ** 2
                        for kk in range(len(self.__data_interp[ii][0])):
                                std[kk] /= len(self.__data_interp[ii]) - 1
                                std[kk] = sqrt(std[kk])

                        self.__means.append(mean)
                        self.__stds.append(std)
				

        def __interpolate_to_times(self):
                self.__data_interp = [ [ s.interpolate(t) for s in self.__data ] for t in self.__times ]

        def plot(self):
                print 'Plotting...'
                
                # plot the data for the first evaluation point
                u1_mean = [ m[0] for m in self.__means ]
                u1_std = [ s[0] for s in self.__stds ]
                
                plt.figure(1)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u1_mean))
                plt.ylabel('Mean value')
                plt.title('Velocity at the grid cell (120,5)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u1_std))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

                # plot the data for the second evaluation point
                u2_mean = [ m[1] for m in self.__means ]
                u2_std = [ s[1] for s in self.__stds ]
                
                plt.figure(2)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u2_mean))
                plt.ylabel('Mean value')
                plt.title('Velocity at the grid cell (64,64)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u2_std))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

                # plot the data for the third evaluation point
                u3_mean = [ m[2] for m in self.__means ]
                u3_std = [ s[2] for s in self.__stds ]
                
                plt.figure(3)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u3_mean))
                plt.ylabel('Mean value')
                plt.title('Velocity at the grid cell (5,120)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u3_std))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

		# the distribution of the reynolds numbers of the samples
		res = [ s.re for s in self.__data ]

		plt.figure(4)
		#n, bins, patches = plt.hist(res, 10, normed=1)
                n, bins, patches = plt.hist(res, 25, normed=1)
		mu = 1500.0
		sigma = 1000.0/6.0
		y = mlab.normpdf(bins, mu, sigma)
		plt.plot(bins, y, 'r--')

                plt.show()


class SimData:

        def __init__(self, sid, re):
                self.sid = sid
                self.re = re
                self.data = {}
                self.__interp_prepared = False

        def prepareInterpolation(self):
                self.__times = sorted(self.data.keys())

                self.__maxtime = max(self.__times)
                self.__mintime = min(self.__times)
#                if self.__maxtime < 30:
#                        print "Faulty simulation data: too small end time: " + str(self.__maxtime) + "! (sid: " + str(self.sid) + ")"
#                if self.__mintime > 0:
#                        #print "Faulty simulation data: too large starting time: " + str(self.__mintime) + "! (sid: " + str(self.sid) + ")"

                self.__interp_prepared = True

        def interpolate(self, time):
                if not self.__interp_prepared:
                        self.prepareInterpolation()

                if time < self.__mintime or time > self.__maxtime:
#                        print 'Warning: interpolation time out of time range!'
                        return tuple([float('NaN')]*3)

                times = self.__times
                index = 0
                while times[index] < time:
                        index += 1
                t1 = times[index-1]
                t2 = times[index]
                w1 = (time-t1)/float(t2-t1)
                w2 = 1.0 - w1

                res = [0.0] * len(self.data[t1])
                for ii in range(len(self.data[t1])):
                        res[ii] = self.data[t1][ii] * w1 + self.data[t2][ii] * w2

                return tuple(res)

