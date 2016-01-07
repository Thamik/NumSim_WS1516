import os
from math import sqrt
from math import pi
from math import exp
from math import isnan
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

	def setGrid(self, grid):
		self.__grid = grid
		print 'Iterpolating data to match the given grid in probability space'
		self.__data_interp_grid = self.__interpolate_to_grid(self.__grid)

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

	def computeTrap(self):
		print 'Computing mean values and standard deviations over time with quadrature...'
		self.__expectQuad = []
		self.__stdsQuad = []
		for t in range(len (self.__times)):
			expect = [0.0] * 3
			for ii in range(len(self.__grid)):
				if (ii == range(len(self.__grid))[0] or ii == range(len(self.__grid))[-1]):
					# first and last element
					for jj in range(len(self.__data_interp_grid[ii].data[0])):
						expect[jj] += 0.5*self.__data_interp_grid[ii].data[t][jj]*self.prob(self.__grid[ii])
	
				else:
					for jj in range(len(self.__data_interp_grid[ii].data[0])):
						expect[jj] += self.__data_interp_grid[ii].data[t][jj]*self.prob(self.__grid[ii])
			
			for kk in range(len(self.__data_interp_grid[0].data[0])):
				h = (self.__grid[-1] - self.__grid[0])/float(len(self.__grid))
				expect[kk] *= h
				#print t, expect[kk]
			
			self.__expectQuad.append(expect)

		for t in range(len(self.__times)):
			std = [0.0] * 3
			for ii in range(len(self.__grid)):
				if (ii == range(len(self.__grid))[0] or ii == range(len(self.__grid))[-1]):
					# first and last element
					for jj in range(len(self.__data_interp_grid[ii].data[0])):
						std[jj] += 0.5*self.prob(self.__grid[ii])*(self.__data_interp_grid[ii].data[t][jj] - self.__expectQuad[t][jj])**2
	
				else:
					for jj in range(len(self.__data_interp_grid[ii].data[0])):
						std[jj] += self.prob(self.__grid[ii])*(self.__data_interp_grid[ii].data[t][jj] - self.__expectQuad[t][jj])**2
			
			for kk in range(len(self.__data_interp_grid[0].data[0])):
				h = (self.__grid[-1] - self.__grid[0])/float(len(self.__grid))
				std[kk] *= h
				std[kk] = sqrt(std[kk])

			self.__stdsQuad.append(std)

	def prob(self, reNum):
		sigma = 1000.0/6.0
		mu = 1500.0
		res = 1/(float(sigma)*sqrt(2*pi)) * exp(-0.5*((reNum-mu)/float(sigma))**2)
		return res

	def computeConvergence(self):
		# only do this all for the first evaluation point/cell
		print 'Computing mean and standard deviations at end time...'

		self.__meanConv = []
		temp_sum_mean = 0.0
		for jj in range(len(self.__data_interp[-2])):
			n = jj + 1
			temp_sum_mean += self.__data_interp[-2][jj][0]
			self.__meanConv.append( (n, temp_sum_mean/n) )
                
                self.__stdConv = []
		temp_sum_std = 0.0
		mean = temp_sum_mean/n
		for jj in range(len(self.__data_interp[-2])):
			n = jj + 1
			temp_sum_std += (self.__data_interp[-2][jj][0] - mean)**2
			self.__stdConv.append( (n, sqrt(temp_sum_std/(n-1))) )

		temp_mean_quad = self.__expectQuad
		temp_stds_quad = self.__stdsQuad
		temp_grid = self.__grid
		self.__meanConvQuad = []
		self.__stdConvQuad = []
		#for n in [ (ii+1)*10 for ii in range(50) ]:
		for n in [50, 100, 200]:
			self.setGrid(np.linspace(1000,2000,n))
			self.computeTrap()
			self.__meanConvQuad.append( (n, self.__expectQuad[-2][1]) )
			self.__stdConvQuad.append( (n, self.__stdsQuad[-2][1]) )

		self.__expectQuad = temp_mean_quad
		self.__stdsQuad = temp_stds_quad
		self.__grid = temp_grid
		

        def __interpolate_to_times(self):
                self.__data_interp = [ [ s.interpolate(t) for s in self.__data ] for t in self.__times ]

	def __interpolate_to_grid(self, grid):
		data_interp_grid = []
		reNumsSort = sorted(range(len(self.__data)), key=lambda k: self.__data[k])
		ind = 0
		for reCurr in grid:
			index = 0
			if (self.__data[reNumsSort[0]].re >= reCurr):
				index = 1
			elif (self.__data[reNumsSort[-1]].re < reCurr):
				index = len(reNumsSort) - 1
			else:
				while self.__data[reNumsSort[index]].re < reCurr:
					index += 1

			re1 = self.__data[reNumsSort[index-1]].re
			re2 = self.__data[reNumsSort[index]].re
			w1 = (reCurr-re1)/float(re2-re1)
			w2 = 1.0 - w1
			tempRes = SimData(ind, reCurr)
			for timeCurr in range(len(self.__times)):
				u1Weigh = self.__data_interp[timeCurr][index-1][0] * w1 + self.__data_interp[timeCurr][index][0] * w2
				u2Weigh = self.__data_interp[timeCurr][index-1][1] * w1 + self.__data_interp[timeCurr][index][1] * w2
				u3Weigh = self.__data_interp[timeCurr][index-1][2] * w1 + self.__data_interp[timeCurr][index][2] * w2
#				if (isnan(u1Weigh)):
#					print '(1) NaN found!!!: ', self.__data_interp[timeCurr][index-1][0], self.__data_interp[timeCurr][index][0]
#				elif (isnan(u2Weigh)):
#					print '(2) NaN found!!!: ', self.__data_interp[timeCurr][index-1][1], self.__data_interp[timeCurr][index][1]
#				elif (isnan(u3Weigh)):
#					print '(3) NaN found!!!: ', self.__data_interp[timeCurr][index-1][2], self.__data_interp[timeCurr][index][2]

				tempRes.data[timeCurr] = (u1Weigh, u2Weigh, u3Weigh)

			data_interp_grid.append(tempRes)
			ind += 1

		return data_interp_grid
			

        def plot(self):
                print 'Plotting...'

		# plot MC results
                
                # plot the data for the first evaluation point
                u1_mean = [ m[0] for m in self.__means ]
                u1_std = [ s[0] for s in self.__stds ]
                
                plt.figure(1)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u1_mean))
                plt.ylabel('Mean value')
                plt.title('MC (' + str(len(self.__data)) + ' samples): Velocity at the grid cell (120,5)')
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
                plt.title('MC (' + str(len(self.__data)) + ' samples): Velocity at the grid cell (64,64)')
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
                plt.title('MC (' + str(len(self.__data)) + ' samples): Velocity at the grid cell (5,120)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u3_std))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

		# the distribution of the reynolds numbers of the samples
		res = [ s.re for s in self.__data ]

		plt.figure(4)
		n, bins, patches = plt.hist(res, 10, normed=1)
		mu = 1500.0
		sigma = 1000.0/6.0
		y = mlab.normpdf(bins, mu, sigma)
		plt.plot(bins, y, 'r--')

		# plot quad results

                # plot the data for the first evaluation point
                u1_expQuad = [ m[0] for m in self.__expectQuad ]
                u1_stdQuad = [ s[0] for s in self.__stdsQuad ]
                
                plt.figure(5)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u1_expQuad))
                plt.ylabel('Mean value')
                plt.title('Quad (' + str(len(self.__grid)) + ' nodes): Velocity at the grid cell (120,5)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u1_stdQuad))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

		# plot the data for the second evaluation point
                u2_expQuad = [ m[1] for m in self.__expectQuad ]
                u2_stdQuad = [ s[1] for s in self.__stdsQuad ]
                
                plt.figure(6)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u2_expQuad))
                plt.ylabel('Mean value')
                plt.title('Quad (' + str(len(self.__grid)) + ' nodes): Velocity at the grid cell (64,64)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u2_stdQuad))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

		# plot the data for the third evaluation point
                u3_expQuad = [ m[2] for m in self.__expectQuad ]
                u3_stdQuad = [ s[2] for s in self.__stdsQuad ]
                
                plt.figure(7)
                plt.subplot(2,1,1)
                plt.plot(np.array(self.__times), np.array(u3_expQuad))
                plt.ylabel('Mean value')
                plt.title('Quad (' + str(len(self.__grid)) + ' nodes): Velocity at the grid cell (5,120)')
                plt.subplot(2,1,2)
                plt.plot(np.array(self.__times), np.array(u3_stdQuad))
                plt.ylabel('Standard deviation')
                plt.xlabel('Time')

		# plot the convergence of mean and std of u1 at the first evaluation point
		# unzip
		n_samples_means, means = map(list, zip(*self.__meanConv))
		n_samples_stds, stds = map(list, zip(*self.__stdConv))

		means_err = map(abs, map(lambda x: x-means[-1], means))
		stds_err = map(abs, map(lambda x: x-stds[-1], stds))

		plt.figure(8)
		plt.plot(np.array(n_samples_means), np.array(means_err))
		plt.ylabel('Mean value')
		plt.title('MC: Mean value at the first evaluation point')
		plt.xlabel('Number of Samples')

		plt.figure(9)
		plt.plot(np.array(n_samples_stds), np.array(stds_err))
		plt.ylabel('Standard deviation')
		plt.title('MC: Standard deviation at the first evaluation point')
		plt.xlabel('Number of Samples')

		# plot the convergence of mean and std of u1 at the first evaluation point
		# unzip
		n_nodes_means, means_quad = map(list, zip(*self.__meanConvQuad))
		n_nodes_stds, stds_quad = map(list, zip(*self.__stdConvQuad))

		meansQuad_err = map(abs, map(lambda x: x-means_quad[-1], means_quad))
		stdsQuad_err = map(abs, map(lambda x: x-stds_quad[-1], stds_quad))

		plt.figure(10)
		plt.plot(np.array(n_nodes_means), np.array(meansQuad_err))
		plt.ylabel('Mean value')
		plt.title('Quad: Mean value at the first evaluation point')
		plt.xlabel('Number of Nodes')

		plt.figure(11)
		plt.plot(np.array(n_nodes_stds), np.array(stdsQuad_err))
		plt.ylabel('Standard deviation')
		plt.title('Quad: Standard deviation at the first evaluation point')
		plt.xlabel('Number of Nodes')

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

	def __lt__(self, other):
		return self.re < other.re
