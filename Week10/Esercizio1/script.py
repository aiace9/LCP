import	numpy as np
import matplotlib.pyplot as plt
import os, sys
import glob

if __name__ == '__main__':
	print 'script partito'
	path = '*.dat'
	files = glob.glob(path)
	D_i = []
	print files
	for file in files:
		D_i.append(np.loadtxt(file))
	
	#inizio calcolo nuove quantita:
	time_steps= D_i[1].shape[0]
	simulations= len(D_i)
	print 'caricamento file completato, sono stati individuati', simulations, 'files'

	D_avg = np.zeros([time_steps,3])

	for i in xrange(0,time_steps):
		cum = np.float32(0.0)
		cum2 = np.float32(0.0)
		for j in xrange(0,simulations):
			cum += D_i[j][i,1]
			cum2+= D_i[j][i,1]**2

		D_avg[i,0] = D_i[1][i,0]
		D_avg[i,1] = cum/(np.float32(simulations))
		D_avg[i,2] = np.sqrt((cum2/simulations) - (cum/simulations)**2)

	np.savetxt('finale.dat1',D_avg)

	plt.clf()
	#Plot delle dispersioni per ogni run
	for i in xrange(1,simulations):
		plt.plot(D_i[i][:,0],D_i[i][:,1],'k.')
	#plot della media dei run
	plt.plot(D_avg[:,0],D_avg[:,1],'r-')
	#plot delle varianze
	plt.plot(D_avg[:,0],D_avg[:,1]+D_avg[:,2],'b-')
	plt.plot(D_avg[:,0],D_avg[:,1]-D_avg[:,2],'b-')
	plt.show()

	print 'end'