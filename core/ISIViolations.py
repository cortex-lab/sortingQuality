def ISIViolations(spikeTrain, minISI, refDur):
	# See matlab version for commentary
	isis = np.diff(spikeTrain)
	nSpikes = len(spikeTrain)
	numViolations = sum(isis<refDur) 
	violationTime = 2*nSpikes*(refDur-minISI)
	totalRate = nSpikes/spikeTrain[-1]
	violationRate = numViolations/violationTime			
	fpRate = violationRate/totalRate
	return fpRate, numViolations