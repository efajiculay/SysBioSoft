from numpy import linspace, interp

def sample_points(edata,nsamp=50000):
	data, Si = edata
	tn = round(data[0][:,0][-1])
	t = linspace(0,tn,nsamp)
	sdata = [t]	
	for dd in data:
		d = []
		n = dd.shape[1]
		for i in range(1,n):
			y = interp(t, dd[:,0], dd[:,i])
			d.append(y)
		sdata.append(d)
	return (sdata, Si)