import sys
from mpi4py import MPI
import numpy as np
import time

comm = MPI.Comm.Get_parent()
myrank = comm.Get_rank()
nprocs = comm.Get_size()

#print comm
#print "proc %s here" % comm.Get_name()
#time.sleep(myrank)

Nrowstotal    = np.int(sys.argv[1])
Ncolumnstotal = np.int(sys.argv[2])
Ncolumnslocal = Ncolumnstotal/nprocs

if Ncolumnstotal%nprocs > myrank:
    Ncolumnslocal+=1

# receive angles and data
angles = np.zeros( (Nrowstotal,), dtype=np.float )
comm.Recv( angles, source=0, tag=myrank*11 )
data = np.zeros( (Nrowstotal, Ncolumnslocal), dtype=np.float )
comm.Recv( data, source=0, tag=myrank*123 )

#print "proc%d" % (myrank),
#print data


### this part is as before: ###

# how many data columns do we have?
Nspots = data.shape[1]

Nphases = 91
phases  = np.linspace( 0, 90, Nphases )

rm = residualmatrix    = np.zeros( (phases.size, Nspots) )
cm = coefficientmatrix = np.zeros( (phases.size, 2, Nspots) )

# init matrix for linear fit
er = np.ones( (angles.size,2) )

for pi,phase in enumerate( phases ):
    # write phase-shifted cos**2 into second column
    er[:,1] = np.cos( 2*(angles-phase)*np.pi/180.0 )
    # perform linear fit 
    f = np.linalg.lstsq( er, data ) 
    residualmatrix[pi,:] = f[1]
    cm[pi,:,:] = f[0][:]

mm = minindices = np.argmin( residualmatrix, axis=0 )
rp = resultingphases = phases[minindices]
I_0 = np.zeros( (Nspots,) )
M_0 = np.zeros( (Nspots,) )
resi = np.zeros( (Nspots,) )

for i in range(Nspots):
    # phases are given by the minimum index, but
    # need correction if second coeff is <0
    if cm[mm[i],1,i] < 0:
        rp[i] -= 90
        # fix that coeff for use below
        cm[mm[i],1,i] *= -1
        
    # I_0 is given by the first coefficient
    I_0[i] = cm[mm[i],0,i]

    # now modulation is given by the ratio of c2/c1
    M_0[i] = cm[mm[i],1,i]/cm[mm[i],0,i]

    # also collect residuals of these minima
    resi[i] = rm[mm[i],i]


### Now, to return data, we assemble it into a single array, ###
### and push that to the root process. ###

results = np.vstack( (rp, I_0, M_0, resi) )
#print 'proc %d sending results.shape=',results.shape
comm.Send( results, dest=0, tag=myrank*5 )


# if myrank==3:
#     print myrank
#     print comm
#     import code
#     code.interact( local=locals() )


comm.Disconnect()



#return rp, I_0, M_0, resi
