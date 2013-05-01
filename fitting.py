import numpy as np


# def err_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot,\
#                                               mod_depth_excitation, phase_excitation,\
#                                               which_error='chi2' ):

#     Fem = fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
#                                                     mod_depth_excitation, phase_excitation, \
#                                                     data_style )
#     Fem  /= np.sum(Fem)    # unless an axis is specified, np.sum() will sum all elements, ...
#     Ftot /= np.sum(Ftot)   # ... regardless of array dimension

#     # compute error
#     if which_error=='chi2':
#         err = np.sum( (Fem-Ftot)**2 )
#     elif which_error=='R2':
#         err = np.sum( (Fem-Ftot)**2 )/np.sum( (Ftot-np.mean(Ftot))**2 )
#     else:
#         ValueError("Unknown value for which_error: %s  (should be 'chi2' or 'R2')" % (which_error))



# def fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
#                                               mod_depth_excitation, phase_excitation, mode ):
def fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
                                              mod_depth_excitation, phase_excitation, mode ):
    mod_depth_funnel = mf = params[0]
#    energy_transfer  = et = params[1]
    theta_funnel     = tf = params[1]
    geometric_ratio  = gr = params[2]

    m_ex  = mod_depth_excitation
    ph_ex = phase_excitation

    # if data_style=='vector':
    #     EX, EM = ex_angles, em_angles 
    # elif data_style=='matrix':    
    #     EX, EM = np.meshgrid( ex_angles, em_angles )

    EX, EM, Ftot = ex_angles.flatten()*np.pi/180.0, em_angles.flatten()*np.pi/180.0, Ftot.flatten()

    # calculate angle between off-axis dipoles in symmetric model
    alpha = 0.5 * np.arccos( .5*(((gr+2)*m_ex)-gr) )

    ph_ii_minus = ph_ex -alpha
    ph_ii_plus  = ph_ex +alpha

    EnNoET  = np.zeros_like( EX )
    EnNoET +=    np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
    EnNoET += gr*np.cos( EX-ph_ex )**2 * np.cos( EM-ph_ex )**2
    EnNoET +=    np.cos( EX-ph_ii_plus )**2 * np.cos( EM-ph_ii_plus )**2

    Fnoet  = EnNoET / (2+gr)
    Fnoet /= np.sum(Fnoet)

    Fet   = .25 * (1+m_ex*np.cos(2*(EX-ph_ex))) * (1+mf*np.cos(2*(EM-theta_funnel-ph_ex)))
    Fet  /= np.sum(Fet)


    # We're trying to work out what et should be, along the lines of
    #    Ftot       = et*Fet + (1-et)*Fnoet 
    # We can rewrite this as follows:
    #    Ftot       = et*Fet + Fnoet - et*Fnoet
    #    Ftot-Fnoet = et*(Fet-Fnoet)
    # This is obviously a linear problem, ie can be thought of as Ax=b,
    # where A=Fet-Fnoet is the 'coefficient matrix' and b=Ftot-Fnoet
    # is the 'inhomogeneity'. We use the linear algebra least-squares
    # solver to do this step.

    A = (Fet-Fnoet).reshape( (Fet.size,1) )
    res = np.linalg.lstsq( A, Ftot-Fnoet )
    et = res[0]

    print et

    Fem = et*Fet + (1-et)*Fnoet 

    if mode=="fitting":
        return res[1]   # this is the residual straight from the least-squares fit
#        return Ftot - (et*Fet + (1-et)*Fnoet)
    elif mode=="display":
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(2,2,1)
        ax.imshow( ex_angles.reshape((181,181)), interpolation='none' )
        ax = fig.add_subplot(2,2,2)
        plt.imshow( em_angles.reshape((181,181)), interpolation='none' )
        ax = fig.add_subplot(2,2,3)
        plt.imshow( Ftot.reshape((181,181)), interpolation='none' )
        ax = fig.add_subplot(2,2,4)
        ax.imshow( Fem.reshape((181,181)), interpolation='none' )

        return et, Fem
    else:
        raise ValueError("Don't know mode %s. Bombing out..." % (mode))
        return





def CosineFitter_mpi_master( angles, data ):

    assert angles.ndim == 1
    assert data.shape[0] == angles.size

    # if data is a 1d-array, then turn it into a 2d with singleton dimension
    if data.ndim==1:
        data = data.reshape( (data.size,1) )
    
    import sys
    from mpi4py import MPI
    # comm = MPI.COMM_WORLD
    # myrank = comm.Get_rank()
    # nprocs = comm.Get_size()

    Nprocs = 4

#    sdata = np.hsplit(data, Nprocs)

    comm = MPI.COMM_SELF
    comm.Set_name('spawner')
    slaveintercomm = comm.Spawn(sys.executable,
                                args=['cosine_fitter_mpi_slave.py', \
                                          str(data.shape[0]), \
                                          str(data.shape[1]) ],
                                maxprocs=Nprocs)

#    print 'spawner'
#    print comm
    
    for n in range(Nprocs):
        localcolumns = range(n,data.shape[1],Nprocs)
#        print 'spawner sending data, len(localcolumns)=',len(localcolumns)
        # send out data
        slaveintercomm.Send( angles, dest=n, tag=n*11 )
        d = data[:,localcolumns].copy()
#        print d
        slaveintercomm.Send( d, dest=n, tag=n*123 )

    results = np.zeros( (4, data.shape[1]) )
    for n in range(Nprocs):        
        localcolumns = range(n,data.shape[1],Nprocs)
#        print "receiving...%d" % (n)
        r = np.zeros( (4,len(localcolumns)) )
#        print 'spawner expecting r.shape=',r.shape
        slaveintercomm.Recv( r, source=n, tag=n*5 )
        results[:,localcolumns] = r   

    slaveintercomm.Disconnect()

    return results[0,:], results[1,:], results[2,:], results[3,:]




def CosineFitter( angles, data ):

    assert angles.ndim == 1
    assert data.shape[0] == angles.size

    # if data is a 1d-array, then turn it into a 2d with singleton dimension
    if data.ndim==1:
        data = data.reshape( (data.size,1) )
    
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

    return rp, I_0, M_0, resi



def generate_fake_data( phase, I, M, sigma=0 ):
    angles = np.linspace(0, 180, 181)
    data   = I*(1+M* np.cos(2*(angles-phase)*np.pi/180.0) )
    if sigma>0:
        data   += np.random.normal( size=(angles.size,), scale=sigma )
    return angles, data
