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



def fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
                                              mod_depth_excitation, phase_excitation, mode ):
#def fit_portrait_single_funnel_symmetric( params, extras ):
    # ex_angles = extras[0]
    # em_angles = extras[1]
    # Ftot      = extras[2]
    # mod_depth_excitation = extras[3] 
    # phase_excitation     = extras[4]
    # mode                 = extras[5]

    md_fu = params[0]
    et    = params[1]
    th_fu = params[2]     # in radians
    gr    = params[3]

    md_ex = mod_depth_excitation
    ph_ex = phase_excitation * np.pi/180.0

    # if data_style=='vector':
    #     EX, EM = ex_angles, em_angles 
    # elif data_style=='matrix':    
    #     EX, EM = np.meshgrid( ex_angles, em_angles )

    N_ex_angles = ex_angles.shape[1]   # because ex_angles was generated via meshgrid
    N_em_angles = ex_angles.shape[0]   # ...

    EX, EM = ex_angles.flatten()*np.pi/180.0, em_angles.flatten()*np.pi/180.0

    # calculate angle between off-axis dipoles in symmetric model
    alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )

    ph_ii_minus = ph_ex -alpha
    ph_ii_plus  = ph_ex +alpha

    EnNoET  = np.zeros_like( EX )
    EnNoET +=    np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
    EnNoET += gr*np.cos( EX-ph_ex )**2 * np.cos( EM-ph_ex )**2
    EnNoET +=    np.cos( EX-ph_ii_plus )**2 * np.cos( EM-ph_ii_plus )**2

    Fnoet  = EnNoET / (2+gr)
#    Fnoet /= np.sum(Fnoet)

    Fet   = .25 * (1+md_ex*np.cos(2*(EX-ph_ex))) * (1+md_fu*np.cos(2*(EM-th_fu-ph_ex)))
#    Fet  /= np.sum(Fet)


    # We're trying to work out what et should be, along the lines of
    #    Ftot       = et*Fet + (1-et)*Fnoet 
    # We can rewrite this as follows:
    #    Ftot       = et*Fet + Fnoet - et*Fnoet
    #    Ftot-Fnoet = et*(Fet-Fnoet)
    # This is obviously a linear problem, ie can be thought of as Ax=b,
    # where A=Fet-Fnoet is the 'coefficient matrix' and b=Ftot-Fnoet
    # is the 'inhomogeneity'. We use the linear algebra least-squares
    # solver to do this step.

#    A = (Fet-Fnoet).reshape( (Fet.size,1) )
#    print mode

    # if mode=='generate data':
    #     et = Ftot
    # else:
    #     Ftot = Ftot.flatten()
    #     res = np.linalg.lstsq( A, Ftot-Fnoet )
    #     et = res[0]
    #     resi = res[1]
    #     print et

    Fem = et*Fet + (1-et)*Fnoet 

    # max-normalize both Ftot and Fem
    Ftot /= np.max(Ftot)
    Fem  /= np.max(Fem)

    resi = np.sum( (Ftot.flatten()-Fem)**2 )

    if mode=="fitting":
        return resi   # this is the residual straight from the least-squares fit
#        return Ftot - (et*Fet + (1-et)*Fnoet)
    elif mode=="display":
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(2,2,1)
        ax.imshow( ex_angles.reshape((N_em_angles,N_ex_angles)), interpolation='none' )
        ax = fig.add_subplot(2,2,2)
        plt.imshow( em_angles.reshape((N_em_angles,N_ex_angles)), interpolation='none' )
        ax = fig.add_subplot(2,2,3)
        plt.imshow( Ftot.reshape((N_em_angles,N_ex_angles)), interpolation='none' )
        ax = fig.add_subplot(2,2,4)
        ax.imshow( Fem.reshape((N_em_angles,N_ex_angles)), interpolation='none' )
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




def CosineFitter_new( angles, data ):

    assert angles.ndim == 1
    assert data.shape[0] == angles.size

    # import matplotlib.pyplot as plt
    # plt.plot(angles, data[:,0],'o:')
    # print data.shape

    # if data is a 1d-array, then turn it into a 2d with singleton dimension
    if data.ndim==1:
        data = data.reshape( (data.size,1) )
    
    # how many data columns do we have?
    Nspots = data.shape[1]

    Nphases = 91
    phases  = np.linspace( 0, 90, Nphases )

    rm = residualmatrix    = np.zeros( (phases.size, Nspots) )
    cm = coefficientmatrix = np.zeros( (phases.size, 2, Nspots) )

    # we'll init all phase-shifts at once, giving
    # a separate column for each phase offset
    # (these next two lines are beautiful, enjoy)
    shiftcos = np.cos( 2*np.pi/180.0*reduce( np.add, np.meshgrid( -phases, angles ) ) )
    er=np.vstack((np.ones(shiftcos.shape),shiftcos)).reshape((shiftcos.shape[0],shiftcos.shape[1]*2), order='F')

    # now er starts with a columns of ones, followed by a phase-shifted cosine,
    # followed by another column of ones, followed by the next phase-shifted cosine, etc ...
    # Why this arrangement? --- first, we have the cosines still available when the fits are done,
    # and second, I can pass contiguous arrays to np.linalg.lstsq(), which should help speed-wise
    # (though i haven't tested this).

    # bestpis   = np.zeros( (Nspots,), dtype=np.int )
    # bestresis = np.ones( (Nspots,) )*np.inf
    # bestfitpa = np.zeros( (2,Nspots) )
    for pi,phase in enumerate( phases ):
        # perform linear fit
        f = np.linalg.lstsq( er[:,2*pi:2*pi+2], data )
        residualmatrix[pi,:] = f[1]
        cm[pi,:,:] = f[0][:]

        # # for which spots is the current fit better than anythings we've seen before?
        # wherebetter = f[1]<bestresis
        # # update the list of best known residuals
        # bestresis[wherebetter] = f[1][wherebetter]
        # # write the phase index into a similar list
        # bestpis[wherebetter]   = pi
        # # as well as the resulting fit parameters
        # bestfitpa[:,wherebetter] = f[0][:,wherebetter]

    mm = minindices = np.argmin( residualmatrix, axis=0 )
    rp = resultingphases = phases[minindices]

    I_0 = np.zeros( (Nspots,) )
    M_0 = np.zeros( (Nspots,) )
    resi = np.zeros( (Nspots,) )
    fit  = np.zeros( (angles.size, Nspots) )

    rawfitpars = np.zeros( (2,Nspots) )

    for i in range(Nspots):

        fit[:,i]     = np.dot( er[:,2*mm[i]:2*mm[i]+2], cm[mm[i],:,i] )
        rawfitpars[:,i] = cm[mm[i],:,i]

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

    return rp, I_0, M_0, resi, fit, rawfitpars


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
    fit  = np.zeros( (angles.size, Nspots) )


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

    return rp, I_0, M_0, resi, [], []




def generate_fake_data( phase, I, M, sigma=0 ):
    angles = np.linspace(0, 180, 181)
    data   = I*(1+M* np.cos(2*(angles-phase)*np.pi/180.0) )
    if sigma>0:
        data   += np.random.normal( size=(angles.size,), scale=sigma )
    return angles, data
