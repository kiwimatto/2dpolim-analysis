import numpy as np


def wrapper_for_de( params, extras ):
    ex_angles = extras[0]
    em_angles = extras[1]
    Ftot      = extras[2]
    mod_depth_excitation = extras[3] 
    phase_excitation     = extras[4]
    mode                 = extras[5]
    use_least_sq         = extras[6]

    result = fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
                                                       mod_depth_excitation, phase_excitation, \
                                                       mode, use_least_sq )
    return result


def fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot, \
                                              mod_depth_excitation, phase_excitation, \
                                              mode, use_least_sq=True ):
    md_fu = params[0]
    th_fu = params[1]
    gr    = params[2]
    if not use_least_sq:
        et    = params[3]

    md_ex = mod_depth_excitation
    ph_ex = phase_excitation

    N_ex_angles = ex_angles.shape[1]   # because ex_angles was generated via meshgrid
    N_em_angles = ex_angles.shape[0]   # ...

    EX, EM = ex_angles.flatten(), em_angles.flatten()

    # calculate angle between off-axis dipoles in symmetric model
    if np.isnan(gr): gr=1.0
    alpha = 0.5 * np.arccos( .5*(((gr+2)*md_ex)-gr) )
    if np.isnan(alpha):
        raise ValueError( "alpha is nan. gr=%2.20e, md_ex=%2.20e" % (gr, md_ex) )

    ph_ii_minus = ph_ex -alpha
    ph_ii_plus  = ph_ex +alpha

    Fnoet  =    np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
    Fnoet += gr*np.cos( EX-ph_ex )**2 * np.cos( EM-ph_ex )**2
    Fnoet +=    np.cos( EX-ph_ii_plus )**2 * np.cos( EM-ph_ii_plus )**2
    Fnoet /= (2.0+gr)

    Fet   = .25 * (1+md_ex*np.cos(2*(EX-ph_ex))) * (1+md_fu*np.cos(2*(EM-th_fu-ph_ex)))

    if use_least_sq:
        # arrange Fet and Fnoet in two columns
        cols = np.hstack([Fet.reshape((Fet.size,1)), Fnoet.reshape((Fnoet.size,1))])

        # Ftot goes into a column...
        Ftot = Ftot.reshape((Ftot.size,1))
        # aaand we solve for coefficients of the two columns Fet and Fnoet,
        # such that Fot \approx a*Fet + b*Fnoet
        res = np.linalg.lstsq( cols, Ftot )
        # why two separate coefficients a and b, isn't it just a=et and b=(1-et)? 
        # Almost.
        # When we computed Fem, with et as part of the fit, we could scale Fem
        # by it's maximum, just as we do with Ftot, before computing the error.
        # Now, however, we do not know how to scale the model. So the coefficients
        # then include a scaling factor: a=A*et and b=A*(1-et).
        # Now we can recover both A=a+b and et=a/(a+b) 
    
        A  = np.sum(res[0])
        et = res[0][0]/A
        resi = res[1]

    else:
        Ftot /= np.max(Ftot)
        Fem   = et*Fet + (1-et)*Fnoet 
        Fem  /= np.max(Fem)
        A     = []
        resi  = np.sum( (Ftot.flatten()-Fem.flatten())**2 )


    if mode=="fitting":
        return resi   
    elif mode=="show_et_and_A":
        return et, A
    elif mode=="display":
        Fem  = A*et*Fet + A*(1-et)*Fnoet 
        Fem  = Fem.reshape((N_em_angles,N_ex_angles))        
        Ftot = Ftot.reshape((N_em_angles,N_ex_angles))
        import matplotlib.pyplot as plt
        fig = plt.figure( figsize=(12,6))
        axll = fig.add_axes([.05,.1,.1,.6])
        axlc = fig.add_axes([.15,.1,.3,.6])
        axlt = fig.add_axes([.15,.7,.3,.1])
        axlc.matshow( Ftot, interpolation='none', origin='bottom', vmin=0, vmax=np.max(Ftot) )
        axll.plot( np.mean( Ftot, axis=1 ), range(N_em_angles) )
        axlt.plot( np.mean( Ftot, axis=0 ) )
        axrl = fig.add_axes([.55,.1,.1,.6])
        axrc = fig.add_axes([.65,.1,.3,.6])
        axrt = fig.add_axes([.65,.7,.3,.1])
        axrc.matshow( np.abs(Ftot-Fem), interpolation='none', origin='bottom', vmin=0, vmax=np.max(Ftot) )
        axrl.plot( np.mean( Fem.reshape((N_em_angles,N_ex_angles)), axis=1 ), range(N_em_angles) )
        axrt.plot( np.mean( Fem.reshape((N_em_angles,N_ex_angles)), axis=0 ) )
        # ax = fig.add_subplot(1,2,2)
        # ax.imshow( Fem.reshape((N_em_angles,N_ex_angles)), interpolation='none' )
        plt.draw()
        return et, Fem
    else:
        raise ValueError("Don't know mode %s. Bombing out..." % (mode))
        return



def CosineFitter_new( angles, data, Nphases=91 ):

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

#    Nphases = 91
    phases  = np.linspace( 0, np.pi/2, Nphases )

    rm = np.zeros( (phases.size, Nspots) )    # residual matrix
    cm = np.zeros( (phases.size, 2, Nspots) ) # coefficient matrix

    # we'll init all phase-shifts at once, giving
    # a separate column for each phase offset
    # (these next two lines are beautiful, enjoy)
    shiftcos = np.cos( 2*reduce( np.add, np.meshgrid( -phases, angles ) ) )
    er=np.vstack((np.ones(shiftcos.shape),shiftcos)).reshape((shiftcos.shape[0],shiftcos.shape[1]*2), order='F')

    # now er starts with a columns of ones, followed by a phase-shifted cosine,
    # followed by another column of ones, followed by the next phase-shifted cosine, etc ...
    # Why this arrangement? --- first, we have the cosines still available when the fits are done,
    # and second, I can pass contiguous arrays to np.linalg.lstsq(), which should help speed-wise
    # (though i haven't tested this).

    for pi,phase in enumerate( phases ):
        # perform linear fit
        f = np.linalg.lstsq( er[:,2*pi:2*pi+2], data )

        # We will mess with the residuals now, because there are some nonsensical solutions
        # which we would like to discard before we pick (from the remaining solutions) the one 
        # with the best residual.
        # If the first coefficient is less than zero, then we have a problem: the fit becomes
        # negative (when the cosine hits its low). 
        firstcoeffislessthanzero = (f[0][0,:] < 0)
        # Additionally, if the absolute of the second coefficient is larger than the first 
        # coefficient, then, again, the fit becomes negative somewhere.
        badsecondcoeff = (np.abs(f[0][1,:]) > f[0][0,:])
        # Then all these fits are nonsensical to us: 
        nonsensical = np.logical_or( firstcoeffislessthanzero, badsecondcoeff )

        # We simply discard these (and hope that the next best fit is still decent) by setting 
        # their residual to infinity.
        f[1][ nonsensical ] = np.inf

        rm[pi,:] = f[1]
        cm[pi,:,:] = f[0][:]

#        print f[1]
        
    mm = np.argmin( rm, axis=0 )         # indices of minima (for each spot) w.r.t. tested phase
    rp = phases[mm]                      # resulting phases

    I_0 = np.zeros( (Nspots,) )
    M_0 = np.zeros( (Nspots,) )
    resi = np.zeros( (Nspots,) )
    fit  = np.zeros( (angles.size, Nspots) )

    rawfitpars = np.zeros( (2,Nspots) )

    for i in range(Nspots):

        rawfitpars[:,i] = cm[mm[i],:,i]

        # phases are given by the minimum index, but
        # need correction if second coeff is <0
        if cm[mm[i],1,i] < 0:
            rp[i] -= np.pi/2
            # fix that coeff for use below
            cm[mm[i],1,i] *= -1
        
        # I_0 is given by the first coefficient
        I_0[i] = cm[mm[i],0,i]

        # now modulation is given by the ratio of c2/c1
        M_0[i] = cm[mm[i],1,i]/cm[mm[i],0,i]

        # also collect residuals of these minima
        resi[i] = rm[mm[i],i]

        # take a look at the fit directly
        fit[:,i]     = np.dot( er[:,2*mm[i]:2*mm[i]+2], cm[mm[i],:,i] )
        try:
            assert np.all(fit[:,i]>=0)   # is it non-negative everywhere?
        except AssertionError:
            # print "############### ---------------- HOUSTON, WE HAVE A PROBLEM ---------------- ###############"
            # print ' residual   : ', rm[mm[i],i] 
            # print ' parameters : ', cm[mm[i],:,i]
            # print ' Setting everything to zero here. '
            # print "!!!!!!!!!!!!!!!!!!!!!!!"
            I_0[i] = 0
            M_0[i] = 0
#            raise AssertionError("The fit should be positive everywhere, but isn't. Something wrong?")


    return rp, I_0, M_0, resi, fit, rawfitpars, mm


def CosineFitter( angles, data ):

    assert angles.ndim == 1
    assert data.shape[0] == angles.size

    # if data is a 1d-array, then turn it into a 2d with singleton dimension
    if data.ndim==1:
        data = data.reshape( (data.size,1) )
    
    # how many data columns do we have?
    Nspots = data.shape[1]

    Nphases = 91
    phases  = np.linspace( 0, np.pi/2, Nphases )

    rm = residualmatrix    = np.zeros( (phases.size, Nspots) )
    cm = coefficientmatrix = np.zeros( (phases.size, 2, Nspots) )

    # init matrix for linear fit
    er = np.ones( (angles.size,2) )

    for pi,phase in enumerate( phases ):
        # write phase-shifted cos**2 into second column
        er[:,1] = np.cos( 2*(angles-phase) )
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
            rp[i] -= np.pi/2
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
    angles = np.linspace(0, np.pi, 181)
    data   = I*(1+M* np.cos(2*(angles-phase)) )
    if sigma>0:
        data   += np.random.normal( size=(angles.size,), scale=sigma )
    return angles, data
