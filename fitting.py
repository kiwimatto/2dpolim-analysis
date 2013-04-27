import numpy as np


def err_portrait_single_funnel_symmetric( params, ex_angles, em_angles, Ftot,\
                                              mod_depth_excitation, phase_excitation,\
                                              which_error='chi2' ):

    Fem = fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, \
                                                    mod_depth_excitation, phase_excitation, \
                                                    data_style)
    Fem  /= np.sum(Fem)    # unless an axis is specified, np.sum() will sum all elements, ...
    Ftot /= np.sum(Ftot)   # ... regardless of array dimension

    # compute error
    if which_error=='chi2':
        err = np.sum( (Fem-Ftot)**2 )
    elif which_error=='R2':
        err = np.sum( (Fem-Ftot)**2 )/np.sum( (Ftot-np.mean(Ftot))**2 )
    else:
        ValueError("Unknown value for which_error: %s  (should be 'chi2' or 'R2')" % (which_error))



def fit_portrait_single_funnel_symmetric( params, ex_angles, em_angles, \
                                              mod_depth_excitation, phase_excitation, \
                                              data_style):
    mod_depth_funnel = mf = params[0]
    energy_transfer  = et = params[1]
    theta_funnel     = tf = params[2]
    geometric_ratio  = gr = params[3]

    m_ex  = mod_depth_excitation
    ph_ex = phase_excitation

    # if data_style=='vector':
    #     EX, EM = ex_angles, em_angles 
    # elif data_style=='matrix':    
    #     EX, EM = np.meshgrid( ex_angles, em_angles )
    EX, EM = ex_angles, em_angles 

    # calculate angle between off-axis dipoles in symmetric model
    alpha = 0.5 * np.arccos( .5*(((gr+2)*m_ex)-gr) )

    ph_ii_minus = ph_ex -alpha
    ph_ii_plus  = ph_ex +alpha

    EnNoET  = np.zeros_like( EX )
    EnNoET +=    np.cos( EX-ph_ii_minus )**2 * np.cos( EM-ph_ii_minus )**2
    EnNoET += gr*np.cos( EX-ph_ex )**2 * np.cos( EM-ph_ex )**2
    EnNoET +=    np.cos( EX-ph_ii_plus )**2 * np.cos( EM-ph_ii_plus )**2

    Fnoet = EnNoET / (2+gr)

    Fet   = .25 * (1+m_ex*np.cos(2*(EX-ph_ex))) * (1+mf*np.cos(2*(EM-theta_funnel-ph_ex)))

    Fem   = et*Fet + (1-et)*Fnoet

    return Fem.flatten()



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
