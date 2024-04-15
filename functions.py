import numpy as np

def strain_1D_kin(ε, E, ypar):
    
    # parameters of the isotropic hardeninc curve
    Y0 = ypar[0]
    H = ypar[1]
    
    # function of the isotropic hardening curve
    def YF(x):
        return Y0+H*x
        
    # number of increments in the loading history
    n = len(ε)
    
    # initializing history variables
    Y = np.array([YF(0)])
    σ = np.array([0])
    λ = np.array([0])
    σ_eq = np.array([0])
    ε_p = np.array([0])
    α = np.array([0])

    for i in range(1,len(ε)):
        
        # assign variables to quantities at the beginning of the increment
        σ_n = σ[i-1]
        ε_n = ε[i-1]
        λ_n = λ[i-1]
        Y_n = Y[i-1]
        ε_pn = ε_p[i-1]
        α_n = α[i-1]

        # strain increment
        dε = ε[i] - ε_n
        
        # trial solution
        dλ = 0
        σ_t = σ_n + E*dε
        σ_1 = σ_t
        Y_1 = Y_n
        λ_1 = λ_n
        ε_p1 = ε_pn
        α_1 = α_n
        Ft = np.abs(σ_t - α_n) - Y0
        
        # stress update for plastic loading
        if Ft > 0:
            dλ = Ft/(E+H)
            λ_1 = λ_n + dλ
            dε_p = dλ*np.sign(σ_t - α_n)
            α_1 = α_1 + H*dε_p
            ε_p1 = ε_pn + dε_p
            Y_1 = YF(λ_1)
            σ_1 = σ_1 - E*dε_p

        # updating the history variables
        σ = np.append(σ, σ_1)
        Y = np.append(Y, Y_1)
        λ = np.append(λ, λ_1)
        σ_eq = np.append(σ_eq, np.abs(σ))
        ε_p = np.append(ε_p, ε_p1)
        α = np.append(α, α_1)
        
    return σ, λ, Y, σ_eq, ε_p, α

def strain_1D_iso(ε, E, ypar):
    '''Returns σ, λ, Y, σ_eq, ε_p, α for inputs:
        - strain: ε
        - Young's modulus: E
        - ypar: [Y0, H] for isotropic hardening model'''
    
    # parameters of the isotropic hardeninc curve
    Y0 = ypar[0]
    H = ypar[1]
    
    # function of the isotropic hardening curve
    def YF(x):
        return Y0+H*x
        
    # number of increments in the loading history
    n = len(ε)
    
    # initializing history variables
    Y = np.array([YF(0)])
    σ = np.array([0])
    λ = np.array([0])
    σ_eq = np.array([0])
    ε_p = np.array([0])

    for i in range(1,len(ε)):
        
        # assign variables to quantities at the beginning of the increment
        σ_n = σ[i-1]
        ε_n = ε[i-1]
        λ_n = λ[i-1]
        Y_n = Y[i-1]
        ε_pn = ε_p[i-1]

        # strain increment
        dε = ε[i] - ε_n
        
        # trial solution
        dλ = 0
        σ_t = σ_n + E*dε
        σ_1 = σ_t
        Y_1 = Y_n
        λ_1 = λ_n
        ε_p1 = ε_pn
        Ft = np.abs(σ_t) - Y_n
        
        # stress update for plastic loading
        if Ft > 0:
            dλ = Ft/(E+H)
            λ_1 = λ_n + dλ
            dε_p = dλ*np.sign(σ_t)
            ε_p1 = ε_pn + dε_p
            Y_1 = YF(λ_1)
            σ_1 = σ_1 - E*dε_p

        # updating the history variables
        σ = np.append(σ, σ_1)
        Y = np.append(Y, Y_1)
        λ = np.append(λ, λ_1)
        σ_eq = np.append(σ_eq, np.abs(σ))
        ε_p = np.append(ε_p, ε_p1)
        
    return σ, λ, Y, σ_eq, ε_p

def disp_to_eps(u, L0):
    '''Returns the eng. strain for u displacement and l initial length.
       Displacement is positive if it stretches l0.'''
    return u / L0

def eps_to_disp(eps, L0):
    '''Returns the displacement for eng. strain and L0 initial length.
       Strain is positive if it stretched l0.'''
    return eps * L0