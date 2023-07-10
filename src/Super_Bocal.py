def supersonic_inlet(M0,gamma):
    """
            Description: This method calculates the on design parameters of a supersonic inlet.

            Arguments:
                M0: Mach number
                gamma: Ratio of specific heats


            Returns: A dictionary containing the list of calculated outputs for each batch.
                A0_At: capture/throat area ratio
                Pty_PTx: pressure for the condition before the shock/pressure for Temperature after shock
                Ac_At: area at the throat (Ac) to the cross-sectional area at the nozzle exit (At) 
                Pty_Pt0: ratio of total pressure at the exit after the shock(Pty) to total pressure at the inlet (Pt0)

            """



    # Eq 3.11 para condição de M0

    ## A0_At = ((A / Astar) * (Pty / P0)) não precisa, fiz uma função

    # Função feita a partir do gráfico 10.23 do livro Elements of Propulsion
    def fun_A0_At(M):
        A0_At = -0.02720395692489004 * M **2 + 0.30287933903754993 * M + 0.7201226260439323
        return A0_At

    A0_At = fun_A0_At(M0) # Resultado da eq.11

    # Eq 3.12 para condição de M0

    # # Função feita a partir do gráfico 10.17 do livro Elements of Propulsion
    def fun_Pty_Ptx(M):
        Pty_Ptx = 0.937222346061776 + 0.41669957451880146 * M - 0.37653137074222387 * M ** 2 + 0.055265789273162676 * M ** 3
        return Pty_Ptx

    Pty_Ptx = fun_Pty_Ptx(M0)

    # Eq 3.13
    A_Astar = (1 / M0) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M0 ** 2)) ** ((gamma + 1) / (2 * (gamma - 1)))

    # Eq 3.12 para condição de M0
    Ac_At = A0_At * ((A_Astar) * (Pty_Ptx))

    # Eq 3.14
    Pty_Pt0 = A0_At*(A_Astar)**-1
    Astar0_Astary = Pty_Pt0
    phoVy_phoV0 = Pty_Pt0
    phoy_pho0 = Pty_Pt0

    # Eq 3.15 são três igualdades 
    Astarx_Astary = Pty_Ptx
    phox_phoy = Pty_Ptx

    return(A0_At,A_Astar,Pty_Ptx,Ac_At,Pty_Pt0)



