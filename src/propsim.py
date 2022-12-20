from ambiance import Atmosphere
import math


class AircraftEngines:

    def __init__(self, height):
        self.atm = Atmosphere(height)
        self.T0 = self.atm.temperature
        self.P0 = self.atm.pressure
        self.a0 = self.atm.speed_of_sound

    def ideal_turbojet(self, M0, gamma, cp, hpr, Tt4, pi_c):
        """
        Description: This method calculates the on design parameters of an ideal turbojet engine.

        Arguments:
            M0: Mach number                             [  -  ]
            gamma: Ratio of specific heats              [kJ/kgK]
            cp: Specific heat at constant pressure      [ J/K ]
            hpr: Low heating value of fuel              [kJ/kg]
            Tt4: Total temperature leaving the burner   [  K  ]
            pi_c: Compressor total pressure ratio       [  -  ]

        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
        """    

        R = (gamma - 1)/gamma * cp
        V0 = self.a0 * M0 

        tau_r = 1 + (gamma - 1)/2 * M0**2

        tau_c = pi_c**((gamma - 1)/gamma)

        tau_lambda = Tt4/self.T0
        f = cp * self.T0/hpr * (tau_lambda - tau_r * tau_c)

        tau_t = 1 - tau_r/tau_lambda * (tau_c - 1)
        V9_a0 = (2/(gamma - 1 )* tau_lambda/(tau_r * tau_c) * (tau_r * tau_c * tau_t - 1))**(1/2)

        F_m0 = self.a0 * (V9_a0 - M0)
        S = f/F_m0 
        eta_T = 1 - 1/(tau_r * tau_c)
        eta_P = 2 * M0/(V9_a0 + M0)
        eta_Total = eta_P * eta_T

        return (F_m0, f, S, eta_T, eta_P, eta_Total)

    def real_turbojet(self, M0, gamma_c, gamma_t, cp_c, cp_t, hpr, Tt4 , pi_c, pi_d_max, pi_b, pi_n,e_c,e_t,eta_b,eta_m,P0_P9):
        """
        Description: This method calculates the on design parameters of an non ideal turbojet engine.

        Arguments:
            M0: Mach number
            gamma_c: Ratio of specific heats in the compressor          [kJ/kgK]
            gamma_t: Ratio of specific heats in the turbine             [kJ/kgK]
            cp_c: Specific heat at constant pressure in the compressor  [ J/K ]
            cp_t: Specific heat at constant pressure in the turbine     [ J/K ]
            hpr: Low heating value of fuel                              [kJ/kg]
            Tt4: Total temperature leaving the burner                   [  K  ]
            pi_c: Compressor total pressure ratio                       [  -  ]
            pi_d_max: Difuser maximum total pressure ratio              [  -  ]
            pi_b: Burner total pressure ratio                           [  -  ]
            pi_n: Nozzle total pressure ratio                           [  -  ]
            e_t: Polytropic compressor efficiency                       [  -  ]
            e_t: Polytropic turbine efficiency                          [  -  ]
            eta_b: Combustor efficiency                                 [  -  ]
            eta_m: Mechanical efficiency                                [  -  ]
            P0_P9: Ratio between Po/P9                                  [  -  ]


        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
        """  

        # O primeiro passo para a solução numerica do motor Turbojet não ideal é o cálculo das constantes gerais dos gases proposta por Clapeyron, para o compressor (Rc) e para a turbina (Rt) em J/(kg.K) a partir dos dados de entrada.
        R_c    = ( (gamma_c - 1)/gamma_c )*(cp_c)
        R_t    = ( (gamma_t - 1)/gamma_t )*(cp_t)

        # Local sound speed (a0) for flight conditions and aircraft speed (V0)
        a0    = math.sqrt(gamma_c*R_c*self.T0)
        V0    = a0*M0

        # Free flow parameters
        tal_r    = 1 + ( (gamma_c - 1)/2 )*(M0**2)
        pi_r     = tal_r**(gamma_c/(gamma_c - 1))

        if (M0 <= 1):
            eta_r  = 1
        elif (M0 > 1 and M0 <= 5):
            eta_r  = 1 - 0.075*((M0 - 1)**(1.35))
        elif (M0 > 5):
            eta_r  = 800/(M0**4 + 935)

        # Diffuser
        pi_d   = pi_d_max*eta_r

        # Compressor
        tal_c    = (pi_c)**((gamma_c - 1)/(gamma_c*e_c))
        eta_c    = (pi_c**((gamma_c - 1)/gamma_c) - 1)/(tal_c - 1)

        # Burner
        tal_lambda    = cp_t*Tt4/(cp_c*self.T0)
        f             = (tal_lambda-tal_r*tal_c)/( (hpr*eta_b/(cp_c*self.T0)) - tal_lambda)

        # Turbine
        tal_t  = 1 - (1/(eta_m*(1+f)))*(tal_r/tal_lambda)*(tal_c-1) 
        pi_t   = tal_t**(gamma_t/((gamma_t-1)*e_t))
        eta_t  = (1 - tal_t)/(1 - tal_t**(1/e_t))

        # Parameters related to the output of the core after the turbine
        Pt9_P9 = (P0_P9)*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n
        M9     = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
        T9_T0  = (tal_lambda*tal_t)/( ((Pt9_P9)**((gamma_t-1)/gamma_t)))*(cp_c/cp_t)   
        V9_a0  = M9*( math.sqrt( gamma_t*R_t*(T9_T0)/(gamma_c*R_c) ) )

        # Results
        F_m0      = a0*((1+f)*V9_a0 - M0 + (1+f)*R_t*T9_T0/(R_c*V9_a0)*(1-P0_P9)/gamma_c)
        S         = f/F_m0
        eta_T   = ((a0**2)*( (1+f)*V9_a0**2 - M0**2 )/(2*f*(hpr)))
        eta_P  = ( 2*V0*F_m0/( (a0**2)*( (1+f)*V9_a0**2 - M0**2 ) ) )
        eta_Total = eta_T*eta_P

        return (F_m0, f, S, eta_T, eta_P, eta_Total)


    def ideal_turbofan(self, M0, gamma, cp, hpr, Tt4, pi_c, pi_f, alpha):
        """
        Description: This method calculates the on design parameters of an ideal turbofan engine.

        Arguments:
            M0: Mach number                             [  -  ]
            gamma: Ratio of specific heats              [kJ/kgK]
            cp: Specific heat at constant pressure      [ J/K ]
            hpr: Low heating value of fuel              [kJ/kg]
            Tt4: Total temperature leaving the burner   [  K  ]   
            pi_c: Compressor total pressure ratio       [  -  ]
            
        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
            FR: Thrust ratio                [6]
        """    

        R = (gamma - 1)/gamma * cp
        V0 = self.a0 * M0 

        tau_r = 1 + (gamma - 1)/2 * M0**2

        tau_f = pi_f**((gamma - 1)/gamma)
        V19_a0 = (2/(gamma - 1) * (tau_r * tau_f - 1))**(1/2)

        tau_c = pi_c**((gamma - 1)/gamma)
        tau_lambda = Tt4/self.T0

        f = cp * self.T0/hpr * (tau_lambda - tau_r * tau_c)

        tau_t = 1 - tau_r/tau_lambda * (tau_c - 1)
        V9_a0 = (2/(gamma - 1) * (tau_lambda - tau_r*(tau_c - 1 + alpha * (tau_f - 1)) - tau_lambda/(tau_r * tau_c)))**(1/2)

        F_m0 = self.a0 * 1/(1 + alpha) * (V9_a0 - M0 + alpha * (V19_a0 - M0))

        f = cp * self.T0/hpr * (tau_lambda - tau_r * tau_c)

        S = f/((1 + alpha) * F_m0)
        eta_T = 1 - 1/(tau_r * tau_c)
        eta_P = 2 * M0 * (V9_a0 - M0 + alpha * (V19_a0 - M0))/((V9_a0 * self.a0)**2/(self.a0**2) - M0**2 + alpha * ((V19_a0*self.a0)**2/(a0**2) - M0**2))
        eta_Total = eta_P * eta_T
        FR = (V9_a0 - M0)/(V19_a0 - M0)

        return (F_m0, f, S, eta_T, eta_P, eta_Total, FR)


    def real_turbofan(self, 
        M0, 
        gamma_c,
        gamma_t,
        cp_c,
        cp_t,
        hpr, 
        Tt4,
        pi_d_max,
        pi_b,
        pi_n,
        pi_fn,
        e_cL,
        e_cH,
        e_f,
        e_tL,
        e_tH,
        eta_b,
        eta_mL,
        eta_mH,
        P0_P9,
        P0_P19,
        tau_n,
        tau_fn,
        pi_cL,
        pi_cH,
        pi_f,
        alfa
        ):

        """
        Description: This method calculates the on design parameters of an ideal turbofan engine.

        Arguments:
            M0: Mach number                             [  -  ]
            gamma: Ratio of specific heats              [kJ/kgK]
            cp: Specific heat at constant pressure      [ J/K ]
            hpr: Low heating value of fuel              [kJ/kg]
            Tt4: Total temperature leaving the burner   [  K  ]
            pi_c: Compressor total pressure ratio       [  -  ]
            
        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
            FR: Thrust ratio                [6]
        """

        R_c = (gamma_c - 1)/gamma_c*cp_c #J/(kg.K)
        R_t = (gamma_t - 1)/gamma_t*cp_t #J/(kg.K)

        a0 = (gamma_c*R_c*T0)**(1/2) #m/s
        V0 = a0*M0 # m/s

        # Calculo dos parametros do escoamento livre
        tau_r = 1 + (gamma_c - 1)/2*M0**2

        pi_r = tau_r**(gamma_c/(gamma_c - 1))

        if M0 <= 1:
            eta_r = 1
        else:
            eta_r = 1 - 0.075*(M0 - 1)**1.35

        # Calculo de parametros do Difusor
        pi_d = pi_d_max*eta_r
        tau_d = pi_d**((gamma_c - 1)/gamma_c)

        # Calculo de parametros do Fan
        tau_f = pi_f**((gamma_c - 1)/(gamma_c*e_f))
        eta_f = (pi_f**((gamma_c - 1)/gamma_c) - 1)/(tau_f - 1)

        # Relação de entalpia entre o queimador e o escoamento livre
        tau_lambda = cp_t*Tt4/(cp_c*T0)

        # Calculo de parametros do compressor
        tau_cL = pi_cL**((gamma_c - 1)/(gamma_c*e_cL))
        tau_cH = pi_cH**((gamma_c - 1)/(gamma_c*e_cH))

        eta_cL = (pi_cL**((gamma_c - 1)/gamma_c) - 1)/(tau_cL - 1)
        eta_cH = (pi_cH**((gamma_c - 1)/gamma_c) - 1)/(tau_cH - 1)

        # Calculo de parametros na turbina
        f = (tau_lambda - tau_r*tau_d*tau_f*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*T0) - tau_lambda) #kgFuel/kgAir
        print("f = {}\n".format(f))

        tau_tH = 1 - (tau_cH - 1)/(1 + f)/tau_lambda*tau_r*tau_d*tau_f*tau_cL*eta_mH
        tau_tL = 1 - ((alfa*(tau_f - 1) + (tau_cL - 1))*eta_mL/(1 + f)*tau_r*tau_d/tau_lambda/tau_tH)

        pi_tH = tau_tH**(gamma_t/((gamma_t - 1)*e_tH))
        pi_tL = tau_tL**(gamma_t/((gamma_t - 1)*e_tL))

        eta_tH = (tau_tH - 1)/(pi_tH**((gamma_t - 1)/gamma_t)-1)
        eta_tL = (tau_tL - 1)/(pi_tL**((gamma_t - 1)/gamma_t)-1)

        # Parametros referentes a saida do core apos a turbina
        Pt9_P9 = P0_P9*pi_r*pi_d*pi_f*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n

        M9 = (2/(gamma_c - 1)*(Pt9_P9**((gamma_c - 1)/gamma_c) - 1))**(1/2)

        Tt9_T0 = cp_c/cp_t*tau_lambda*tau_tH*tau_tL*tau_n

        T9_T0 = Tt9_T0/Pt9_P9**((gamma_t - 1)/gamma_t)

        V9_a0 = M9*(R_t*gamma_t/(R_c*gamma_c)*T9_T0)**(1/2)

        # Parametros referentes a saida do bypass apos o fan
        Pt19_P19 = P0_P19*pi_r*pi_d*pi_f*pi_fn

        M19 = (2/(gamma_c - 1)*(Pt19_P19**((gamma_c - 1)/gamma_c) - 1))**(1/2)

        Tt19_T0 = tau_r*tau_d*tau_f*tau_fn

        T19_T0 = Tt19_T0/Pt19_P19**((gamma_c - 1)/gamma_c)

        V19_a0 = M19*(T19_T0)**(1/2)

        # Outputs da analise do motor
        # Calculo dos parametros de desempenho do motor
        FF_m0 = alfa/(1 + alfa)*a0*(V19_a0 - M0 + T19_T0/V19_a0*(1 - P0_P19)/gamma_c) #N/(kg/s)
        FC_m0 = 1/(1 + alfa)*a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gamma_c) #N/(kg/s)
        F_m0 = FF_m0  + FC_m0 #N/(kg/s)

        S = f/((1 + alfa)*F_m0) #(kgFuel/s)/N

        FR = FF_m0/FC_m0

        # Calculo dos parametros de eficiencia global do motor
        eta_T = a0*a0*((1 + f)*V9_a0*V9_a0 + alfa*(V19_a0*V19_a0)- (1 + alfa)*M0*M0)/(2*f*hpr)

        eta_P = 2*M0*((1 + f)*V9_a0 + alfa*V19_a0 - (1 + alfa)*M0)/((1 + f)*(V9_a0**2) + alfa*V19_a0**2 - (1 + alfa)*M0**2)

        eta_Total = eta_P*eta_T

        return (F_m0, f, S, eta_T, eta_P, eta_Total, FR)




    def real_turbofan_off_design(self, 
        # Escolhas
        M0, 
        gamma_c, 
        gamma_t, 
        cp_c, 
        cp_t, 
        hpr, 
        Tt4, 
        pi_d_max, 
        pi_b, 
        pi_c,
        pi_tH, 
        pi_n,
        pi_fn,
        tau_tH,
        eta_f,
        eta_cL,
        eta_cH,
        eta_b,
        eta_mL,
        eta_mH,
        eta_tL,

        # Referência
        M0_R,
        T0_R,
        P0_R,
        tau_r_R,
        tau_lambda_R,
        pi_r_R,
        Tt4_R,
        pi_d_R,
        pi_f_R,
        pi_cH_R,
        pi_cL_R,
        pi_tL_R,
        tau_f_R,
        tau_tL_R,
        alfa_R,
        M9_R,
        M19_R,
        m0_R

        ):

        tau_cH_R = pi_cH_R**((gamma_c - 1)/(gamma_c))
        tau_cL_R = pi_cL_R**((gamma_c - 1)/(gamma_c))

        # Equações
        R_c = (gamma_c - 1)/gamma_c*cp_c #J/(kg.K)
        R_t = (gamma_t - 1)/gamma_t*cp_t #J/(kg.K)
        a0 = (gamma_c*R_c*T0)**(1/2) #m/s
        V0 = a0*M0 #m/s
        tau_r = 1 + (gamma_c - 1)/2*M0**2
        pi_r = tau_r**(gamma_c/(gamma_c - 1))
        if M0 <= 1:
            eta_r = 1
        else:
            eta_r = 1 - 0.075*(M0 - 1)**1.35

        pi_d = pi_d_max*eta_r
        tau_lambda = cp_t*Tt4/(cp_c*T0)

        teste = 10
        while teste > 0.0001:
            tau_tL = tau_tL_R
            tau_f = tau_f_R
            tau_cL = tau_cL_R
            pi_tL = pi_tL_R
            pi_cL = pi_cL_R
            tau_cH = 1 + Tt4/T0/(Tt4_R/T0_R)*(tau_f_R*tau_cL_R)/(tau_r*tau_cL)*(tau_cH_R - 1)
            pi_cH = (1 + eta_cH*(tau_cH - 1))**(gamma_c/(gamma_c - 1))
            pi_f = (1 + (tau_f - 1)*eta_f)**(gamma_c/(gamma_c - 1))
            Pt19_P0 = pi_r*pi_d*pi_f*pi_fn
            if Pt19_P0 < ((gamma_c + 1)/2)**(gamma_c/(gamma_c - 1)):
                Pt19_P19 = Pt19_P0
            else:
                Pt19_P19 = ((gamma_c + 1)/2)**(gamma_c/(gamma_c - 1))
            M19 = (2/(gamma_c - 1)*(Pt19_P19**((gamma_c - 1)/gamma_c) - 1))**(1/2)
            Pt9_P0 = pi_r*pi_d*pi_f*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n
            if Pt9_P0 < ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1)):
                Pt9_P9 = Pt9_P0
            else:
                Pt9_P9 = ((gamma_t + 1)/2)**(gamma_t/(gamma_t - 1))

            M9 = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
            MFP_M19 = M19*(gamma_c/R_c)**(1/2)*(1 + (gamma_c - 1)/2*M19**2)**((gamma_c + 1)/(2*(1 - gamma_c)))
            MFP_M19_R = M19_R*(gamma_c/R_c)**(1/2)*(1 + (gamma_c - 1)/2*M19_R**2)**((gamma_c + 1)/(2*(1 - gamma_c)))
            alfa = alfa_R*pi_cL_R*pi_cH_R/pi_f_R/(pi_cL*pi_cH/pi_f)*((tau_lambda)/(tau_r*tau_f)/(tau_lambda_R/(tau_r_R*tau_f_R)))**(1/2)*MFP_M19/MFP_M19_R
            tau_f = 1 + (tau_f_R - 1)*((1 - tau_tL)/(1 - tau_tL_R)*(tau_lambda/tau_r)/(tau_lambda_R/tau_r_R)*(tau_cL_R - 1 + alfa_R*(tau_f_R - 1))/(tau_cL_R - 1 + alfa*(tau_f_R - 1)))
            tau_cL = 1 + (tau_f - 1)*(tau_cL_R - 1)/(tau_f_R - 1)
            pi_cL = (1 + eta_cL*(tau_cL - 1))**(gamma_c/(gamma_c - 1))
            tau_tL = 1 - eta_tL*(1 - pi_tL**((gamma_t - 1)/gamma_t))
            MFP_M9 = M9*(gamma_t/R_t)**(1/2)*(1 + (gamma_t - 1)/2*M9**2)**((gamma_t + 1)/(2*(1 - gamma_t)))
            MFP_M9_R = M9_R*(gamma_t/R_t)**(1/2)*(1 + (gamma_t - 1)/2*M9_R**2)**((gamma_t + 1)/(2*(1 - gamma_t)))
            pi_tL = pi_tL_R*(tau_tL/tau_tL_R)**(1/2)*MFP_M9_R/MFP_M9
            teste = abs(tau_tL - tau_tL_R)
            tau_tL_R = tau_tL
            tau_f_R = tau_f
            tau_cL_R = tau_cL
            pi_tL_R = pi_tL

        m0 = m0_R*(1 + alfa)/(1 + alfa_R)*P0*pi_r*pi_d*pi_cL*pi_cH/(P0_R*pi_r_R*pi_d_R*pi_cL_R*pi_cH_R)*(Tt4_R/Tt4)**(1/2) #kg/s
        f = (tau_lambda - tau_r*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*T0) - tau_lambda) #kgFuel/kgAir
        T9_T0 = tau_lambda*tau_tH*tau_tL/(Pt9_P9**((gamma_t - 1)/gamma_t))*cp_c/cp_t
        V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
        T19_T0 = tau_r*tau_f/(Pt19_P19**((gamma_c - 1)/gamma_c))
        V19_a0 = M19*(T19_T0)**(1/2)
        P19_P0 = Pt19_P0/(1 + (gamma_t - 1)/2*M19**2)**(gamma_t/(gamma_t - 1))
        P9_P0 = Pt9_P0/(1 + (gamma_c - 1)/2*M9**2)**(gamma_c/(gamma_c - 1))
        FF_m0 = alfa/(1 + alfa)*a0*(V19_a0 - M0 + T19_T0/V19_a0*(1 - 1/P19_P0)/gamma_c) #N/(kg/s)
        FC_m0 = 1/(1 + alfa)*a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - 1/P9_P0)/gamma_c) #N/(kg/s)
        F_m0 = FF_m0  + FC_m0 #N/(kg/s)
        S = f/((1 + alfa)*F_m0) #(kgFuel/s)/N
        F = m0*F_m0 #N
        N_NR_fan = (T0*tau_r/(T0_R*tau_r_R)*(pi_f**((gamma_c - 1)/gamma_c) - 1)/(pi_f_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
        N_NR_H = (T0*tau_r*tau_cL/(T0_R*tau_r_R*tau_cL_R)*(pi_cH**((gamma_c - 1)/gamma_c) - 1)/(pi_cH_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
        eta_T = a0**2*((1 + f)*V9_a0**2 + alfa*(V19_a0**2)- (1 + alfa)*M0**2)/(2*f*hpr)
        eta_P = 2*V0*(1 + alfa)*F_m0/(a0**2*((1 + f)*V9_a0**2 + alfa*V19_a0**2 - (1 + alfa)*M0**2))
        eta_Total = eta_P*eta_T
        mf = S*F #kgFuel/s vazao de combustivel
        AF = 1/f #kgAir/kgFuel
        mC = m0*1/(1 + alfa) #kgAir/s vazao de ar pelo Core
        mF = m0*alfa/(1 + alfa) #kgAir/s vazao de ar pelo Fan


    def ideal_ramjet(self, M0, gamma, cp, hpr, Tt4):
        """
        Description: This method calculates the on design parameters of an ramjet turbojet engine.

        Arguments:
            M0: Mach number                             [  -  ]
            gamma: Ratio of specific heats              [kJ/kgK]
            cp: Specific heat at constant pressure      [ J/K ]
            hpr: Low heating value of fuel              [kJ/kg]
            Tt4: Total temperature leaving the burner   [  K  ]

        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
        """

        R = (gamma - 1)/gamma*cp # J/(kg.K)

        a0 = (gamma*R*self.T0)**(1/2) #m/s

        tau_r = 1 + (gamma - 1)/2*M0**2

        tau_lambda = Tt4/self.T0
        V9_a0 = M0 * (tau_lambda/tau_r)**0.5

        F_m0 = a0 * (V9_a0 - M0)
        f = (cp * self.T0)/hpr * (tau_lambda - tau_r)
        S = f/F_m0

        eta_T = 1 - 1/(tau_r)
        eta_P = 2/((tau_lambda/tau_r)**0.5 + 1)
        eta_Total = eta_P*eta_T

        return (F_m0, f, S, eta_T, eta_P, eta_Total)

def main():
    test = AircraftEngines(5000)

    test.ideal_turbojet()

if __name__ == '__main__': main()