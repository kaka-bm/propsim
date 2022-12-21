from ambiance import Atmosphere
import math


class AircraftEngines:

    def __init__(self, height):
        self.atm = Atmosphere(height)
        self.T0 = self.atm.temperature[0]
        self.P0 = self.atm.pressure[0]
        self.a0 = self.atm.speed_of_sound[0]

    def ideal_turbojet(self, M0, gamma, cp, hpr, Tt4, pi_c, batch_size=1, min_pi_c=0.001, max_pi_c=40):
        """
        Description: This method calculates the on design parameters of an ideal turbojet engine.

        Arguments:
            M0: Mach number                                     [  -  ]
            gamma: Ratio of specific heats                      [kJ/kgK]
            cp: Specific heat at constant pressure              [ J/K ]
            hpr: Low heating value of fuel                      [kJ/kg]
            Tt4: Total temperature leaving the burner           [  K  ]
            pi_c: Compressor total pressure ratio               [  -  ]
            batch_size: Number of points to iterate pi_c        [  -  ]
            min_pi_c: Min value for pi_c (only used in batch)   [  -  ]
            max_pi_c: Max value for pi_c (only used in batch)   [  -  ]

        Returns: A dictionary containing the list of calculated outputs for each batch.
            pi_c: Compressor total pressure ratio   
            F_m0: Specific Thrust                   
            f: Fuel Air ratio                       
            S: Specific fuel consumption            
            eta_T: Thermal efficiency               
            eta_P: Propulsive efficiency            
            eta_Total: Total efficiency             
        """

        output = {
            'pi_c': [],
            'F_m0': [],
            'f': [],
            'S': [],
            'eta_T': [],
            'eta_P': [],
            'eta_Total': []
        }

        pi_c_increase = 1

        if batch_size <= 0:
            return output
        elif batch_size == 1:
            max_pi_c = pi_c  
        else:
            pi_c = min_pi_c
            pi_c_increase = (max_pi_c - min_pi_c)/batch_size

        while pi_c <= max_pi_c:
        
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

            if math.isnan(F_m0) or math.isnan(S) or math.isnan(f) or math.isnan(eta_P) or math.isnan(eta_T) or math.isnan(eta_Total):
                pi_c += pi_c_increase
                continue

            output['pi_c'].append(pi_c)
            output['F_m0'].append(F_m0)
            output['f'].append(f)
            output['S'].append(S)
            output['eta_T'].append(eta_T)
            output['eta_P'].append(eta_P)
            output['eta_Total'].append(eta_Total)

            pi_c += pi_c_increase
        
        return output

    def real_turbojet(self, M0, gamma_c, gamma_t, cp_c, cp_t, hpr, Tt4 , pi_c, pi_d_max, pi_b, pi_n, e_c, e_t, eta_b, eta_m, P0_P9, batch_size=1, min_pi_c=0.001, max_pi_c=40):
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
            batch_size: Number of points to iterate pi_c                [  -  ]
            min_pi_c: Min value for pi_c (only used in batch)           [  -  ]
            max_pi_c: Max value for pi_c (only used in batch)           [  -  ]

        Returns: A dictionary containing the list of calculated outputs for each batch.
            pi_c: Compressor total pressure ratio   
            F_m0: Specific Thrust                   
            f: Fuel Air ratio                       
            S: Specific fuel consumption            
            eta_T: Thermal efficiency               
            eta_P: Propulsive efficiency            
            eta_Total: Total efficiency   
        """  

        output = {
            'pi_c': [],
            'F_m0': [],
            'f': [],
            'S': [],
            'eta_T': [],
            'eta_P': [],
            'eta_Total': []
        }

        pi_c_increase = 1

        if batch_size <= 0:
            return output
        elif batch_size == 1:
            max_pi_c = pi_c  
        else:
            pi_c = min_pi_c
            pi_c_increase = (max_pi_c - min_pi_c)/batch_size

        while pi_c <= max_pi_c:

            # Gas constants
            R_c    = ((gamma_c - 1)/gamma_c ) * (cp_c)
            R_t    = ((gamma_t - 1)/gamma_t ) * (cp_t)

            # Free stream
            a0    = math.sqrt(gamma_c*R_c*self.T0)
            V0    = a0 * M0

            # Free flow parameters
            tal_r    = 1 + ((gamma_c - 1)/2 ) * (M0**2)
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

            # Auxiliar
            Pt9_P9 = (P0_P9)*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n
            M9     = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
            T9_T0  = (tal_lambda*tal_t)/( ((Pt9_P9)**((gamma_t-1)/gamma_t)))*(cp_c/cp_t)   
            V9_a0  = M9*( math.sqrt( gamma_t*R_t*(T9_T0)/(gamma_c*R_c) ) )

            # Results
            F_m0      = a0*((1+f)*V9_a0 - M0 + (1+f)*R_t*T9_T0/(R_c*V9_a0)*(1-P0_P9)/gamma_c)
            S         = f/F_m0
            eta_T     = ((a0**2)*( (1+f)*V9_a0**2 - M0**2 )/(2*f*(hpr)))
            eta_P     = ( 2*V0*F_m0/( (a0**2)*( (1+f)*V9_a0**2 - M0**2 ) ) )
            eta_Total = eta_T*eta_P

            if math.isnan(F_m0) or math.isnan(S) or math.isnan(f) or math.isnan(eta_P) or math.isnan(eta_T) or math.isnan(eta_Total):
                pi_c += pi_c_increase
                continue

            output['pi_c'].append(pi_c)
            output['F_m0'].append(F_m0)
            output['f'].append(f)
            output['S'].append(S)
            output['eta_T'].append(eta_T)
            output['eta_P'].append(eta_P)
            output['eta_Total'].append(eta_Total)

            pi_c += pi_c_increase

        return output

    def real_turbojet_off_design(self,
        M0,
        Tt4,
        P0_P9,
        # Constantes
        gamma_c,
        cp_c,
        gamma_t,
        cp_t,
        hpr,
        pi_d_max,
        pi_b,
        pi_t,
        pi_n,
        tau_t,
        eta_c,
        eta_b,
        eta_m,

        # Condições de referência
        M0_R,
        T0_R,
        P0_R,
        tau_r_R,
        pi_r_R,
        Tt4_R,
        pi_d_R,
        pi_c_R,
        tau_c_R,
        Pt9_P9_R,

        #  Inputs extras do Rolls-Royce Nene
        m0_R = 50 # kg/s (?)
        ):
        
        T0 = self.T0
        P0 = self.P0

        Tt2_R = T0_R*tau_r_R


        #  Equations
        R_c = (gamma_c - 1)/gamma_c*cp_c # J/(kg.K)
        R_t = (gamma_t - 1)/gamma_t*cp_t # J/(kg.K)
        a0 = (gamma_c*R_c*T0)**(1/2) # m/s
        V0 = a0*M0
        tau_r = 1 + (gamma_c - 1)/2*M0**2
        pi_r = tau_r**(gamma_c/(gamma_c - 1))
        if M0 <= 1:
            eta_r = 1
        else:
            eta_r = 1 - 0.075*(M0 - 1)**1.35

        pi_d = pi_d_max*eta_r
        Tt2 = T0*tau_r
        tau_c = 1 + (tau_c_R - 1)*Tt4/Tt2/(Tt4_R/Tt2_R)
        pi_c = (1 + eta_c*(tau_c - 1))**(gamma_c/(gamma_c - 1))
        tau_lambda = cp_t*Tt4/(cp_c*T0)
        f = (tau_lambda - tau_r*tau_c)/(hpr*eta_b/(cp_c*T0) - tau_lambda) # kgFuel/kgAir
        m0 = m0_R*P0*pi_r*pi_d*pi_c/(P0_R*pi_r_R*pi_d_R*pi_c_R)*(Tt4_R/Tt4)**(1/2) # kg/s
        Pt9_P9 = P0_P9*pi_r*pi_d*pi_c*pi_b*pi_t*pi_n
        M9 = (2/(gamma_t - 1)*(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
        T9_T0 = tau_lambda*tau_t/(Pt9_P9**((gamma_t - 1)/gamma_t))*cp_c/cp_t
        V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
        F_m0 = a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - P0_P9)/gamma_c) # N/(kg/s)
        F = F_m0*m0 # N
        S = f/F_m0 # (kgFuel/s)/N
        eta_T = a0**2*((1 + f)*V9_a0**2 - M0**2)/(2*f*hpr)
        eta_P = 2*V0*F_m0/(a0**2*((1 + f)*V9_a0**2 - M0**2))
        eta_Total = eta_P*eta_T
        N_NR = (T0*tau_r/(T0_R*tau_r_R)*(pi_c**((gamma_c - 1)/gamma_c) - 1)/(pi_c_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
        A9_A9R = (Pt9_P9/Pt9_P9_R)**((gamma_t + 1)/(2*gamma_t))*((Pt9_P9_R**((gamma_t - 1)/gamma_t) - 1)/(Pt9_P9**((gamma_t - 1)/gamma_t) - 1))**(1/2)
        mc2_mc2_R = pi_c/pi_c_R*((Tt4_R/Tt2_R)/(Tt4/Tt2))**(1/2) # vazão mássica corrigida no compressor

        #  Outputs extras
        V9 = V9_a0*a0 # m/s
        AF = 1/f  # kgAir/kgFuel
        Pt4 = P0*pi_r*pi_d*pi_c*pi_b # Pa
        Pt9 = P0*pi_r*pi_d*pi_c*pi_b*pi_t # Pa
        T9 = T0*T9_T0 # K
        
        

    def ideal_turbofan(self, M0, gamma, cp, hpr, Tt4, pi_c, pi_f, alpha, batch_size=1, min_pi_c=0.001, max_pi_c=40):
        """
        Description: This method calculates the on design parameters of an ideal turbofan engine.

        Arguments:
            M0: Mach number                                             [  -  ]
            gamma: Ratio of specific heats                              [kJ/kgK]
            cp: Specific heat at constant pressure                      [ J/K ]
            hpr: Low heating value of fuel                              [kJ/kg]
            Tt4: Total temperature leaving the burner                   [  K  ]   
            pi_c: Compressor total pressure ratio                       [  -  ]
            batch_size: Number of points to iterate pi_c                [  -  ]
            min_pi_c: Min value for pi_c (only used in batch)           [  -  ]
            max_pi_c: Max value for pi_c (only used in batch)           [  -  ]
        
        Returns: A dictionary containing the list of calculated outputs for each batch.
            pi_c: Compressor total pressure ratio   
            F_m0: Specific Thrust                   
            f: Fuel Air ratio                       
            S: Specific fuel consumption            
            eta_T: Thermal efficiency               
            eta_P: Propulsive efficiency            
            eta_Total: Total efficiency
            FR: Thrust ratio   
        """    
        output = {
            'pi_c': [],
            'F_m0': [],
            'f': [],
            'S': [],
            'eta_T': [],
            'eta_P': [],
            'eta_Total': [],
            'FR': []
        }

        pi_c_increase = 1

        if batch_size <= 0:
            return output
        elif batch_size == 1:
            max_pi_c = pi_c  
        else:
            pi_c = min_pi_c
            pi_c_increase = (max_pi_c - min_pi_c)/batch_size

        while pi_c <= max_pi_c:

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
            eta_P = 2 * M0 * (V9_a0 - M0 + alpha * (V19_a0 - M0))/((V9_a0 * self.a0)**2/(self.a0**2) - M0**2 + alpha * ((V19_a0*self.a0)**2/(self.a0**2) - M0**2))
            eta_Total = eta_P * eta_T
            FR = (V9_a0 - M0)/(V19_a0 - M0)

            if math.isnan(F_m0) or math.isnan(S) or math.isnan(f) or math.isnan(eta_P) or math.isnan(eta_T) or math.isnan(eta_Total):
                pi_c += pi_c_increase
                continue

            output['pi_c'].append(pi_c)
            output['F_m0'].append(F_m0)
            output['f'].append(f)
            output['S'].append(S)
            output['eta_T'].append(eta_T)
            output['eta_P'].append(eta_P)
            output['eta_Total'].append(eta_Total)
            output['FR'].append(FR)

            pi_c += pi_c_increase

        return output


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
        alpha,
        batch_size=1, min_pi_c=0.001, max_pi_c=40 ):

        """
        Description: This method calculates the on design parameters of a real twin spool turbofan engine.

        Arguments:
            M0: Mach number 
            gamma_c: Ratio of specific heats at the compressor
            gamma_t: Ratio of specific heats at the turbine
            cp_c: Specific heat at constant pressure at the compressor
            cp_t: Specific heat at constant pressure at the turbine
            hpr: Low heating value of fuel 
            Tt4: Total temperature leaving the burner
            pi_d_max: Maximum total pressure ratio at the difuser
            pi_b: Total pressure ratio at the burner
            pi_n: Total pressure ratio at the nozzle
            pi_fn: Total pressure ratio between the fan and the nozzle
            e_cL: Politropic efficiency of the low pressure compressor
            e_cH: : Politropic efficiency of the high pressure compressor
            e_f: Politropic efficiency of the fan   
            e_tL: Politropic efficiency of the low pressure turbine
            e_tH: Politropic efficiency of the high pressure turbine
            eta_b: Burner efficiency
            eta_mL: Mechanical efficiency for low pressure
            eta_mH: Mechanical efficiency for high pressure
            P0_P9: Pressure ratio between the freestream and the gas generator outlet
            P0_P19: Pressure ratio between the freestream and the fan outlet
            tau_n: Total temperature ratio at the nozzle
            tau_fn: Total temperature ratio between the fan and the nozzle
            pi_cL: Total pressure ratio at the low pressure compressor
            pi_cH: Total pressure ratio at the high pressure compressor
            pi_f: Total pressure ratio at the fan
            alpha: Bypass ratio
        
        Returns: A dictionary containing the list of calculated outputs for each batch.
            F_m0: Specific Thrust                   
            f: Fuel Air ratio                       
            S: Specific fuel consumption            
            eta_T: Thermal efficiency               
            eta_P: Propulsive efficiency            
            eta_Total: Total efficiency
            FR: Thrust ratio 
        """
        output = {
            'F_m0': [],
            'f': [],
            'S': [],
            'eta_T': [],
            'eta_P': [],
            'eta_Total': [],
            'FR': []
        }


        R_c = (gamma_c - 1)/gamma_c*cp_c #J/(kg.K)
        R_t = (gamma_t - 1)/gamma_t*cp_t #J/(kg.K)

        a0 = (gamma_c*R_c*self.T0)**(1/2) #m/s
        V0 = a0*M0 # m/s

        # Free stream parameters
        tau_r = 1 + (gamma_c - 1)/2*M0**2

        pi_r = tau_r**(gamma_c/(gamma_c - 1))

        if M0 <= 1:
            eta_r = 1
        else:
            eta_r = 1 - 0.075*(M0 - 1)**1.35

        # Diffuser parameters
        pi_d = pi_d_max*eta_r
        tau_d = pi_d**((gamma_c - 1)/gamma_c)

        # Fan parameters
        tau_f = pi_f**((gamma_c - 1)/(gamma_c*e_f))
        eta_f = (pi_f**((gamma_c - 1)/gamma_c) - 1)/(tau_f - 1)

        # Enthalpy
        tau_lambda = cp_t*Tt4/(cp_c*self.T0)

        # Compressor parameters
        tau_cL = pi_cL**((gamma_c - 1)/(gamma_c*e_cL))
        tau_cH = pi_cH**((gamma_c - 1)/(gamma_c*e_cH))

        eta_cL = (pi_cL**((gamma_c - 1)/gamma_c) - 1)/(tau_cL - 1)
        eta_cH = (pi_cH**((gamma_c - 1)/gamma_c) - 1)/(tau_cH - 1)

        # Turbine parameters
        f = (tau_lambda - tau_r*tau_d*tau_f*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*self.T0) - tau_lambda) 

        tau_tH = 1 - (tau_cH - 1)/(1 + f)/tau_lambda*tau_r*tau_d*tau_f*tau_cL*eta_mH
        tau_tL = 1 - ((alpha*(tau_f - 1) + (tau_cL - 1))*eta_mL/(1 + f)*tau_r*tau_d/tau_lambda/tau_tH)

        pi_tH = tau_tH**(gamma_t/((gamma_t - 1)*e_tH))
        pi_tL = tau_tL**(gamma_t/((gamma_t - 1)*e_tL))

        eta_tH = (tau_tH - 1)/(pi_tH**((gamma_t - 1)/gamma_t)-1)
        eta_tL = (tau_tL - 1)/(pi_tL**((gamma_t - 1)/gamma_t)-1)

        # Engine core parameters
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

        # Engine performance parameters
        FF_m0 = alpha/(1 + alpha)*a0*(V19_a0 - M0 + T19_T0/V19_a0*(1 - P0_P19)/gamma_c) #N/(kg/s)
        FC_m0 = 1/(1 + alpha)*a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t/R_c*T9_T0/V9_a0*(1 - P0_P9)/gamma_c) #N/(kg/s)
        F_m0 = FF_m0  + FC_m0 #N/(kg/s)

        S = f/((1 + alpha)*F_m0) 

        FR = FF_m0/FC_m0

        eta_T = a0*a0*((1 + f)*V9_a0*V9_a0 + alpha*(V19_a0*V19_a0)- (1 + alpha)*M0*M0)/(2*f*hpr)
        eta_P = 2*M0*((1 + f)*V9_a0 + alpha*V19_a0 - (1 + alpha)*M0)/((1 + f)*(V9_a0**2) + alpha*V19_a0**2 - (1 + alpha)*M0**2)
        eta_Total = eta_P*eta_T

        output['F_m0'].append(F_m0)
        output['f'].append(f)
        output['S'].append(S)
        output['eta_T'].append(eta_T)
        output['eta_P'].append(eta_P)
        output['eta_Total'].append(eta_Total)
        output['FR'].append(FR)

        return output


    def real_turbofan_off_design(self, 
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

        # On-design references
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
        alpha_R,
        M9_R,
        M19_R,
        m0_R
        ):

        tau_cH_R = pi_cH_R**((gamma_c - 1)/(gamma_c))
        tau_cL_R = pi_cL_R**((gamma_c - 1)/(gamma_c))

        R_c = (gamma_c - 1)/gamma_c*cp_c 
        R_t = (gamma_t - 1)/gamma_t*cp_t 
        a0 = (gamma_c*R_c*self.T0)**(1/2) 
        V0 = a0*M0 
        tau_r = 1 + (gamma_c - 1)/2*M0**2
        pi_r = tau_r**(gamma_c/(gamma_c - 1))

        if M0 <= 1:
            eta_r = 1
        else:
            eta_r = 1 - 0.075*(M0 - 1)**1.35

        pi_d = pi_d_max*eta_r
        tau_lambda = cp_t*Tt4/(cp_c*self.T0)

        teste = 10
        while teste > 0.0001:
            tau_tL = tau_tL_R
            tau_f = tau_f_R
            tau_cL = tau_cL_R
            pi_tL = pi_tL_R
            pi_cL = pi_cL_R
            tau_cH = 1 + Tt4/self.T0/(Tt4_R/T0_R)*(tau_f_R*tau_cL_R)/(tau_r*tau_cL)*(tau_cH_R - 1)
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
            alpha = alpha_R*pi_cL_R*pi_cH_R/pi_f_R/(pi_cL*pi_cH/pi_f)*((tau_lambda)/(tau_r*tau_f)/(tau_lambda_R/(tau_r_R*tau_f_R)))**(1/2)*MFP_M19/MFP_M19_R
            tau_f = 1 + (tau_f_R - 1)*((1 - tau_tL)/(1 - tau_tL_R)*(tau_lambda/tau_r)/(tau_lambda_R/tau_r_R)*(tau_cL_R - 1 + alpha_R*(tau_f_R - 1))/(tau_cL_R - 1 + alpha*(tau_f_R - 1)))
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

        m0 = m0_R*(1 + alpha)/(1 + alpha_R)*self.P0*pi_r*pi_d*pi_cL*pi_cH/(P0_R*pi_r_R*pi_d_R*pi_cL_R*pi_cH_R)*(Tt4_R/Tt4)**(1/2) 
        f = (tau_lambda - tau_r*tau_cL*tau_cH)/(hpr*eta_b/(cp_c*self.T0) - tau_lambda) 
        T9_T0 = tau_lambda*tau_tH*tau_tL/(Pt9_P9**((gamma_t - 1)/gamma_t))*cp_c/cp_t
        V9_a0 = M9*(gamma_t*R_t/(gamma_c*R_c)*T9_T0)**(1/2)
        T19_T0 = tau_r*tau_f/(Pt19_P19**((gamma_c - 1)/gamma_c))
        V19_a0 = M19*(T19_T0)**(1/2)
        P19_P0 = Pt19_P0/(1 + (gamma_t - 1)/2*M19**2)**(gamma_t/(gamma_t - 1))
        P9_P0 = Pt9_P0/(1 + (gamma_c - 1)/2*M9**2)**(gamma_c/(gamma_c - 1))
        FF_m0 = alpha/(1 + alpha)*a0*(V19_a0 - M0 + T19_T0/V19_a0*(1 - 1/P19_P0)/gamma_c) 
        FC_m0 = 1/(1 + alpha)*a0*((1 + f)*V9_a0 - M0 + (1 + f)*R_t*T9_T0/(R_c*V9_a0)*(1 - 1/P9_P0)/gamma_c) 
        F_m0 = FF_m0  + FC_m0 
        S = f/((1 + alpha)*F_m0) 
        N_NR_fan = (self.T0*tau_r/(T0_R*tau_r_R)*(pi_f**((gamma_c - 1)/gamma_c) - 1)/(pi_f_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
        N_NR_H = (self.T0*tau_r*tau_cL/(T0_R*tau_r_R*tau_cL_R)*(pi_cH**((gamma_c - 1)/gamma_c) - 1)/(pi_cH_R**((gamma_c - 1)/gamma_c) - 1))**(1/2)
        eta_T = a0**2*((1 + f)*V9_a0**2 + alpha*(V19_a0**2)- (1 + alpha)*M0**2)/(2*f*hpr)
        eta_P = 2*V0*(1 + alpha)*F_m0/(a0**2*((1 + f)*V9_a0**2 + alpha*V19_a0**2 - (1 + alpha)*M0**2))
        eta_Total = eta_P*eta_T
        F = F_m0 * m0
        mf = S*F 
        AF = 1/f 
        mC = m0*1/(1 + alpha)
        mF = m0*alpha/(1 + alpha) 

        output = {
            'F': [F],
            'm0': [m0],
            'f': [f],
            'S': [S],
            'eta_T': [eta_T],
            'eta_P': [eta_P],
            'eta_Total': [eta_Total],
            'alpha': [alpha],
            'pi_f': [pi_f],
            'pi_cH': [pi_cH],
            'pi_tL': [pi_tL],
            'tau_f': [tau_f],
            'tau_cH': [tau_cH],
            'tau_tL': [tau_tL],
            'M9': [M9],
            'M19': [M19],
            'N_fan': [N_NR_fan],
            'N_hp_spool': [N_NR_H]
        }

        return output


    def ideal_ramjet(self, M0, gamma, cp, hpr, Tt4):
        """
        Description: This method calculates the on design parameters of an ramjet turbojet engine.

        Arguments:
            M0: Mach number                             [  -  ]
            gamma: Ratio of specific heats              [kJ/kgK]
            cp: Specific heat at constant pressure      [ J/K ]
            hpr: Low heating value of fuel              [kJ/kg]
            Tt4: Total temperature leaving the burner   [  K  ]

        Returns: A dictionary containing the list of calculated outputs.
            F_m0: Specific Thrust                   
            f: Fuel Air ratio                       
            S: Specific fuel consumption            
            eta_T: Thermal efficiency               
            eta_P: Propulsive efficiency            
            eta_Total: Total efficiency
        """

        output = {
            'F_m0': [],
            'f': [],
            'S': [],
            'eta_T': [],
            'eta_P': [],
            'eta_Total': [],
            'FR': []
        }

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

        output['F_m0'].append(F_m0)
        output['f'].append(f)
        output['S'].append(S)
        output['eta_T'].append(eta_T)
        output['eta_P'].append(eta_P)
        output['eta_Total'].append(eta_Total)

        return output

def main():
    engines = AircraftEngines(12500)

    out = engines.ideal_turbofan(M0=0.7, gamma=1.4, cp=1004, hpr=42.8*10**6, Tt4=1850, pi_c=10, pi_f=2, alpha=5, batch_size=3, min_pi_c=7, max_pi_c=15)

    print(out)

if __name__ == '__main__': main()