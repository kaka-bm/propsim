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
            M0: Mach number
            gamma: Ratio of specific heats
            cp: Specific heat at constant pressure
            hpr: Low heating value of fuel
            Tt4: Total temperature leaving the burner
            pi_c: Compressor total pressure ratio

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
            gamma_c: Ratio of specific heats in the compressor
            gamma_t: Ratio of specific heats in the turbine
            cp_c: Specific heat at constant pressure in the compressor
            cp_t: Specific heat at constant pressure in the turbine
            hpr: Low heating value of fuel [kJ/kg]
            Tt4: Total temperature leaving the burner
            pi_c: Compressor total pressure ratio
            pi_d_max: Difuser maximum total pressure ratio
            pi_b: Burner total pressure ratio
            pi_n: Nozzle total pressure ratio
            e_t: Polytropic compressor efficiency
            e_t: Polytropic turbine efficiency
            eta_b: Combustor efficiency
            eta_m: Mechanical efficiency
            P0_P9: Ratio between Po/P9


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
            M0: Mach number
            gamma: Ratio of specific heats
            cp: Specific heat at constant pressure
            hpr: Low heating value of fuel
            Tt4: Total temperature leaving the burner
            pi_c: Compressor total pressure ratio
            
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


    def real_turbofan(self, M0, gamma, cp, hpr, Tt4, pi_c, pi_f, alpha):
        """
        Description: This method calculates the on design parameters of an ideal turbofan engine.

        Arguments:
            M0: Mach number
            gamma: Ratio of specific heats
            cp: Specific heat at constant pressure
            hpr: Low heating value of fuel
            Tt4: Total temperature leaving the burner
            pi_c: Compressor total pressure ratio
            
        Returns: A tuple containing the following outputs
            F_m0: Specific Thrust           [0]
            f: Fuel Air ratio               [1]
            S: Specific fuel consumption    [2]
            eta_T: Thermal efficiency       [3]
            eta_P: Propulsive efficiency    [4]
            eta_Total: Total efficiency     [5]
            FR: Thrust ratio                [6]
        """   

    def ideal_ramjet(self, M0, gamma, cp, hpr, Tt4):
        """
        Description: This method calculates the on design parameters of an ramjet turbojet engine.

        Arguments:
            M0: Mach number
            gamma: Ratio of specific heats [kJ/kgK]
            cp: Specific heat at constant pressure
            hpr: Low heating value of fuel
            Tt4: Total temperature leaving the burner [kJ/kg]

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