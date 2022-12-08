from ambiance import Atmosphere


class AircraftEngines:

    def __init__(self, height):
        __atm = Atmosphere(height)
        T0 = __atm.temperature
        P0 = __atm.pressure
        a0 = __atm.speed_of_sound

    def ideal_turbojet(self, M0, gamma, cp, hpr, Tt4, pi_c):
        """
        Description: This method calculates the on design parameters of an ideal turbojet engine.

        Arguments:
            M0: Mach number
            gamma: Ratio of specific heats
            cp: Specific heat at constant pressure
            hpr: Low heating value of fuel
            Tt4:
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


    def ideal_turbofan(self, M0, gamma, cp, hpr, Tt4, pi_c, pi_f, alpha):
        """
        Description: This method calculates the on design parameters of an ideal turbofan engine.

        Arguments:
            M0: Mach number
            gamma: Ratio of specific heats
            cp: Specific heat at constant pressure
            hpr:
            Tt4: Total temperature leaving the 
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

def main():
    test = AircraftEngines(5000)

    test.ideal_turbojet()

if __name__ == '__main__': main()