from ambiance import Atmosphere
import math

class sub_bocal:

    def teste_MFP(gamma,M0):
        if gamma == 1.3:
            MFP = 0.03254898 * M0**3 - 0.1088489 * M0**2 + -0.30371691 * M0 + 1.16783628
        elif gamma == 1.4:
            MFP = 0.1983 * M0**3 - 1.0612 * M0**2 + 1.5190 * M0 + 0.0183
        else:
            MFP = 0.03254898 * M0**3 - 0.1088489 * M0**2 + -0.30371691 * M0 + 1.16783628
        return MFP


    def __init__(self, height):
            self.atm = Atmosphere(height)
            self.T0 = self.atm.temperature[0]
            self.P0 = self.atm.pressure[0]
            self.a0 = self.atm.speed_of_sound[0]

    def Subsonic_inlet(self,mdot0,gamma,Me,Ms,R,T,M0=0.8):

        import math
        if gamma == 1.3:
            MFP = 0.03254898 * M0**3 - 0.1088489 * M0**2 + -0.30371691 * M0 + 1.16783628
        elif gamma == 1.4:
            MFP = 0.1983 * M0**3 - 1.0612 * M0**2 + 1.5190 * M0 + 0.0183
        else:
            MFP = 0.03254898 * M0**3 - 0.1088489 * M0**2 + -0.30371691 * M0 + 1.16783628

        MFP_at_M_08 = MFP



        # Eq 3.3
        Dt = math.sqrt((4 / math.pi) * ((mdot0 * math.sqrt(self.T0)) / self.P0) * (1 / MFP_at_M_08))

        # Eq 3.1
        At = math.pi * (Dt ** 2) / 4

        # Eq 3.4
        # Informações que estão no artigo
        Tref = 288.2  # K
        Pref = 101.3  # kPa

        mdotc0 = mdot0 * math.sqrt(self.T0 / Tref) / (self.P0 / Pref)

        # Eq 3.5
        Dt = 0.07413 * math.sqrt(mdotc0)

        # Eq 3.6
        Ve = Me

        # Eq 3.7
        Vs = Ms * self.a0

        # Eq 3.8
        atmm = self.atm
        As = Ms / (Vs*1.225)

        # Eq 3.9
        rho = (self.P0 * M0) / (R * T)

        # Eq 3.10
        theta = math.asin((Vs - Ve) / Ve)


        return Dt,mdotc0,As,rho,theta
