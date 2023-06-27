import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def nozzle_design(mdot8, Pt8, Tt8, Pt9_Pt8, gamma, MFP, CD, M0, P0):
    output = {
    'A9_A8': [],
    'V9i': [],
    'V9': [],
    'CV': [],
    'Cfg': [],
    'Fg': []
    }
    R = 1716
    
    upper = 4
    lower = 1
    step = 0.1
    for i in range(int((upper - lower) / step) + 1):
        A9_A8 = lower + i * step
        
        Pt9 = Pt9_Pt8*Pt8

        ## Eq 1
        A8e = (mdot8*math.sqrt(Tt8))/(Pt8*MFP)
        A8 = CD*A8e
        A9 = A8*A9_A8

        ## 2
        A_A_star_9i = A9/(CD*A8)

        ## 4
        P9i_Pt9 = (1+((gamma - 1)/2)*M0**2)**(-1) #!!!!!!
        P9i_Pt9 = 0.099

        ## 5
        P9i = P9i_Pt9*Pt8
        ## 6
        V9i = math.sqrt(R * Tt8) * math.sqrt((2 * gamma / (gamma - 1)) * (1 - (P9i / Pt8)**((gamma - 1) / gamma)))
        ## 7
        A_Astar_9 = (Pt9 / Pt8) * (A9 / A8) * (1 / CD)

        ## 8 utilizando método da bisseção tb para encontrar M9
        ## A_Astar_9 = (2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M9**2)) ** ((gamma + 1) / (2 * (gamma - 1))
        #M9 = get_M9(M9_lower = 0.1, M9_upper = 5.0, tolerance = 1e-6)  
        M92 = ( (A_Astar_9)**( 2*(gamma-1)/(gamma+1) )*( (2/(gamma+1))**-1 ) - 1 ) / ( (gamma-1)/2  )
        M9 = math.sqrt(M92)  
        M9 = 2.146

        ## 9
        P9_Pt9 =  (1 + ((gamma - 1) / 2) * M9**2)**-(gamma/(gamma-1))
        P9_Pt9 = 0.1025

        ##10
        P9 = P9_Pt9 * Pt8
        ##11                                                                     
        V9 = math.sqrt(R * Tt8) * math.sqrt((2 * gamma / (gamma - 1)) * (1 - (P9 / Pt8)**((gamma - 1) / gamma)))
        ##12
        CV = V9 / V9i
        ##13                                                                     
        Cfg = CD * CV * math.sqrt((1 - (P9i / Pt8)**((gamma - 1) / gamma)) / (1 - (P0 / Pt8)**((gamma - 1) / gamma))) * (1 + ((gamma - 1) / (2 * gamma)) * ((1 - (P0 / P9)) / ((Pt9 / P9)**((gamma - 1) / gamma) - 1)))
        ##14
        Fg = mdot8*V9/32.174 + (P9-P0)*A9

        output['V9i'].append(V9i)
        
        output['A9_A8'].append(A9_A8)
        output['V9'].append(V9)
        output['CV'].append(CV)
        output['Cfg'].append(Cfg)
        output['Fg'].append(Fg)
        
    df = pd.DataFrame(output)
    sns.set_style("darkgrid")  # Define o fundo do gráfico com grid
    sns.set_palette(palette = sns.dark_palette("navy", reverse=True))  # Define a paleta de cores como tons de azul
    plt.title("Gross Thrust Vs Nozzle Area Ratio")  # Define o título do gráfic
    sns.scatterplot(data=df, x="A9_A8", y="Fg")
    plt.show()

    return()