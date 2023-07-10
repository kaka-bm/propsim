import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def convergent_divergent_nozzle_design(mdot8, Pt8, Tt8, Pt9_Pt8, gamma, CD, is_supersonic, P0):
    output = {
    'A9_A8': [],
    'V9i': [],
    'V9': [],
    'CV': [],
    'Cfg': [],
    'Fg': []
    }
    R = 1716
    MFP=0.5224
    
    areas = [1, 1.25 , 1.5, 1.75, 2, 2.23, 2.5, 2.78, 3, 3.3, 3.5, 3.74, 4, 4.5, 5, 6, 7, 8, 9, 10 ]
    machs = [1, 1.59, 1.83, 2, 2.16, 2.26, 2.38, 2.48, 2.56, 2.62, 2.7, 2.76, 2.84, 2.94, 3.04, 3.22, 3.36, 3.48, 3.6, 3.7]
    pressoes = [0.54, 0.25, 0.169, 1.29, 0.1, 0.085, 0.069, 0.059, 0.052, 0.047, 0.041, 0.037, 0.033, 0.028, 0.023,
                0.017, 0.014, 0.011, 0.0099, 0.0085 ]
    for i in range(len(areas)):
        A9_A8 = areas[i]      
        
        Pt9 = Pt9_Pt8*Pt8

        ## Eq 1
        A8e = (mdot8*math.sqrt(Tt8))/(Pt8*MFP)
        A8 = CD*A8e
        A9 = A8*A9_A8

        ## 2
        A_A_star_9i = (1/CD)*A9_A8

        ## 4
        P9i_Pt9 = (1 + ((gamma - 1) / 2) * machs[i]**2)**-(gamma/(gamma-1)) - 0.0052
        #P9i_Pt9 = pressoes[i] - 0.0052

        ## 5
        P9i = P9i_Pt9*Pt8
        ## 6
        V9i = math.sqrt(R * Tt8) * math.sqrt(abs((2 * gamma / (gamma - 1)) * (1 - (P9i / Pt8)**((gamma - 1) / gamma))))
        ## 7
        A_Astar_9 = (Pt9 / Pt8) * (A9_A8) * (1 / CD)
        A_Astar_9 = areas[i]

        ## 8
        M92 = ( (A_Astar_9)**( 2*(gamma-1)/(gamma+1) )*( (2/(gamma+1))**-1 ) - 1 ) / ( (gamma-1)/2  )
        M9 = math.sqrt(M92)  
        M9 = machs[i]

        ## 9
        P9_Pt9 =  (1 + ((gamma - 1) / 2) * M9**2)**-(gamma/(gamma-1))
        P9_Pt9 = (1 + ((gamma - 1) / 2) * machs[i]**2)**-(gamma/(gamma-1)) #pressoes[i]

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
    
    df_otimo = df[df.Fg == df['Fg'].max()]
    print('Razão de área ótima para maior empuxo e respectivos valores de velocidade, CV, Cfg e Fg')
    display(df_otimo)
    #print('Valor ótimo da razão de área para maior empuxo: ')
    
    fig, axes = plt.subplots(1, 4, figsize=(20, 4))
    sns.set_palette(palette = sns.dark_palette("navy", reverse=True))  # Define a paleta de cores como tons de azul
    #plt.title("Gross Thrust, Velocity, CV and Cfg Vs Nozzle Area Ratio")  # Define o título do gráfic
    sns.lineplot(ax=axes[0],data=df, x="A9_A8", y="Fg")
    axes[0].set_title("Gross Thrust Vs Nozzle Area Ratio")
    sns.lineplot(ax=axes[1],data=df, x="A9_A8", y="V9")
    axes[1].set_title("Velocity Vs Nozzle Area Ratio") 
    sns.lineplot(ax=axes[2],data=df, x="A9_A8", y="CV")
    axes[2].set_title("CV Vs Nozzle Area Ratio") 
    sns.lineplot(ax=axes[3],data=df, x="A9_A8", y="CV")
    axes[3].set_title("Cfg Vs Nozzle Area Ratio") 
    
    print('Gráficos com variação dos valores em função da razão das áreas')
    
    plt.show()