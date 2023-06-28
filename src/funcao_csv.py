import pandas as pd


# Definir função para obter valores correspondentes
def obter_valores(A_Ast):
    A_Ast = A_Ast * 100000
    # Especificar o caminho para o arquivo de texto
    caminho_arquivo = "C:\\Users\\carlo\\OneDrive - Universidade Federal de Minas Gerais\\2023-1\\Propulsão 2\\artigo\\código\\versao4.csv"

    # Importar o arquivo de texto como DataFrame usando o pandas
    df = pd.read_csv(caminho_arquivo,
                     delimiter=";")  # Dependendo do formato do arquivo, talvez seja necessário ajustar o delimitador

    # Renomear as colunas
    novos_nomes = {
        'col1': "M",
        "col2": 'T/Tt',
        "col3": 'P/Pt',
        "col4": 'pho/phot',
        "col5": 'A/A*',
        "col6": "MFP",
        # Adicione mais pares de colunas antigas e novas conforme necessário
    }
    df = df.rename(columns=novos_nomes)

    # Filtrar o DataFrame com base no valor de A/A*, caso o valor não seja exato esse será aproximado
    df_filtrado = df[df["A/A*"] == A_Ast]

    if df_filtrado.empty:
        print("Valor de A/A* não encontrado. Por favor, readeque o valor de entrada e rode a função novamente")
        return

    # Obter os valores correspondentes nas outras colunas
    valores_M = df_filtrado["M"].values
    valores_T_Tt = df_filtrado["T/Tt"].values
    valores_P_Pt = df_filtrado["P/Pt"].values
    valores_pho_phot = df_filtrado["pho/phot"].values
    valores_A_Ast = df_filtrado["A/A*"].values
    valores_MFP = df_filtrado["MFP"].values

    # Retornar os valores correspondentes
    return valores_M, valores_T_Tt, valores_P_Pt, valores_pho_phot, valores_A_Ast, valores_MFP


# Testando a função para valor de razão de área
valor_A_Ast = 1.45947
valores_M, valores_T_Tt, valores_P_Pt, valores_pho_phot, valores_A_Ast, valores_MFP = obter_valores(valor_A_Ast)

# Exibir os valores correspondentes
print("Valores de M:", valores_M)
print("Valores de T/Tt:", valores_T_Tt)
print("Valores de P/Pt:", valores_P_Pt)
print("Valores de pho/phot:", valores_pho_phot)
print("Valores de A/A*:", valor_A_Ast)
print("Valores de MFP:", valores_MFP)
