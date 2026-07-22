import random
import copy
import matplotlib.pyplot as plt

# =====================================================================
# 1. LEITURA E ESCRITA DE FICHEIROS FASTA
# =====================================================================
def ler_fasta_alinhado(caminho_ficheiro):
    cabecalhos = []
    sequencias = []
    seq_atual = ""
    
    try:
        with open(caminho_ficheiro, 'r') as f:
            for linha in f:
                linha = linha.strip()
                if linha.startswith(">"):
                    cabecalhos.append(linha)
                    if seq_atual:
                        sequencias.append(seq_atual)
                        seq_atual = ""
                elif linha:
                    seq_atual += linha
            if seq_atual:
                sequencias.append(seq_atual)
        return cabecalhos, sequencias
    except FileNotFoundError:
        print(f"Erro: Não foi possível encontrar o ficheiro '{caminho_ficheiro}'.")
        return [], []

def salvar_resultado_fasta(cabecalhos, alinhamento, nome_arquivo_saida):
    with open(nome_arquivo_saida, 'w') as f:
        for i in range(len(alinhamento)):
            f.write(f"{cabecalhos[i]}\n")
            f.write(f"{alinhamento[i]}\n")
    print(f"\n[SUCESSO] O ficheiro final foi salvo como: '{nome_arquivo_saida}'")

# =====================================================================
# 2. ARTEFATOS PARA O RELATÓRIO (GRÁFICO E VISUALIZAÇÃO)
# =====================================================================
def gerar_grafico_convergencia(historico_scores, geracoes, nome_arquivo="grafico_convergencia_AG.png"):
    """Gera o gráfico de convergência provando a fuga do ótimo local."""
    plt.figure(figsize=(10, 6))
    plt.plot(range(geracoes), historico_scores, marker='', color='b', linestyle='-', linewidth=2)
    plt.title("Evolução do Score WSP - Refinamento por Algoritmo Genético", fontsize=14)
    plt.xlabel("Gerações", fontsize=12)
    plt.ylabel("Fitness (Score WSP)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Linha base do KAlign/Clustal (estagnação original)
    plt.axhline(y=historico_scores[0], color='r', linestyle='--', label='Alinhamento Base (Estagnação)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(nome_arquivo, dpi=300)
    print(f"[SUCESSO] Gráfico de convergência salvo como: '{nome_arquivo}'")

def exibir_antes_depois(alinhamento_base, alinhamento_final, janela=60):
    """Imprime um recorte visual das sequências para a análise qualitativa do relatório."""
    print("\n" + "="*60)
    print(" ANÁLISE QUALITATIVA: ANTES X DEPOIS (Para o Relatório)")
    print("="*60)
    print("-> Alinhamento Base (KAlign - Com Gaps Sub-ótimos):")
    for seq in alinhamento_base[:4]: # Mostra as 4 primeiras sequências
        print(f"   {seq[:janela]}...")
        
    print("\n-> Alinhamento Final (Algoritmo Genético - Refinado):")
    for seq in alinhamento_final[:4]:
        print(f"   {seq[:janela]}...")
    print("="*60 + "\n")

# =====================================================================
# 3. O CORAÇÃO DO ALGORITMO GENÉTICO
# =====================================================================
def calcular_fitness(alinhamento):
    score = 0
    tamanho = len(alinhamento[0])
    num_seqs = len(alinhamento)
    
    for col in range(tamanho):
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                l1 = alinhamento[i][col]
                l2 = alinhamento[j][col]
                if l1 != '-' and l2 != '-' and l1 == l2:
                    score += 1
                elif l1 == '-' and l2 == '-':
                    pass
                else:
                    score -= 0.5
    return score

def crossover_genetico(pai_1, pai_2):
    """ Crossover Seguro: Preserva a integridade da fita de ADN original. """
    tamanho = min(len(pai_1[0]), len(pai_2[0]))
    if tamanho < 3: return pai_1, pai_2
    
    ponto_corte = random.randint(1, tamanho - 1)
    filho_1, filho_2 = [], []
    
    for i in range(len(pai_1)):
        # --- FILHO 1 ---
        esq_1 = pai_1[i][:ponto_corte]
        # Conta quantas letras reais (sem gaps) tem na metade esquerda
        letras_esq_1 = len(esq_1.replace('-', ''))
        
        # Vai ao Pai 2 e extrai EXATAMENTE o resto das letras que faltam
        letras_vistas = 0
        corte_p2 = 0
        if letras_esq_1 > 0:
            for j, char in enumerate(pai_2[i]):
                if char != '-': letras_vistas += 1
                if letras_vistas == letras_esq_1:
                    corte_p2 = j + 1
                    break
        dir_2 = pai_2[i][corte_p2:]
        filho_1.append(esq_1 + dir_2)
        
        # --- FILHO 2 ---
        esq_2 = pai_2[i][:ponto_corte]
        letras_esq_2 = len(esq_2.replace('-', ''))
        
        letras_vistas = 0
        corte_p1 = 0
        if letras_esq_2 > 0:
            for j, char in enumerate(pai_1[i]):
                if char != '-': letras_vistas += 1
                if letras_vistas == letras_esq_2:
                    corte_p1 = j + 1
                    break
        dir_1 = pai_1[i][corte_p1:]
        filho_2.append(esq_2 + dir_1)
        
    # Preenche o final com gaps para garantir que a matriz do alinhamento fica quadrada
    max_len_1 = max(len(s) for s in filho_1)
    filho_1 = [s.ljust(max_len_1, '-') for s in filho_1]
    
    max_len_2 = max(len(s) for s in filho_2)
    filho_2 = [s.ljust(max_len_2, '-') for s in filho_2]
    
    return filho_1, filho_2

def mutacao(alinhamento, taxa_mutacao=0.05):
    """ Desliza um 'gap' para a posição ao lado de forma biologicamente segura. """
    alinhamento_mutado = copy.deepcopy(alinhamento)
    for i in range(len(alinhamento_mutado)):
        if random.random() < taxa_mutacao:
            seq = list(alinhamento_mutado[i])
            
            # Tenta encontrar uma posição segura para deslizar o gap (limite de 10 tentativas)
            for _ in range(10):
                pos = random.randint(0, len(seq) - 2)
                
                # Só autoriza a troca se um for letra e o outro for gap '-'
                if (seq[pos] == '-' and seq[pos+1] != '-') or (seq[pos] != '-' and seq[pos+1] == '-'):
                    seq[pos], seq[pos+1] = seq[pos+1], seq[pos]
                    break # Troca feita com sucesso, sai do loop
                    
            alinhamento_mutado[i] = "".join(seq)
    return alinhamento_mutado

# =====================================================================
# 4. CICLO DE EXECUÇÃO PRINCIPAL
# =====================================================================
def refinar_alinhamento_ag(caminho_base, arquivo_saida, geracoes=100, tam_populacao=20):
    print(f"Lendo o alinhamento base: '{caminho_base}'...")
    cabecalhos, base_seqs = ler_fasta_alinhado(caminho_base)
    if not base_seqs: return

    populacao = [base_seqs]
    for _ in range(tam_populacao - 1):
        populacao.append(mutacao(base_seqs, taxa_mutacao=0.15))
        
    print(f"Iniciando o processo evolutivo ({geracoes} gerações)...")
    
    historico_scores = [] # Guarda os scores para o gráfico
    
    for geracao in range(geracoes):
        populacao_avaliada = [(calcular_fitness(ind), ind) for ind in populacao]
        populacao_avaliada.sort(reverse=True, key=lambda x: x[0])
        
        melhor_score = populacao_avaliada[0][0]
        historico_scores.append(melhor_score)
        
        if geracao % 10 == 0 or geracao == geracoes - 1:
            print(f"  -> Geração {geracao:03d} | Melhor WSP: {melhor_score:.2f}")
            
        nova_populacao = [populacao_avaliada[0][1], populacao_avaliada[1][1]]
        
        while len(nova_populacao) < tam_populacao:
            pai_1 = random.choice(populacao_avaliada[:5])[1]
            pai_2 = random.choice(populacao_avaliada[:5])[1]
            filho_1, filho_2 = crossover_genetico(pai_1, pai_2)
            nova_populacao.extend([mutacao(filho_1), mutacao(filho_2)])
            
        populacao = nova_populacao[:tam_populacao]
        
    melhor_alinhamento_final = populacao_avaliada[0][1]
    
    # Geração dos artefatos para o relatório
    gerar_grafico_convergencia(historico_scores, geracoes)
    exibir_antes_depois(base_seqs, melhor_alinhamento_final)
    salvar_resultado_fasta(cabecalhos, melhor_alinhamento_final, arquivo_saida)

    print("\n=== PRÓXIMO PASSO (AVALIAÇÃO OFICIAL) ===")
    print(f"No seu terminal do WSL, rode o seguinte comando para obter o Q-Score:")
    print(f"qscore -test {arquivo_saida} -ref in/BB11001 -truncname")

# =====================================================================
# INÍCIO DO PROGRAMA
# =====================================================================
if __name__ == "__main__":
    FICHEIRO_ENTRADA = "bb11001_kalign_base.fasta"
    FICHEIRO_SAIDA = "bb11001_resultado_AG.fasta"
    
    refinar_alinhamento_ag(
        caminho_base=FICHEIRO_ENTRADA, 
        arquivo_saida=FICHEIRO_SAIDA, 
        geracoes=150,       
        tam_populacao=30    
    )