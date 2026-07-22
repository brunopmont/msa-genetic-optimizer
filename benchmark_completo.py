import os
import subprocess
import random
import copy
import re
import time  # NOVA IMPORTAÇÃO PARA CRONOMETRAGEM

# =====================================================================
# 1. FUNÇÕES DO ALGORITMO GENÉTICO (Segurança Biológica Garantida)
# =====================================================================
def ler_fasta_alinhado(caminho):
    cabecalhos, sequencias, seq_atual = [], [], ""
    try:
        with open(caminho, 'r') as f:
            for linha in f:
                linha = linha.strip()
                if linha.startswith(">"):
                    cabecalhos.append(linha)
                    if seq_atual: sequencias.append(seq_atual); seq_atual = ""
                elif linha: seq_atual += linha
            if seq_atual: sequencias.append(seq_atual)
        return cabecalhos, sequencias
    except FileNotFoundError:
        return [], []

def salvar_resultado_fasta(cabecalhos, alinhamento, saida):
    with open(saida, 'w') as f:
        for i in range(len(alinhamento)):
            f.write(f"{cabecalhos[i]}\n{alinhamento[i]}\n")

def calcular_fitness(alinhamento):
    score, tamanho, num_seqs = 0, len(alinhamento[0]), len(alinhamento)
    for col in range(tamanho):
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                l1, l2 = alinhamento[i][col], alinhamento[j][col]
                if l1 != '-' and l2 != '-' and l1 == l2: score += 1
                elif l1 != '-' or l2 != '-': score -= 0.5
    return score

def crossover_genetico(pai_1, pai_2):
    tamanho = min(len(pai_1[0]), len(pai_2[0]))
    if tamanho < 3: return pai_1, pai_2
    ponto_corte = random.randint(1, tamanho - 1)
    filho_1, filho_2 = [], []
    
    for i in range(len(pai_1)):
        esq_1 = pai_1[i][:ponto_corte]
        letras_esq_1 = len(esq_1.replace('-', ''))
        vistas, corte_p2 = 0, 0
        if letras_esq_1 > 0:
            for j, char in enumerate(pai_2[i]):
                if char != '-': vistas += 1
                if vistas == letras_esq_1: corte_p2 = j + 1; break
        filho_1.append(esq_1 + pai_2[i][corte_p2:])
        
        esq_2 = pai_2[i][:ponto_corte]
        letras_esq_2 = len(esq_2.replace('-', ''))
        vistas, corte_p1 = 0, 0
        if letras_esq_2 > 0:
            for j, char in enumerate(pai_1[i]):
                if char != '-': vistas += 1
                if vistas == letras_esq_2: corte_p1 = j + 1; break
        filho_2.append(esq_2 + pai_1[i][corte_p1:])
        
    m1 = max(len(s) for s in filho_1)
    filho_1 = [s.ljust(m1, '-') for s in filho_1]
    m2 = max(len(s) for s in filho_2)
    filho_2 = [s.ljust(m2, '-') for s in filho_2]
    return filho_1, filho_2

def mutacao(alinhamento, taxa=0.05):
    mutado = copy.deepcopy(alinhamento)
    for i in range(len(mutado)):
        if random.random() < taxa:
            seq = list(mutado[i])
            for _ in range(10):
                pos = random.randint(0, len(seq) - 2)
                if (seq[pos] == '-' and seq[pos+1] != '-') or (seq[pos] != '-' and seq[pos+1] == '-'):
                    seq[pos], seq[pos+1] = seq[pos+1], seq[pos]
                    break
            mutado[i] = "".join(seq)
    return mutado

def extrair_qscore(saida_terminal):
    match_q = re.search(r'Q=([0-9.]+)', saida_terminal)
    match_tc = re.search(r'TC=([0-9.]+)', saida_terminal)
    q = match_q.group(1) if match_q else "Erro"
    tc = match_tc.group(1) if match_tc else "Erro"
    return q, tc

# =====================================================================
# 2. MOTOR DE AUTOMAÇÃO (Força Bruta + Cronometragem)
# =====================================================================
if __name__ == "__main__":
    
    instancias = [
        "BB11001", "BB11010", "BB12012", "BB12022", "BB20010", 
        "BB20030", "BB30016", "BB30024", "BB40011", "BB50014"
    ]
    
    resultados_finais = []
    
    pasta_logs = "logs_individuais"
    pasta_base = "base"
    pasta_ag = "genetic-alg"
    pasta_resultados = "resultados"
    
    os.makedirs(pasta_logs, exist_ok=True)
    os.makedirs(pasta_base, exist_ok=True)
    os.makedirs(pasta_ag, exist_ok=True)
    os.makedirs(pasta_resultados, exist_ok=True)
    
    print("Iniciando o Benchmark Automatizado Completo...\n")
    print(f"Estrutura montada: Pastas prontas. Executando em FORÇA TOTAL com Cronometragem.\n")
    
    for seq in instancias:
        print(f"\n[{seq}] Processando instância...")
        
        # Inicia o cronômetro para esta instância
        tempo_inicio = time.time()
        
        caminho_in = f"bench/bench1.0/bali3/in/{seq}"
        caminho_ref = f"bench/bench1.0/bali3/ref/{seq}"
        
        arquivo_base = os.path.join(pasta_base, f"{seq}_base.fasta")
        arquivo_ag = os.path.join(pasta_ag, f"{seq}_ag.fasta")
        arquivo_log_execucao = os.path.join(pasta_logs, f"log_execucao_{seq}.txt")
        
        # 1. KAlign: Gera o Seed Alignment
        subprocess.run(["kalign", "-i", caminho_in, "-o", arquivo_base], capture_output=True)
        
        # 2. Algoritmo Genético
        cabecalhos, base_seqs = ler_fasta_alinhado(arquivo_base)
        if base_seqs:
            num_seqs = len(base_seqs)
            tamanho = len(base_seqs[0])
            
            print(f"  -> Dimensões: {num_seqs} fitas x {tamanho} colunas.")
            
            geracoes_atuais = 200
            pop_atual = 40
                
            print(f"  -> Iniciando AG ({geracoes_atuais} gerações). Acompanhe o log!")
            
            with open(arquivo_log_execucao, "w") as f_log:
                f_log.write(f"=== INICIANDO EXECUÇÃO FORÇA TOTAL: {seq} ===\n")
                f_log.write(f"Dimensões: {num_seqs} fitas x {tamanho} colunas\n")
                f_log.write(f"Parâmetros: {geracoes_atuais} Gerações | Pop: {pop_atual}\n\n")
                
                populacao = [base_seqs] + [mutacao(base_seqs, 0.40) for _ in range(pop_atual - 1)]
                
                for g in range(geracoes_atuais):
                    populacao_avaliada = [(calcular_fitness(ind), ind) for ind in populacao]
                    populacao_avaliada.sort(reverse=True, key=lambda x: x[0])
                    nova_populacao = [populacao_avaliada[0][1], populacao_avaliada[1][1]]
                    
                    if g % 5 == 0:
                        linha_log = f"Geração {g:03d}/{geracoes_atuais} | Melhor WSP: {populacao_avaliada[0][0]:.2f}\n"
                        print("     ..." + linha_log.strip())
                        f_log.write(linha_log)
                        f_log.flush()
                        
                    while len(nova_populacao) < pop_atual:
                        metade_pop = max(2, pop_atual // 2)
                        p1 = random.choice(populacao_avaliada[:metade_pop])[1]
                        p2 = random.choice(populacao_avaliada[:metade_pop])[1]
                        f1, f2 = crossover_genetico(p1, p2)
                        
                        nova_populacao.extend([mutacao(f1, 0.10), mutacao(f2, 0.10)])
                        
                    populacao = nova_populacao[:pop_atual]
                
                f_log.write("\n=== ALGORITMO GENÉTICO FINALIZADO ===\n")
            
            melhor_alinhamento = populacao_avaliada[0][1]
            salvar_resultado_fasta(cabecalhos, melhor_alinhamento, arquivo_ag)
            
            # 3. Avaliação no QScore
            comando_qscore = ["./qscore", "-seqdiffwarn", "-ignoretestcase", "-test", arquivo_ag, "-ref", caminho_ref, "-truncname"]
            processo = subprocess.run(comando_qscore, capture_output=True, text=True)
            q, tc = extrair_qscore(processo.stdout)
            
            # Para o cronômetro e calcula a duração
            tempo_fim = time.time()
            tempo_total_segundos = tempo_fim - tempo_inicio
            
            # Salva na lista com formatação de 2 casas decimais
            resultados_finais.append((seq, q, tc, f"{tempo_total_segundos:.2f}"))
            
            print(f"  -> Concluído em {tempo_total_segundos:.2f} segundos! Q-Score: {q} | TC: {tc}")
            
            caminho_resultado_ind = os.path.join(pasta_resultados, f"resultado_final_{seq}.txt")
            with open(caminho_resultado_ind, "w") as f_ind:
                f_ind.write(f"Instância: {seq}\nQ-Score: {q}\nTC-Score: {tc}\nTempo de Execução (s): {tempo_total_segundos:.2f}\n")
        else:
            print(f"  -> Erro: Falha ao ler {arquivo_base}.")

    # =====================================================================
    # 3. EXPORTAÇÃO DA TABELA CONSOLIDADA (AGORA COM TEMPO)
    # =====================================================================
    print("\n" + "="*70)
    print(" TABELA FINAL DE RESULTADOS (MODO FORÇA BRUTA)")
    print("="*70)
    print("| Instância | Q-Score (AG) | TC-Score (AG) | Tempo (s) |")
    print("| :--- | :--- | :--- | :--- |")
    for res in resultados_finais:
        print(f"| **{res[0]}** | {res[1]} | {res[2]} | {res[3]} |")
    print("="*70)
    
    caminho_tabela_geral = os.path.join(pasta_resultados, "resultados_gerais_benchmark.txt")
    with open(caminho_tabela_geral, "w") as f_geral:
        f_geral.write("TABELA DE COMPARAÇÃO DE RESULTADOS - AG FORÇA TOTAL\n")
        f_geral.write("="*70 + "\n")
        f_geral.write(f"{'Instância':<15} | {'Q-Score':<10} | {'TC-Score':<10} | {'Tempo (s)':<15}\n")
        f_geral.write("-" * 70 + "\n")
        for res in resultados_finais:
            f_geral.write(f"{res[0]:<15} | {res[1]:<10} | {res[2]:<10} | {res[3]:<15}\n")
        f_geral.write("="*70 + "\n")
        
    print(f"\n[SUCESSO] Pipeline concluído. O compilado final está seguro na pasta '{pasta_resultados}/'.")