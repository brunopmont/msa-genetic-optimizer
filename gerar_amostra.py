import os

def extrair_primeiras_fitas(caminho, num_fitas=4, tamanho_corte=45):
    """Lê um FASTA e extrai os primeiros 'tamanho_corte' caracteres das primeiras 'num_fitas'."""
    linhas_formatadas = []
    try:
        with open(caminho, 'r') as f:
            nome_atual = ""
            seq_atual = ""
            fitas_lidas = 0
            
            for linha in f:
                linha = linha.strip()
                if linha.startswith(">"):
                    if seq_atual and fitas_lidas < num_fitas:
                        # Formata o nome para ficar alinhado (ex: Fita 1:)
                        linhas_formatadas.append(f"Fita {fitas_lidas+1}: {seq_atual[:tamanho_corte]}")
                        fitas_lidas += 1
                        if fitas_lidas == num_fitas: break
                    nome_atual = linha[1:]
                    seq_atual = ""
                else:
                    seq_atual += linha
                    
            # Captura a última fita se o arquivo acabar
            if seq_atual and fitas_lidas < num_fitas:
                linhas_formatadas.append(f"Fita {fitas_lidas+1}: {seq_atual[:tamanho_corte]}")
                
    except FileNotFoundError:
        return [f"Erro: Arquivo não encontrado -> {caminho}"]
        
    return linhas_formatadas

if __name__ == "__main__":
    # Escolha a instância que quer exibir no artigo. A BB11001 é ótima por ser curta.
    instancia_alvo = "BB11001" 
    
    # Caminhos baseados na sua estrutura de pastas
    caminho_bruto = f"bench/bench1.0/bali3/in/{instancia_alvo}"
    caminho_ag = f"genetic-alg/{instancia_alvo}_ag.fasta"
    
    print("\n" + "="*60)
    print(" BLOCO 1: ANTES")
    print("="*60)
    antes = extrair_primeiras_fitas(caminho_bruto, num_fitas=4, tamanho_corte=45)
    for linha in antes:
        # Colocando espaços entre as letras para ficar igual ao exemplo didático
        linha_espacada = linha.split(":")[0] + ":" + "".join([f" {c}" for c in linha.split(":")[1]])
        print(linha_espacada)
        
    print("\n" + "="*60)
    print(" BLOCO 2: DEPOIS")
    print("="*60)
    depois = extrair_primeiras_fitas(caminho_ag, num_fitas=4, tamanho_corte=45)
    for linha in depois:
        linha_espacada = linha.split(":")[0] + ":" + "".join([f" {c}" for c in linha.split(":")[1]])
        print(linha_espacada)
    print("="*60 + "\n")