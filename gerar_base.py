import subprocess
import os

def gerar_alinhamento_kalign(arquivo_entrada, arquivo_saida):
    """
    Chama a ferramenta KAlign instalada no sistema para gerar o alinhamento base.
    """
    # Verifica se o arquivo cru do BAliBase existe na pasta
    if not os.path.exists(arquivo_entrada):
        print(f"[ERRO] O arquivo '{arquivo_entrada}' não foi encontrado.")
        print("Verifique se o nome está correto e se ele está na mesma pasta do script.")
        return

    print(f"A iniciar o KAlign para o arquivo cru: '{arquivo_entrada}'...")
    
    # Monta o comando que seria digitado no terminal: kalign -i entrada -o saida
    comando = ["kalign", "-i", arquivo_entrada, "-o", arquivo_saida]
    
    try:
        # Pede para o Python rodar o comando no WSL de forma invisível
        subprocess.run(comando, check=True, capture_output=True, text=True)
        print(f"[SUCESSO] Alinhamento base gerado com perfeição!")
        print(f"O ficheiro foi salvo como: '{arquivo_saida}'")
        print("Agora você já pode rodar o seu Algoritmo Genético usando este arquivo.")
        
    except subprocess.CalledProcessError as e:
        print("\n[ERRO] O KAlign falhou ao processar as sequências.")
        print("Detalhe do erro:", e.stderr)
    except FileNotFoundError:
        print("\n[ERRO] Comando 'kalign' não encontrado pelo Python.")
        print("Certifique-se de que rodou o comando 'sudo apt-get install kalign' no terminal.")

# =====================================================================
# INÍCIO DO SCRIPT
# =====================================================================
if __name__ == "__main__":
    # 1. Nome do arquivo cru que você baixou do BAliBase (pasta 'in')
    # Mude para "in/BB11001" se o arquivo estiver dentro da sub-pasta "in"
    DIRETORIO_CRU = "bench/bench1.0/bali3/in/"
    ARQUIVO_CRU = "BB11001" 
    
    # 2. O nome do arquivo que será gerado com os gaps (O seu "Seed Alignment")
    ARQUIVO_ALINHADO = "bb11001_kalign_base.fasta"
    
    gerar_alinhamento_kalign(os.path.join(DIRETORIO_CRU, ARQUIVO_CRU), ARQUIVO_ALINHADO)