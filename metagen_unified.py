#!/usr/bin/env python3

import subprocess
import shutil
from pathlib import Path

print("================================================================================")
print("                 OLÁ! ESTE É O PIPELINE DE METAGENÔMICA VIRAL                   ")
print("================================================================================")
print("Antes de iniciarmos, escolha o o(s) software(s) de classificação taxonômica de sua preferência.")
print("1) Kraken2")
print("2) DIAMOND")
print("3) Kraken2 & Diamond")

while True:
    escolha = input("Enter your choice (1, 2, or 3): ").strip()
    if escolha in ["1", "2", "3"]:
        break
    print("Invalid option. Please enter 1, 2, or 3.")

MAMBA = "/home/thitopaz/y/micromamba --root-prefix /home/thitopaz/micromamba run -n"

for p in ["1_raw_data", "2_clean", "3_contigs", "4_kraken", "5_daa"]:
    Path(p).mkdir(parents=True, exist_ok=True)

arquivos = sorted(list(Path(".").glob("*.fastq.gz")) + list(Path(".").glob("*.fq.gz")))
amostras = {}
for arq in arquivos:
    amostras.setdefault(arq.name.split("_")[0], []).append(arq)

for _id, files in amostras.items():
    print(f"\n[INFO] ========================================================")
    print(f"[INFO] Iniciando processamento da amostra: {_id}")
    print(f"[INFO] ========================================================")
    
    if len(files) == 2:
        r1, r2 = sorted(files)
        r1.rename(f"1_raw_data/{r1.name}")
        r2.rename(f"1_raw_data/{r2.name}")
        
        # 1. Controle de Qualidade (Fastp) - PE
        print(f" -> [ETAPA] Executando Controle de Qualidade (fastp) [Paired-End]...")
        subprocess.run(f'{MAMBA} validation_room fastp -i "1_raw_data/{r1.name}" -I "1_raw_data/{r2.name}" -o "2_clean/{_id}_clean_R1.fastq.gz" -O "2_clean/{_id}_clean_R2.fastq.gz" -h "2_clean/{_id}_report.html" -j "2_clean/{_id}_report.json" -w 12', shell=True)
        
        # 2. Remoção de rRNA (SortMeRNA) - PE
        print(f" -> [ETAPA] Executando Remoção de rRNA (SortMeRNA) [Paired-End]...")
        smr_cmd = (
            f'{MAMBA} validation_room sortmerna '
            f'--ref "/home/thitopaz/databases/smr_v4.3_default_db.fasta" '
            f'--reads "2_clean/{_id}_clean_R1.fastq.gz" --reads "2_clean/{_id}_clean_R2.fastq.gz" '
            f'--workdir "2_clean/sortmerna_{_id}" --aligned "2_clean/{_id}_rRNA" '
            f'--other "2_clean/{_id}_non_rRNA" --fastx --paired_in --out2 --threads 12 -v'
        )
        subprocess.run(smr_cmd, shell=True)
        
        # 3. Sobrescrever arquivos substituindo de forma instantânea (já que saíram em .gz)
        print(f" -> [ETAPA] Atualizando arquivos de saída com as leituras não-rRNA...")
        
        # Define os caminhos gerados pelo SortMeRNA baseados na sua imagem
        smr_out1 = Path(f"2_clean/{_id}_non_rRNA_fwd.fq.gz")
        smr_out2 = Path(f"2_clean/{_id}_non_rRNA_rev.fq.gz")
        
        target1 = Path(f"2_clean/{_id}_clean_R1.fastq.gz")
        target2 = Path(f"2_clean/{_id}_clean_R2.fastq.gz")
        
        if smr_out1.exists() and smr_out2.exists():
            # Substitui os arquivos do fastp pelos do sortmerna de forma limpa
            shutil.move(str(smr_out1), str(target1))
            shutil.move(str(smr_out2), str(target2))
        else:
            print(f"[ERRO] Arquivos de saída do SortMeRNA não foram encontrados para {_id}!")

        # Limpeza do diretório de trabalho do SortMeRNA (PE)
        print(f" -> [ETAPA] Apagando diretório de trabalho temporário do SortMeRNA...")
        shutil.rmtree(f"2_clean/sortmerna_{_id}", ignore_errors=True)
        
        # 4. Montagem de Contigs (MEGAHIT) - PE
        print(f" -> [ETAPA] Executando Montagem de Contigs (MEGAHIT) [Paired-End]...")
        subprocess.run(f'{MAMBA} assembly_room megahit -1 "2_clean/{_id}_clean_R1.fastq.gz" -2 "2_clean/{_id}_clean_R2.fastq.gz" -o "3_contigs/{_id}_assembly" -t 12', shell=True)
    
        fasta = Path(f"3_contigs/{_id}_assembly/final.contigs.fa")
        
        if fasta.exists():
            print(f" -> [ETAPA] Organizando arquivos e limpando diretório do MEGAHIT...")
            
            fasta.rename(f"3_contigs/{_id}_contigs.fasta")
            subprocess.run('mkdir -p "3_1_logs"', shell=True)
            subprocess.run(f'mv "3_contigs/{_id}_assembly/log" "3_1_logs/{_id}_megahit.log"', shell=True)
            subprocess.run(f'rm -rf "3_contigs/{_id}_assembly"', shell=True)
    
    
    
    else:
        se = files[0]
        se.rename(f"1_raw_data/{se.name}")
        
        # 1. Controle de Qualidade (Fastp) - SE
        print(f" -> [ETAPA] Executando Controle de Qualidade (fastp) [Single-End]...")
        subprocess.run(f'{MAMBA} validation_room fastp -i "1_raw_data/{se.name}" -o "2_clean/{_id}_clean_SE.fastq.gz" -h "2_clean/{_id}_report.html" -j "2_clean/{_id}_report.json" -w 12', shell=True)
        
        # 2. Remoção de rRNA (SortMeRNA) - SE
        print(f" -> [ETAPA] Executando Remoção de rRNA (SortMeRNA) [Single-End]...")
        smr_cmd = (
            f'{MAMBA} validation_room sortmerna '
            f'--ref "/home/thitopaz/databases/smr_v4.3_default_db.fasta" '
            f'--reads "2_clean/{_id}_clean_SE.fastq.gz" '
            f'--workdir "2_clean/sortmerna_{_id}" --aligned "2_clean/{_id}_rRNA" '
            f'--other "2_clean/{_id}_non_rRNA" --fastx --threads 12 -v'
        )
        subprocess.run(smr_cmd, shell=True)
        
        # 3. Sobrescrever arquivo single-end
        print(f" -> [ETAPA] Atualizando arquivo de saída com a leitura não-rRNA...")
        smr_out_se = Path(f"2_clean/{_id}_non_rRNA.fq.gz")
        target_se = Path(f"2_clean/{_id}_clean_SE.fastq.gz")
        
        if smr_out_se.exists():
            shutil.move(str(smr_out_se), str(target_se))
        else:
            print(f"[ERRO] Arquivo de saída do SortMeRNA não foi encontrado para {_id}!")

        # Limpeza do diretório de trabalho do SortMeRNA (SE)
        print(f" -> [ETAPA] Apagando diretório de trabalho temporário do SortMeRNA...")
        shutil.rmtree(f"2_clean/sortmerna_{_id}", ignore_errors=True)
        
        # 4. Montagem de Contigs (MEGAHIT) - SE
        print(f" -> [ETAPA] Executando Montagem de Contigs (MEGAHIT) [Single-End]...")
        subprocess.run(f'{MAMBA} assembly_room megahit -r "2_clean/{_id}_clean_SE.fastq.gz" -o "3_contigs/{_id}_assembly" -t 12', shell=True)
        
        fasta = Path(f"3_contigs/{_id}_assembly/final.contigs.fa")
        
        if fasta.exists():
            print(f" -> [ETAPA] Organizando arquivos e limpando diretório do MEGAHIT...")
            
            fasta.rename(f"3_contigs/{_id}_contigs.fasta")
            subprocess.run('mkdir -p "3_1_logs"', shell=True)
            subprocess.run(f'mv "3_contigs/{_id}_assembly/log" "3_1_logs/{_id}_megahit.log"', shell=True)
            subprocess.run(f'rm -rf "3_contigs/{_id}_assembly"', shell=True)

        
    # 5. Classificação Taxonômica dos Contigs Gerados
    
        
        if escolha in ["1", "3"]:
            print(f" -> [ETAPA] Executando Classificação Taxonômica com Kraken2 (Confiança: 0.05)...")
            subprocess.run(f'{MAMBA} classification_room kraken2 --db "/home/thitopaz/databases/kraken2_db" --threads 12 --use-names --confidence 0.05 --output "4_kraken/{_id}_kraken.out" --report "4_kraken/{_id}_kraken_report.txt" "3_contigs/{_id}_contigs.fasta"', shell=True)
        
        if escolha in ["2", "3"]:
            print(f" -> [ETAPA] Executando Alinhamento de Proteínas com DIAMOND (blastx)...")
            subprocess.run(f'{MAMBA} classification_room diamond blastx -d "/home/thitopaz/databases/viral_protein.dmnd" -q "3_contigs/{_id}_contigs.fasta" -o "5_daa/{_id}_contigs.daa" -t 12 --sensitive --tmpdir . -f 100', shell=True)
        
            print(f" -> [ETAPA] Convertendo resultado DAA para XML (Formato BLAST)...")
            subprocess.run(f'{MAMBA} classification_room diamond view -a "5_daa/{_id}_contigs.daa" -f 5 -o "5_daa/{_id}_contigs.xml"', shell=True)
    
            print(f" -> [ETAPA] Gerando gráfico interativo no Krona...")
            subprocess.run(f'{MAMBA} classification_room ktImportBLAST "5_daa/{_id}_contigs.xml" -o "5_daa/{_id}_krona.html"', shell=True)

    else:
        print(f"[AVISO] Arquivo final.contigs.fa não foi gerado para a amostra {_id}. Etapa de classificação ignorada.")

print("\n================================================================================")
print("                       PIPELINE FINALIZADO COM SUCESSO!                         ")
print("================================================================================")

