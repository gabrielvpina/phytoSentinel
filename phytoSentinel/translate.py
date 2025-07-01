from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys

def safe_translate(seq, table):
    """
    Traduz uma sequência de DNA em proteína, códons inválidos viram 'X'.
    """
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = str(seq[i:i+3])
        try:
            aa = str(Seq(codon).translate(table=table))
        except TranslationError:
            aa = 'X'
        except Exception:
            aa = 'X'
        protein.append(aa)
    return Seq("".join(protein))

def translate_genome_6frames_with_ambiguity(input_fasta_file, output_fasta_file):
    """
    Traduz um arquivo FASTA genômico completo em todas as 6 frames de leitura possíveis.
    Códons ambíguos são traduzidos como 'X'.
    """
    print(f"Iniciando a tradução de '{input_fasta_file}' em 6 quadros de leitura...")

    genetic_code = CodonTable.unambiguous_dna_by_id[1]
    output_records = []
    num_sequences_processed = 0

    try:
        for record in SeqIO.parse(input_fasta_file, "fasta"):
            seq_obj = record.seq
            seq_id = record.id
            num_sequences_processed += 1

            print(f"Processando sequência: {seq_id} (tamanho: {len(seq_obj)} bp)")

            # Fita direta (frames +1, +2, +3)
            for i in range(3):
                translated_seq = safe_translate(seq_obj[i:], table=genetic_code)
                output_records.append(SeqRecord(
                    translated_seq,
                    id=f"{seq_id}_frame_{i+1}",
                    description=f"Translation of {seq_id} (forward) in frame {i+1}"
                ))

            # Fita reversa complementar (frames -1, -2, -3)
            rev_comp_seq_obj = seq_obj.reverse_complement()
            for i in range(3):
                translated_seq_rc = safe_translate(rev_comp_seq_obj[i:], table=genetic_code)
                output_records.append(SeqRecord(
                    translated_seq_rc,
                    id=f"{seq_id}_rev_frame_{i+1}",
                    description=f"Translation of {seq_id} (reverse complement) in frame {i+1}"
                ))

        # Escreve no arquivo de saída
        SeqIO.write(output_records, output_fasta_file, "fasta")
        print(f"\nTradução concluída com sucesso! {len(output_records)} sequências de proteína salvas em '{output_fasta_file}'.")

    except FileNotFoundError:
        print(f"Erro: O arquivo '{input_fasta_file}' não foi encontrado.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Ocorreu um erro inesperado durante a tradução: {e}", file=sys.stderr)
        sys.exit(1)

# Execução principal
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python translate.py <arquivo_genoma.fasta> <arquivo_saida.fasta>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    translate_genome_6frames_with_ambiguity(input_file, output_file)
