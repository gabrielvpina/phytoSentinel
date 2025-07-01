from Bio.Seq import Seq
from Bio.Data import CodonTable

def translate_genome_6frames(fasta_file, output_file):
    genetic_code = CodonTable.unambiguous_dna_by_id[1] # table 1

    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_id = ""
        current_seq = ""
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    # process previous seq
                    seq_obj = Seq(current_seq)
                    for i in range(3): # Forward frames
                        prot = seq_obj[i:].translate(table=genetic_code, to_stop=False, gap="-")
                        f_out.write(f">{current_id}_frame_{i+1}\n{str(prot)}\n")

                    # Reverse frames
                    rev_comp_seq = seq_obj.reverse_complement()
                    for i in range(3):
                        prot_rc = rev_comp_seq[i:].translate(table=genetic_code, to_stop=False, gap="-")
                        f_out.write(f">{current_id}_rev_frame_{i+1}\n{str(prot_rc)}\n")
                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line

        # Processar a última sequência no arquivo
        if current_seq:
            seq_obj = Seq(current_seq)
            for i in range(3): # Forward frames
                prot = seq_obj[i:].translate(table=genetic_code, to_stop=False, gap="-")
                f_out.write(f">{current_id}_frame_{i+1}\n{str(prot)}\n")

            # Reverse frames
            rev_comp_seq = seq_obj.reverse_complement()
            for i in range(3):
                prot_rc = rev_comp_seq[i:].translate(table=genetic_code, to_stop=False, gap="-")
                f_out.write(f">{current_id}_rev_frame_{i+1}\n{str(prot_rc)}\n")


translate_genome_6frames('test.fasta', 'test_6frames.fasta')