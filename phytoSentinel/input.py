import os, subprocess, shutil
from Bio import SeqIOnumberSeqs



def validateFasta(inputContig):
    """
    check the fasta file - only full nucleotide FASTA 
    """
    try:

        records = SeqIO.parse(inputContig, "fasta")

        # empty file
        first_record = next(records, None)
        if first_record is None:
            print("\nError: The Fasta file is empty or don't have valid content.")
            return False

        # validate first record
        if not all(c in "ACGTN" for c in first_record.seq.upper()):
            print(f"\nErro: The sequence '{first_record.id}' have invalid characters.")
            return False

        # validate all records
        for record in records:
            if not all(c in "ACGTN" for c in record.seq.upper()):
                print(f"\nError: The sequence '{record.id}' have invalid characters.")
                return False
        return True

    except FileNotFoundError:
        print("Error: Fasta file not found.")
        return False
    except Exception as e:
        print(f"Error to process FASTA file: {e}")
        return False


def createOutDir(inputContig, vvFolder):
    # directory w/ results
    if not os.path.exists(vvFolder):
        os.makedirs(vvFolder)
        infile = os.path.normpath(inputContig)
        outfile = os.path.join(vvFolder, inputContig)

    shutil.copy(infile, outfile)


def countFasta(vvFolder):
    # count sequences in fasta
    input_fasta = None
    for fasta in os.listdir(vvFolder):
        if fasta.endswith(".fasta"):
            input_fasta = os.path.join(vvFolder, fasta)
            break

    numberSeqs = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))

    return numberSeqs