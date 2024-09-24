# rna_to_protein.py

# Define the genetic code dictionary
genetic_code = {
    'AUG': 'M', 'UGG': 'W', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y',
    'UGU': 'C', 'UGC': 'C', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
    'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

def translate_rna_to_protein(rna_seq):
    """
    Translate an RNA sequence into its corresponding protein sequence.

    Args:
        rna_seq (str): A string of RNA bases (A, U, C, G).

    Returns:
        str: The translated protein sequence, or a message if the RNA sequence is invalid.
    """
    # Ensure the sequence is uppercase
    rna_seq = rna_seq.upper().strip()

    # Check if the sequence is empty
    if not rna_seq:
        return "Error: RNA sequence is empty."

    # Check that the sequence contains at least one complete codon
    if len(rna_seq) < 3:
        return "Error: RNA sequence too short for translation."

    # Split RNA sequence into codons (3 bases each)
    protein_seq = []
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]

        # Skip incomplete codons at the end
        if len(codon) != 3:
            print(f"Incomplete codon ignored: {codon}")
            continue

        amino_acid = genetic_code.get(codon, '-')

        # Stop translation at Stop codon
        if amino_acid == 'Stop':
            break

        protein_seq.append(amino_acid)

    return ''.join(protein_seq)

def main():
    # Example usage
    rna_seq = input("Enter an RNA sequence: ").strip()
    protein = translate_rna_to_protein(rna_seq)
    print("Protein sequence:", protein)

if __name__ == "__main__":
    main()
