import matplotlib.pyplot as plt
# DNA Properties Tool
print('''\n Hey Wassup!! Choose the action you want to perform with your DNA / Protein sequence:\n
-------> 1> AT-GC Content of DNA Sequence
-------> 2> Hydrophobicity Profile of Amino Acid Sequence
-------> 3> Six Reading Frames of Amino Acid Sequences from DNA Sequence
-------> 4> Molecular weight of Protein Sequence \n''')
input_method = input('Click the respective number(1-4) of the function you want to apply: ')
#Getting six reading Frames-GACATTGTGAACAGTAAAAAAGTCCATGCAATGCGCAAGGAGCAGAAGAGGAAGCAGGGCAAGCAGCGCTCCATGGGCTCTCCCATGGACTACTCTCCTCTGCCCATCGACAAGCATGAGCCTGAATTTGGTCCATGCAGAAGAAAACTGGATGGG
def six_reading_frames():
    dna_seq = input("Add your DNA sequence here to get Six reading frames OF Amino Acid Sequences: ").upper()
    codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
    print("Six Reading Frames of Amino acid Sequences are: ")
    aa_seqs=['','','','','','']
    for i in range(0,len(dna_seq),3):
        if dna_seq[i:i+3] in codon_table:
            aa_seqs[0]+=codon_table[dna_seq[i:i+3]]
    for i in range(1,len(dna_seq),3):
        if dna_seq[i:i+3] in codon_table:
            aa_seqs[1]+=codon_table[dna_seq[i:i+3]]
    for i in range(2,len(dna_seq),3):
        if dna_seq[i:i+3] in codon_table:
            aa_seqs[2]+=codon_table[dna_seq[i:i+3]]
    # For complementary sequence
    comp_seq=''
    for char in dna_seq[::-1]:
        if char == 'A':
            comp_seq+='T'
        elif char == 'T':
            comp_seq+='A'
        elif char == 'G':
            comp_seq+='C'
        elif char == 'C':
            comp_seq+='G'
    for i in range(0,len(comp_seq),3):
        if comp_seq[i:i+3] in codon_table:
            aa_seqs[3]+=codon_table[comp_seq[i:i+3]]
    for i in range(1,len(comp_seq),3):
        if comp_seq[i:i+3] in codon_table:
            aa_seqs[4]+=codon_table[comp_seq[i:i+3]]
    for i in range(2,len(comp_seq),3):
        if comp_seq[i:i+3] in codon_table:
            aa_seqs[5]+=codon_table[comp_seq[i:i+3]]
    print(f'''
          Three Forward Frames: \n
          sequence 1 ---> {aa_seqs[0]}
          sequence 2 ---> {aa_seqs[1]}
          sequence 3 ---> {aa_seqs[2]} \n
          Three backward Frames: \n
          sequence 4 ---> {aa_seqs[3]}
          sequence 5 ---> {aa_seqs[4]}
          sequence 6 ---> {aa_seqs[5]}''')

#Getting AT-GC content of DNA sequence
def at_gc_content():
    dna_seq = input("Add your DNA sequence here to get AT and GC contents: ").upper()
    counter={'A':0,'T':0,'G':0,'C':0}
    for char in dna_seq:
        if char in counter:
            counter[char]+=1
    print(counter)
    at_content = round(((counter['A']+counter['T'])/len(dna_seq))*100, 2)
    gc_content = round(((counter['G']+counter['C'])/len(dna_seq))*100,2)
    print(f"AT Content: {at_content}%\nGC Content: {gc_content}%")
#Getting Hydrophobicity Profile of Amino Acid Sequence
def hydrophobicity_profile():
    hydrophobicity_scale = {
        'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
        'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
        'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
        'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
    }
    aa_seq = input("Add your Amino Acid sequence here to get Hydrophobicity Profile: ").upper()
    hydrophobicity_values = [hydrophobicity_scale.get(aa, 0) for aa in aa_seq]
    plt.figure(figsize=(10, 5))
    plt.plot(range(len(hydrophobicity_values)), hydrophobicity_values, marker='o',  color='b')
    plt.xticks(range(len(hydrophobicity_values)), list(aa_seq), rotation=45)
    plt.axhline(0, color='gray', linestyle='--')
    plt.xlabel('Amino Acids')
    plt.ylabel('Hydrophobicity Value')
    plt.title('Hydrophobicity Profile of Amino Acid Sequence')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    print("Hydrophobicity Profile:", hydrophobicity_values)
#Getting Molecular weight of Protein Sequence
def molecular_weight():
    aa_seq = input("Add your Protein sequence here to get Molecular weight: ").upper()
    aa_weights = {
        'A': 89.09, 'C': 121.16, 'D': 133.10, 'E': 147.13,
        'F': 165.19, 'G': 75.07, 'H': 155.16, 'I': 131.17,
        'K': 146.19, 'L': 131.17, 'M': 149.21, 'N': 132.12,
        'P': 115.13, 'Q': 146.15, 'R': 174.20, 'S': 105.09,
        'T': 119.12, 'V': 117.15, 'W': 204.23, 'Y': 181.19
    }
    mw = sum(aa_weights.get(aa, 0) for aa in aa_seq) - 18*(len(aa_seq)-1)  # Subtracting 18 for water molecules (H2O) in peptide bonds
    mw = round(mw, 2)  # Rounding to two decimal places
    print(f"Molecular Weight of the Protein Sequence: {mw} Da")

# six_reading_frames()
# at_gc_content()
# hydrophobicity_profile()
# molecular_weight()
if input_method == '1':
    at_gc_content()
elif input_method == '2':
    hydrophobicity_profile()
elif input_method == '3':
    six_reading_frames()
elif input_method == '4':
    molecular_weight()
