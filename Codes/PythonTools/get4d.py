# https://www.biostars.org/p/157527/

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

def translate(seq):
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''

    for I in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break

    return peptide

codon_dict={}
for aa in set(amino_acids):
    temp=[]
    for c in codons:
       if translate(c)==aa:
           temp.append(c)
    if len(temp)==4:
        codon_dict[aa]=temp

print(codon_dict)
