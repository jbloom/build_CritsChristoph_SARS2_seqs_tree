"""Filter alignment and mask specifies sites."""


import Bio.SeqIO
import Bio.Seq

import pandas as pd

mask_sites = snakemake.params.mask_sites
max_unmasked_n = snakemake.params.max_unmasked_n
drop_strains = snakemake.params.drop_strains

metadata = pd.read_csv(snakemake.input.metadata, sep="\t").assign(
    strain=lambda x: x["strain"].str.replace("/", "_")
)
strains = set(metadata["strain"])

retained_seqs = []
found_strains = set()
for seq in Bio.SeqIO.parse(snakemake.input.alignment, "fasta"):
    strain = seq.id.replace("/", "_")
    assert strain in strains
    if (strain in drop_strains) or (seq.id in drop_strains):
        print(f"Dropping {strain} as it is in `drop_strains`")
        continue
    s = list(str(seq.seq))
    for i in mask_sites:
        s[i - 1] = "N"
    s = "".join(s)
    unmasked_n = sum(nt in {"N", "n"} for (i, nt) in enumerate(seq.seq) if i + 1 not in mask_sites)
    if unmasked_n > max_unmasked_n:
        print(f"Dropping {strain} because it has {unmasked_n} ambiguous nucleotides at unmasked sites")
        continue
    retained_seqs.append((strain, s))
    found_strains.add(strain)

assert len(retained_seqs) == len(found_strains)

print(f"Retaining {len(retained_seqs)} sequences.")
with open(snakemake.output.alignment, "w") as f:
    for head, seq in retained_seqs:
        f.write(f">{head}\n{seq}\n")

metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)
