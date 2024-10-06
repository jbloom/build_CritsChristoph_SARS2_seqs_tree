"""``snakemake`` file that runs pipeline."""

# some configuration
reference = "hCoV-19/Wuhan/Hu-1/2019"  # Wuhan-Hu-1 reference
outgroup = "hCoV-19_Wuhan_WH04_2020"  # Crits-Christoph et al root at this A lineage in Fig 1A
mask_sites = list(range(1, 223)) + list(range(29700, 29904))  # mask sites near termini
max_unmasked_n = 1000  # max number of N nucleotides outside masked sites before strain dropped
drop_strains = []  # problematic dropped not caught by other filters


rule all:
    """Target rule w final outputs."""
    input:
        "auspice/build_CritsChristoph_SARS2_seqs_tree.json",


rule prep_sequences_and_metadata:
    """Prepare sequences and metadata."""
    input:
        seqs=[
            f"data/{seqset}_sequences.fa"
            for seqset in ["genbank", "gisaid", "custom", "ngdc"]
        ],
        accession_info="data/seq_metadata.csv",
    output:
        sequences="results/unaligned_sequences.fa",
        metadata="results/metadata_unfiltered.tsv",
    log:
        notebook="results/prep_sequences_and_metadata.ipynb",
    notebook:
        "notebooks/prep_sequences_and_metadata.py.ipynb"


rule align:
    """Make an alignment."""
    input:
        unaligned_seqs=rules.prep_sequences_and_metadata.output.sequences,
    output:
        alignment="results/alignment_unfiltered.fa",
    params:
        reference=reference,
    threads: 4
    shell:
        """
        augur align \
            -s {input.unaligned_seqs} \
            -o {output.alignment} \
            --reference-name {params.reference} \
            --fill-gaps \
            --nthreads {threads}
        """


rule filter_alignment:
    """Filter the alignment and mask specified sites."""
    input:
        alignment=rules.align.output.alignment,
        metadata=rules.prep_sequences_and_metadata.output.metadata,
    output:
        alignment="results/alignment_filtered.fa",
        metadata="results/metadata_filtered.tsv",
    params:
        mask_sites=mask_sites,
        max_unmasked_n=max_unmasked_n,
        drop_strains=drop_strains,
    script:
        "scripts/filter_alignment.py"


rule build_tree:
    """Build the phylogenetic tree."""
    input:
        alignment=rules.filter_alignment.output.alignment,
    output:
        tree="results/initial_tree.nwk",
        exclude_sites="results/exclude_sites.txt",
    params:
        outgroup=outgroup,
        exclude_sites="\n".join(map(str, mask_sites)),
        blmin=0.0000001,  # minimum branch length for zero-length branches
    shell:
        """
        printf "{params.exclude_sites}" > {output.exclude_sites}
        augur tree \
            -a {input.alignment} \
            --method iqtree \
            --tree-builder-args "\-o {params.outgroup} \-ninit 100 \-n 5 \-me 0.01 \-seed 1 \-blmin {params.blmin}" \
            --output {output.tree}
        """


rule collapse_zero_length_branches:
    """Collapse zero-length branches."""
    input:
        tree=rules.build_tree.output.tree,
    output:
        tree="results/collapsed_zero_branches_tree.nwk",
    params:
        blmin=rules.build_tree.params.blmin,
    script:
        "scripts/collapse_zero_length_branches.py"


rule refine_tree:
    """Run the ``augur`` refine command."""
    input:
        alignment=rules.filter_alignment.output.alignment,
        tree=rules.collapse_zero_length_branches.output.tree,
        metadata=rules.prep_sequences_and_metadata.output.metadata,
    output:
        tree="results/refined_tree.nwk",
        node_data="results/branch_lengths.json",
    shell:
        """
        augur refine \
            -a {input.alignment} \
            -t {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns strain \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --keep-root \
            --verbosity 2
        """

rule ancestral:
    """Mutations and ancestral sequences."""
    input:
        tree=rules.refine_tree.output.tree,
        alignment=rules.filter_alignment.output.alignment,
    output:
        node_data="results/nt-muts.json",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous
        """


rule export_tree:
    """Export the tree to JSON for Nextstrain auspice."""
    input:
        tree=rules.refine_tree.output.tree,
        node_data=rules.refine_tree.output.node_data,
        nt_muts=rules.ancestral.output.node_data,
        metadata=rules.filter_alignment.output.metadata,
        description="data/description.md",
    output:
        tree="auspice/build_CritsChristoph_SARS2_seqs_tree.json",
        auspice_config="results/auspice_config.json",
    shell:
        """
        printf '{{"colorings": [{{"key": "date", "title": "date", "type": "temporal"}}]}}' > {output.auspice_config}
        augur export v2 \
            --tree {input.tree} \
            --output {output.tree} \
            --auspice-config {output.auspice_config} \
            --node-data {input.node_data} {input.nt_muts} \
            --metadata {input.metadata} \
            --metadata-id-columns strain \
            --color-by-metadata date \
            --metadata-columns accession source \
            --include-root-sequence-inline \
            --panels tree \
            --title 'Phylogeny of SARS-CoV-2 sequences from Crits-Christoph et al (2024)' \
            --description {input.description}
        """
