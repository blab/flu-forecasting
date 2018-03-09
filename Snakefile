
BUILD = "flu_h3n2_ha_15y_17v"

rule all:
    input: "dist/augur/builds/flu/auspice/%s_meta.json" % BUILD

rule augur_process:
    input: "dist/augur/builds/flu/prepared/{build}.json"
    output: "dist/augur/builds/flu/auspice/{build}_meta.json"
    conda: "envs/anaconda.python2.yaml"
    shell: """cd dist/augur/builds/flu && python flu.process.py -j ../../../../{input} --no_mut_freqs"""

rule augur_prepare:
    input: sequences="dist/fauna/data/h3n2_ha.fasta", titers="dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    output: "dist/augur/builds/flu/prepared/{build}.json"
    conda: "envs/anaconda.python2.yaml"
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v 17 --sequences ../../../../{input.sequences} --titers ../../../../{input.titers} \
  --file_prefix {wildcards.build} --lineage h3n2 --segment ha --sampling even --time_interval 2002-10-01 2017-04-01"""

rule download_sequences_and_titers:
    output: "dist/fauna/data/h3n2_ha.fasta", "dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    conda: "envs/anaconda.python2.yaml"
    shell: "cd dist/fauna && python download_all.py --virus flu --flu_lineages h3n2 --segments ha --sequences --titers"
