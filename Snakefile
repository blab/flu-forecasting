
BUILD = "flu_h3n2_ha_15y_17v"
YEAR_RANGES = ("2002-2017",)
VIRUSES = ("17",)
NUMBER_OF_SAMPLES = 1
SAMPLES = range(NUMBER_OF_SAMPLES)

def _get_start_year_from_range(wildcards):
    return wildcards["year_range"].split("-")[0]

def _get_end_year_from_range(wildcards):
    return wildcards["year_range"].split("-")[1]

rule all:
    input: expand("dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)

rule run_fitness_model:
    input: "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json"
    output: "models/{year_range}/{viruses}/{sample}.json"
    conda: "envs/anaconda.python2.yaml"
    shell: "python "


rule augur_process:
    input: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    output: "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json"
    conda: "envs/anaconda.python2.yaml"
    shell: """cd dist/augur/builds/flu && python flu.process.py -j ../../../../{input} --no_mut_freqs --tree_method fasttree"""

rule augur_prepare:
    input: sequences="dist/fauna/data/h3n2_ha.fasta", titers="dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    output: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    conda: "envs/anaconda.python2.yaml"
    params: start_year=_get_start_year_from_range, end_year=_get_end_year_from_range
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v {wildcards.viruses} --sequences ../../../../{input.sequences} --titers ../../../../{input.titers} \
  --file_prefix flu_h3n2_ha_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample} --lineage h3n2 --segment ha --sampling even --time_interval {params.start_year}-10-01 {params.end_year}-04-01"""

rule download_sequences_and_titers:
    output: "dist/fauna/data/h3n2_ha.fasta", "dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    conda: "envs/anaconda.python2.yaml"
    shell: "cd dist/fauna && python download_all.py --virus flu --flu_lineages h3n2 --segments ha --sequences --titers"

rule clean:
    run:
        for year_range in YEAR_RANGES:
            for virus in VIRUSES:
                shell("rm -f dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{virus}v.json".format(year_range=year_range, virus=virus))
                shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
