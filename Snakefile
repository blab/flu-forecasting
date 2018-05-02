# Imports.
import pandas as pd

localrules: clean, augur_prepare

# Load configuration parameters.
configfile: "config.json"

YEAR_RANGES = config["year_ranges"]
VIRUSES = config["viruses"]
PREDICTORS = config["predictors"]
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = range(NUMBER_OF_SAMPLES)

def _get_start_year_from_range(wildcards):
    return wildcards["year_range"].split("-")[0]

def _get_end_year_from_range(wildcards):
    return wildcards["year_range"].split("-")[1]

def _get_predictor_list(wildcards):
    return " ".join(wildcards["predictors"].split("-"))

rule all:
    input: "model_accuracy.tab", "model_parameters.tab"

rule aggregate_model_parameters:
    input: expand("model_parameters/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_parameters.tab"
    run:
        df = pd.concat([pd.read_table(i) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule aggregate_model_accuracy:
    input: expand("model_accuracy/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_accuracy.tab"
    run:
        df = pd.concat([pd.read_table(i) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule summarize_model:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    output:
        accuracy="model_accuracy/{year_range}/{viruses}/{predictors}/{sample}.tab",
        parameters="model_parameters/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input[0], "r") as fh:
            model = json.load(fh)

        accuracy = model["accuracy"]
        df = pd.DataFrame({
            "year_range": [wildcards.year_range],
            "viruses": [wildcards.viruses],
            "predictors": [wildcards.predictors],
            "sample": [wildcards.sample],
            "rho_rel": [accuracy["rho_rel"]],
            "clade_error": [accuracy["clade_error"]]
        })
        df.to_csv(output["accuracy"], sep="\t", index=False)

        df = pd.DataFrame(model["params"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output["parameters"], sep="\t", index=False)

rule run_fitness_model:
    input:
        tree="dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        frequencies="dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_frequencies.json",
        titers="dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    params: predictor_list=_get_predictor_list
    conda: "envs/anaconda.python2.yaml"
    shell: "python fit_model.py {input.tree} {input.frequencies} {output} {params.predictor_list} --titers {input.titers}"

rule augur_process:
    input: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    output: "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
            "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_frequencies.json"
    conda: "envs/anaconda.python2.yaml"
    shell: """cd dist/augur/builds/flu && python flu.process.py -j ../../../../{input} --no_mut_freqs --tree_method raxml"""

rule augur_prepare:
    input: sequences="dist/fauna/data/h3n2_ha.fasta"
    output: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    conda: "envs/anaconda.python2.yaml"
    params: start_year=_get_start_year_from_range, end_year=_get_end_year_from_range
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v {wildcards.viruses} --sequences ../../../../{input.sequences} \
  --file_prefix flu_h3n2_ha_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample} --lineage h3n2 --segment ha --time_interval {params.start_year}-10-01 {params.end_year}-04-01 --sampling even"""

rule download_sequences_and_titers:
    output: "dist/fauna/data/h3n2_ha.fasta", "dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    conda: "envs/anaconda.python2.yaml"
    shell: "cd dist/fauna && python download_all.py --virus flu --flu_lineages h3n2 --segments ha --sequences --titers"

rule clean:
    run:
        print(config.get("clean"))
        if config.get("clean") is not None:
            items = config.get("clean").split(",")
            for year_range in YEAR_RANGES:
                for virus in VIRUSES:
                    for item in items:
                        if item == "all":
                            shell("rm -f dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                        elif item == "auspice":
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
                        elif item == "tree":
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*tree*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/processed/flu_h3n2_ha_{year_range}y_{virus}v*newick*".format(year_range=year_range, virus=virus))
                            shell("rm -f dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{virus}v*.json".format(year_range=year_range, virus=virus))
