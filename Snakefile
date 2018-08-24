# Imports.
import json
import pandas as pd

# Set snakemake directory
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

localrules: clean, download_sequences_and_titers, augur_prepare, summarize_model, aggregate_model_parameters, aggregate_model_accuracy, aggregate_tree_plots

# Load configuration parameters.
configfile: "config.json"

YEAR_RANGES = config["year_ranges"]
VIRUSES = config["viruses"]
PREDICTORS = config["predictors"]
NUMBER_OF_SAMPLES = config["number_of_samples"]
SAMPLES = range(NUMBER_OF_SAMPLES)

def _get_start_date_from_range(wildcards):
    return "%s-10-01" % wildcards["year_range"].split("-")[0]

def _get_end_date_from_range(wildcards):
    return "%s-04-01" % wildcards["year_range"].split("-")[1]

def _get_predictor_list(wildcards):
    return " ".join(wildcards["predictors"].split("-"))

rule all:
    input: "model_accuracy.tab", "model_parameters.tab", "trees.pdf", "models.tab", "model_fold_change.pdf"

rule trees:
    input: "trees.pdf"

rule aggregate_tree_plots:
    input: expand("figures/trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES)
    output: "trees.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

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

rule aggregate_models:
    input: expand("models/{year_range}/{viruses}/{predictors}/{sample}.tab", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "models.tab"
    run:
        df = pd.concat([pd.read_table(i) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule aggregate_model_fold_change:
    input: expand("model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf", year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule plot_model_fold_change:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    output: "model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf"
    conda: "envs/anaconda.python2.yaml"
    shell: "python plot_model_fold_change.py {input} {output} {wildcards.year_range} {wildcards.viruses} {wildcards.predictors} {wildcards.sample}"

rule convert_model_json_to_tsv:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input[0], "r") as fh:
            model_json = json.load(fh)

        df = pd.DataFrame(model_json["data"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample

        df.to_csv(output[0], sep="\t", header=True, index=False)

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
            "correlation_rel": [accuracy["correlation_rel"]],
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
        frequencies="frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        titers="dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    output: "models/{year_range}/{viruses}/{predictors}/{sample}.json"
    params: predictor_list=_get_predictor_list
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/fitness_model_{year_range}y_{viruses}v_{sample}_{predictors}.txt"
    log: "logs/fitness_model_{year_range}y_{viruses}v_{sample}_{predictors}.log"
    shell: "python fit_model.py {input.tree} {input.frequencies} {output} {params.predictor_list} --titers {input.titers} &> {log}"

rule estimate_frequencies:
    input:
        tree="dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        weights="region_weights.json"
    output: "frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/estimate_frequencies_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/estimate_frequencies_{year_range}y_{viruses}v_{sample}.log"
    shell: """python scripts/frequencies.py {input.tree} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--start-date {params.start_date} \
--end-date {params.end_date} \
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes &> {log}"""

rule plot_tree:
    input: "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
    output: "figures/trees/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.pdf"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/plot_tree_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/plot_tree_{year_range}y_{viruses}v_{sample}.log"
    shell: """cd dist/augur/scripts && python plot_tree.py {SNAKEMAKE_DIR}/{input} {SNAKEMAKE_DIR}/{output} &> {SNAKEMAKE_DIR}/{log}"""

rule augur_process:
    input: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    output: "dist/augur/builds/flu/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/augur_process_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/augur_process_{year_range}y_{viruses}v_{sample}.log"
    shell: """cd dist/augur/builds/flu && python flu.process.py -j ../../../../{input} --no_mut_freqs --no_tree_freqs --tree_method raxml --export_translations &> {SNAKEMAKE_DIR}/{log}"""

rule augur_prepare:
    input: sequences="dist/fauna/data/h3n2_ha.fasta"
    output: "dist/augur/builds/flu/prepared/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/augur_prepare_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/augur_prepare_{year_range}y_{viruses}v_{sample}.log"
    params: start_date=_get_start_date_from_range, end_date=_get_end_date_from_range
    shell: """cd dist/augur/builds/flu && python flu.prepare.py -v {wildcards.viruses} --sequences ../../../../{input.sequences} \
  --file_prefix flu_h3n2_ha_{wildcards.year_range}y_{wildcards.viruses}v_{wildcards.sample} --lineage h3n2 --segment ha --time_interval {params.start_date} {params.end_date} --sampling even -r 12y &> {SNAKEMAKE_DIR}/{log}"""

rule download_sequences_and_titers:
    output: "dist/fauna/data/h3n2_ha.fasta", "dist/fauna/data/h3n2_public_hi_cell_titers.tsv"
    conda: "envs/anaconda.python2.yaml"
    benchmark: "benchmarks/download_sequences_and_titers.txt"
    log: "logs/download_sequences_and_titers.log"
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
