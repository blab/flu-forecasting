"""
Rules for building fitness models from augur builds.
"""

rule estimate_frequencies:
    message:
        """
        Estimating frequencies for {input.tree}
          - narrow bandwidth: {params.narrow_bandwidth}
          - wide bandwidth: {params.wide_bandwidth}
          - proportion wide: {params.proportion_wide}
        """
    input:
        tree=rules.refine.output.tree,
        metadata=rules.parse.output.metadata,
        weights="data/region_weights.json"
    output:
        frequencies="frequencies/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.json"
    params:
        narrow_bandwidth=config["frequencies"]["narrow_bandwidth"],
        wide_bandwidth=config["frequencies"]["wide_bandwidth"],
        proportion_wide=config["frequencies"]["proportion_wide"],
        pivot_frequency=config["frequencies"]["pivot_frequency"],
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_frequencies_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/estimate_frequencies_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.log"
    shell: """python scripts/frequencies.py {input.tree} {input.metadata} {output} \
--narrow-bandwidth {params.narrow_bandwidth} \
--wide-bandwidth {params.wide_bandwidth} \
--proportion-wide {params.proportion_wide} \
--pivot-frequency {params.pivot_frequency} \
--start-date {params.start_date} \
--end-date {params.end_date} \
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes \
--censored &> {log}"""

rule run_fitness_model:
    input:
        ha_tree="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        ha_metadata="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_meta.json",
        ha_sequences="auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_seq.json",
        #na_tree="auspice/flu_h3n2_na_{year_range}y_{viruses}v_{sample}_tree.json",
        frequencies="frequencies/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        titers=rules.download_titers.output.titers,
        dms="data/dms-h3n2-preferences-rescaled.csv",
        masks="config/ha_masks.tsv"
    output:
        model="models/{year_range}/{viruses}/{predictors}/{sample}.json",
        tip_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/tips_{sample}.tsv",
        clade_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/clades_{sample}.tsv",
        validation_data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/validation_{sample}.tsv"
    params:
        predictor_list=_get_predictor_list,
        min_freq=config["fitness_model"]["min_freq"],
        max_freq=config["fitness_model"]["max_freq"],
        min_training_window=config["fitness_model"]["min_training_window"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.txt"
    log: "logs/fitness_model_{year_range}y_{viruses}v_{sample}/{predictors}.log"
    shell: "python3 src/fit_model.py {input.ha_tree} {input.ha_metadata} {input.ha_sequences} {input.frequencies} {output.model} {params.predictor_list} --titers {input.titers} --dms {input.dms} --masks {input.masks} --tip-data-frame {output.tip_data_frame} --clade-data-frame {output.clade_data_frame} --validation-data-frame {output.validation_data_frame} --min-freq {params.min_freq} --max-freq {params.max_freq} -v --min-training-window {params.min_training_window} &> {log}"

rule summarize_model:
    input:
        rules.run_fitness_model.output.model
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
            "mcc": [accuracy["mcc"]],
            "clade_error": [accuracy["clade_error"]]
        })
        df.to_csv(output["accuracy"], sep="\t", index=False, na_rep="NaN")

        df = pd.DataFrame(model["params"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output["parameters"], sep="\t", index=False)

rule convert_model_json_to_tsv:
    input:
        model=rules.run_fitness_model.output.model
    output:
        model="models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    run:
        with open(input["model"], "r") as fh:
            model_json = json.load(fh)

        df = pd.DataFrame(model_json["data"])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample

        df.to_csv(output["model"], sep="\t", header=True, index=False)

rule label_model_validation_table:
    input: rules.run_fitness_model.output.validation_data_frame
    output:
        data_frame="model_data_frames/{year_range}/{viruses}/{predictors}/labeled_validation_{sample}.tsv"
    run:
        df = pd.read_table(input[0])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["predictors"] = wildcards.predictors
        df["sample"] = wildcards.sample
        df.to_csv(output[0], sep="\t", header=True, index=False, na_rep="N/A")

rule aggregate_models:
    input: expand(rules.convert_model_json_to_tsv.output.model, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "models.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="null")

rule aggregate_model_validation:
    input: expand(rules.label_model_validation_table.output.data_frame, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "model_validation.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False, na_values="N/A") for i in input], ignore_index=True, sort=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="N/A")

rule aggregate_model_parameters:
    input:
        expand(rules.summarize_model.output.parameters, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output:
        parameters="model_parameters.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule aggregate_model_accuracy:
    input:
        expand(rules.summarize_model.output.accuracy, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output:
        accuracy="model_accuracy.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False) for i in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="null")

#
# Model parameter and accuracy plots
#

rule plot_model_parameters:
    input: rules.aggregate_model_parameters.output.parameters
    output: "figures/model_parameters.pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python scripts/plot_parameters.py {input} {output}"

rule plot_frequency_correlation:
    input: rules.aggregate_model_accuracy.output.accuracy
    output:
        correlation="figures/frequency_correlation.pdf",
        mcc="figures/mcc.pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python scripts/plot_accuracy.py {input} {output.correlation} {output.mcc}"

#
# Model fold change plots
#

rule plot_model_fold_change:
    input: "models/{year_range}/{viruses}/{predictors}/{sample}.tab"
    output:
        faceted="figures/faceted_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf",
        combined="figures/combined_model_fold_change/{year_range}/{viruses}/{predictors}/{sample}.pdf"
    conda: "../envs/anaconda.python3.yaml"
    shell: "python3 scripts/plot_model_fold_change.py {input} {output.faceted} {output.combined} {wildcards.year_range} {wildcards.viruses} {wildcards.predictors} {wildcards.sample}"

rule aggregate_combined_model_fold_change:
    input:
        expand(rules.plot_model_fold_change.output.combined, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "figures/combined_model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule aggregate_faceted_model_fold_change:
    input: expand(rules.plot_model_fold_change.output.faceted, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, predictors=PREDICTORS)
    output: "figures/faceted_model_fold_change.pdf"
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"
