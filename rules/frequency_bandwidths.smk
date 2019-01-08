"""
Rules to calculate fitness model accuracy by frequency bandwidth.
"""

rule estimate_frequencies_by_bandwidth:
    message:
        """
        Estimating frequencies by bandwidth for {input.tree}
          - narrow bandwidth: {wildcards.bandwidth}
        """
    input:
        tree=rules.refine.output.tree,
        metadata=rules.parse.output.metadata,
        weights="data/region_weights.json"
    output:
        frequencies="results/frequencies_by_bandwidth/{bandwidth}/flu_{lineage}_{segment}_{year_range}y_{viruses}v_{sample}.json"
    params:
        pivot_frequency=config["frequencies"]["pivot_frequency"],
        start_date=_get_start_date_from_range,
        end_date=_get_end_date_from_range
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/estimate_frequencies_by_bandwidth_{lineage}_{segment}_{bandwidth}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/estimate_frequencies_by_bandwidth_{lineage}_{segment}_{bandwidth}_{year_range}y_{viruses}v_{sample}.log"
    shell: """python3 scripts/frequencies.py {input.tree} {input.metadata} {output} \
--narrow-bandwidth {wildcards.bandwidth} \
--proportion-wide 0.0 \
--pivot-frequency {params.pivot_frequency} \
--start-date {params.start_date} \
--end-date {params.end_date} \
--weights {input.weights} \
--weights-attribute region \
--include-internal-nodes \
--censored &> {log}"""

rule run_fitness_model_by_frequency_bandwidth:
    input:
        ha_tree="results/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_tree.json",
        ha_metadata="results/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_meta.json",
        ha_sequences="results/auspice/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}_seq.json",
        frequencies="results/frequencies_by_bandwidth/{bandwidth}/flu_h3n2_ha_{year_range}y_{viruses}v_{sample}.json",
        masks="config/ha_masks.tsv"
    output:
        model="results/models_by_bandwidth/{bandwidth}/{year_range}/{viruses}/{sample}.json",
        validation_data_frame="results/model_data_frames_by_bandwidth/{bandwidth}/{year_range}/{viruses}/validation_{sample}.tsv"
    params:
        predictor_list="lbi",
        min_freq=config["fitness_model"]["min_freq"],
        max_freq=config["fitness_model"]["max_freq"],
        min_training_window=config["fitness_model"]["min_training_window"]
    conda: "../envs/anaconda.python3.yaml"
    benchmark: "benchmarks/fitness_model_by_bandwidth_{bandwidth}_{year_range}y_{viruses}v_{sample}.txt"
    log: "logs/fitness_model_by_bandwidth_{bandwidth}_{year_range}y_{viruses}v_{sample}.log"
    shell: "python3 src/fit_model.py {input.ha_tree} {input.ha_metadata} {input.ha_sequences} {input.frequencies} {output.model} {params.predictor_list} --masks {input.masks} --validation-data-frame {output.validation_data_frame} --min-freq {params.min_freq} --max-freq {params.max_freq} -v --min-training-window {params.min_training_window} &> {log}"

rule label_model_validation_by_bandwidth_table:
    input:
        rules.run_fitness_model_by_frequency_bandwidth.output.validation_data_frame
    output:
        data_frame="results/model_data_frames_by_bandwidth/{bandwidth}/{year_range}/{viruses}/labeled_validation_{sample}.tsv"
    run:
        df = pd.read_table(input[0])
        df["year_range"] = wildcards.year_range
        df["viruses"] = wildcards.viruses
        df["bandwidth"] = wildcards.bandwidth
        df["sample"] = wildcards.sample
        df.to_csv(output["data_frame"], sep="\t", header=True, index=False, na_rep="N/A")

rule aggregate_model_validation_by_bandwidth:
    input: expand(rules.label_model_validation_by_bandwidth_table.output.data_frame, year_range=YEAR_RANGES, viruses=VIRUSES, sample=SAMPLES, bandwidth=BANDWIDTHS)
    output: "results/model_validation_by_bandwidth.tab"
    run:
        df = pd.concat([pd.read_table(i, keep_default_na=False, na_values="N/A") for i in input], ignore_index=True, sort=True)
        df.to_csv(output[0], sep="\t", index=False, na_rep="N/A")
