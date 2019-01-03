"""
Rules to filter passaged viruses except for reference strains.
"""

rule filter_passaged_viruses:
    message:
        """
        Filtering to exclude passaged viruses
        """
    input:
        sequences = "data/h3n2_ha.fasta",
        metadata = "filtering_metadata.tsv",
        excluded = "outliers/non_reference_passaged_strains.txt"
    output:
        sequences = "data/h3n2_ha_unpassaged.fasta"
    conda: "envs/anaconda.python3.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.excluded} \
            --output {output.sequences}
        """

rule build_simple_metadata_for_filtering:
    input: "data/h3n2_ha.fasta"
    output: "filtering_metadata.tsv"
    shell: """grep "^>" {input} | sed 's/>//' | sort | uniq | awk '{{ if (NR == 1) {{ print "strain" }} print }}' > {output}"""

rule remove_reference_strains_from_passaged_strains:
    message: "Removing reference strains from list of passaged strains"
    input:
        passaged="outliers/passaged_strains.txt",
        references="data/reference_viruses.txt"
    output: "outliers/non_reference_passaged_strains.txt"
    run:
        with open(input["references"], "r") as fh:
            references = set([line.rstrip() for line in fh])

        with open(input["passaged"], "r") as fh:
            with open(output[0], "w") as oh:
                for line in fh:
                    # Parse out strain name from pipe-delimited values.
                    strain, rest = line.split("|", 1)

                    # If this exact strain is not a reference, write it out.
                    if strain not in references:
                        oh.write(line)

rule find_passaged_strains:
    input: "data/h3n2_ha.fasta"
    output: "outliers/passaged_strains.txt"
    shell: """grep "^>" {input} | grep -v -E "\|(cell|unpassaged)\|" | sed 's/>//' | sort | uniq > {output}"""
