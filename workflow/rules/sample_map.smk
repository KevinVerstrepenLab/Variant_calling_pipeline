# Build a GATK sample map: "<sample>\t<path_to_gvcf>"
rule make_sample_map:
    input:
        gvcfs = expand(f"{GVCF}/{{sample}}.g.vcf.gz", sample=SAMPLES),
        idxs  = expand(f"{GVCF}/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES)
    output:
        map = f"{JOINT}/samples.map"
    run:
        import os
        os.makedirs(JOINT, exist_ok=True)
        with open(output.map, "w") as fh:
            for s in SAMPLES:
                fh.write(f"{s}\t{GVCF}/{s}.g.vcf.gz\n")
