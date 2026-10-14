# Calling-BAM filter validation

The default calling BAM retains reads with MAPQ >= 30, excludes flags 3844, and requires the proper-pair flag only for paired-end samples. Single-end samples never receive `-f 0x2`.

Before using these defaults for a new species or strain panel, compare at least three representative isolates (low divergence, high divergence, and unusual coverage/ploidy) under:

1. the default settings;
2. MAPQ 20 versus 30;
3. proper-pair filtering enabled versus disabled for paired libraries;
4. soft-clipped bases enabled versus disabled in HaplotypeCaller.

Record mapped read counts, retained calling-BAM depth, callable genome fraction, SNP/INDEL counts, Ti/Tv where meaningful, genotype missingness, and concordance at trusted loci. Keep the strict defaults only when they do not systematically remove evidence in divergent, repetitive, aneuploid, or structurally variable regions.
