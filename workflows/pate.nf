/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pate_pipeline'

// Reference preparation
include { PREPARE_LOCUS_REFERENCE } from '../modules/local/prepare_locus_reference/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'

// Read mapping
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { PICARD_FASTQTOSAM } from '../modules/local/picard_fastqtosam/main'
include { PICARD_MERGEBAMALIGNMENT } from '../modules/local/picard_mergebamalignment/main'
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PERLOCUS } from '../modules/nf-core/samtools/index/main'

// Variant calling and filtering
include { GATK4_HAPLOTYPECALLER } from '../modules/local/gatk4_haplotypecaller/main'
include { GATK4_SELECTVARIANTS as SELECT_SNP_BIALLELIC } from '../modules/local/gatk4_selectvariants/main'
include { GATK4_VARIANTANNOTATOR } from '../modules/local/gatk4_variantannotator/main'
include { GATK4_VARIANTFILTRATION } from '../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS as SELECT_FILTERED_BIALLELIC } from '../modules/local/gatk4_selectvariants/main'
include { GATK4_FASTAALTERNATEREFERENCEMAKER } from '../modules/local/gatk4_fastaalternatereferencemaker/main'

// Phasing
include { BAMTOOLS_SPLIT } from '../modules/local/bamtools_split/main'
include { HPOPG_PHASE } from '../modules/local/hpopg_phase/main'

// Allele extraction
include { BUILD_PHASED_CONSENSUS } from '../modules/local/build_phased_consensus/main'
include { COLLECT_PHASING_STATS } from '../modules/local/collect_phasing_stats/main'

// QC subworkflows
include { BAM_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Prepare per-sample reference from per-locus FASTA files
    //
    ch_references_dir = channel.fromPath(params.references, checkIfExists: true)

    ch_sample_refs = ch_samplesheet.map { meta, reads -> [ meta, reads ] }
        .combine(ch_references_dir)
        .map { meta, reads, refs_dir -> [ meta, refs_dir ] }

    PREPARE_LOCUS_REFERENCE (
        ch_sample_refs
    )
    ch_versions = ch_versions.mix(PREPARE_LOCUS_REFERENCE.out.versions.first())

    // Per-sample reference FASTA channel
    ch_sample_reference = PREPARE_LOCUS_REFERENCE.out.fasta

    //
    // MODULE: Index per-sample reference with BWA
    //
    BWA_INDEX (
        ch_sample_reference
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions_bwa.first())

    //
    // MODULE: Create FASTA index for per-sample reference
    //
    SAMTOOLS_FAIDX (
        ch_sample_reference,
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions_samtools.first())

    //
    // MODULE: Create sequence dictionary for per-sample reference
    //
    GATK4_CREATESEQUENCEDICTIONARY (
        ch_sample_reference
    )
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions_gatk4.first())

    //
    // Build a per-sample reference bundle: [meta, fasta, fai, dict]
    // This ensures all downstream modules get the correct per-sample reference
    //
    ch_ref_bundle = ch_sample_reference
        .join(SAMTOOLS_FAIDX.out.fai, by: [0])
        .join(GATK4_CREATESEQUENCEDICTIONARY.out.dict, by: [0])
        // [meta, fasta, fai, dict]

    //
    // MODULE: Map reads to per-sample reference with BWA-MEM
    // Join reads with their sample's BWA index and reference FASTA
    //
    ch_bwa_all = ch_samplesheet
        .join(BWA_INDEX.out.index, by: [0])
        .join(ch_sample_reference, by: [0])
        .multiMap { meta, reads, index, fasta ->
            reads: [ meta, reads ]
            index: [ meta, index ]
            fasta: [ meta, fasta ]
        }

    BWA_MEM (
        ch_bwa_all.reads,
        ch_bwa_all.index,
        ch_bwa_all.fasta,
        true // sort_bam
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    //
    // MODULE: Convert FASTQ to unmapped BAM (needed for MergeBamAlignment)
    //
    PICARD_FASTQTOSAM (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions_picard.first())

    //
    // MODULE: Merge aligned and unmapped BAMs
    // Join aligned BAM + unmapped BAM + per-sample reference by meta
    //
    ch_merge_all = BWA_MEM.out.bam
        .join(PICARD_FASTQTOSAM.out.bam, by: [0])
        .join(ch_sample_reference, by: [0])
        .multiMap { meta, aligned_bam, unmapped_bam, fasta ->
            bams:  [ meta, aligned_bam, unmapped_bam ]
            fasta: [ meta, fasta ]
        }

    PICARD_MERGEBAMALIGNMENT (
        ch_merge_all.bams,
        ch_merge_all.fasta
    )
    ch_versions = ch_versions.mix(PICARD_MERGEBAMALIGNMENT.out.versions_picard.first())

    //
    // Conditional: Mark duplicates if requested
    //
    if (params.remove_duplicates) {
        ch_markdup_all = PICARD_MERGEBAMALIGNMENT.out.bam
            .join(ch_ref_bundle, by: [0])
            .multiMap { meta, bam, fasta, fai, dict ->
                bam:   [ meta, bam ]
                fasta: [ meta, fasta ]
                fai:   [ meta, fai ]
            }

        PICARD_MARKDUPLICATES (
            ch_markdup_all.bam,
            ch_markdup_all.fasta,
            ch_markdup_all.fai
        )
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions_picard.first())
        ch_bam_for_calling = PICARD_MARKDUPLICATES.out.bam
    } else {
        ch_bam_for_calling = PICARD_MERGEBAMALIGNMENT.out.bam
    }

    //
    // MODULE: Index BAMs
    //
    SAMTOOLS_INDEX (
        ch_bam_for_calling
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions_samtools.first())

    // Create BAM + BAI channel
    ch_bam_bai = ch_bam_for_calling
        .join(SAMTOOLS_INDEX.out.bai, by: [0])

    //
    // SUBWORKFLOW: BAM stats via samtools
    // Join BAM+BAI with per-sample reference
    //
    ch_bamstats_all = ch_bam_bai
        .join(ch_sample_reference, by: [0])
        .multiMap { meta, bam, bai, fasta ->
            bam_bai: [ meta, bam, bai ]
            fasta:   [ meta, fasta ]
        }

    BAM_STATS_SAMTOOLS (
        ch_bamstats_all.bam_bai,
        ch_bamstats_all.fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]})

    //
    // MODULE: Variant calling with GATK4 HaplotypeCaller
    // Note: --ploidy is injected from meta.ploidy via ext.args in modules.config
    // Join BAM+BAI with per-sample reference bundle
    //
    ch_hc_all = ch_bam_bai
        .join(ch_ref_bundle, by: [0])
        .multiMap { meta, bam, bai, fasta, fai, dict ->
            bam_bai: [ meta, bam, bai ]
            fasta:   [ meta, fasta ]
            fai:     [ meta, fai ]
            dict:    [ meta, dict ]
        }

    GATK4_HAPLOTYPECALLER (
        ch_hc_all.bam_bai,
        ch_hc_all.fasta,
        ch_hc_all.fai,
        ch_hc_all.dict
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions_gatk4.first())

    //
    // MODULE: Select SNP biallelic variants
    //
    ch_hc_vcf_tbi = GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, by: [0])

    SELECT_SNP_BIALLELIC (
        ch_hc_vcf_tbi
    )
    ch_versions = ch_versions.mix(SELECT_SNP_BIALLELIC.out.versions_gatk4.first())

    //
    // MODULE: Annotate variants with AlleleFraction
    // Join VCF+TBI with per-sample reference bundle
    //
    ch_va_all = SELECT_SNP_BIALLELIC.out.vcf
        .join(SELECT_SNP_BIALLELIC.out.tbi, by: [0])
        .join(ch_ref_bundle, by: [0])
        .multiMap { meta, vcf, tbi, fasta, fai, dict ->
            vcf_tbi: [ meta, vcf, tbi ]
            fasta:   [ meta, fasta ]
            fai:     [ meta, fai ]
            dict:    [ meta, dict ]
        }

    GATK4_VARIANTANNOTATOR (
        ch_va_all.vcf_tbi,
        ch_va_all.fasta,
        ch_va_all.fai,
        ch_va_all.dict
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTANNOTATOR.out.versions_gatk4.first())

    //
    // MODULE: Apply hard filters
    // Join VCF+TBI with per-sample reference bundle
    //
    ch_vf_all = GATK4_VARIANTANNOTATOR.out.vcf
        .join(GATK4_VARIANTANNOTATOR.out.tbi, by: [0])
        .join(ch_ref_bundle, by: [0])
        .multiMap { meta, vcf, tbi, fasta, fai, dict ->
            vcf_tbi: [ meta, vcf, tbi ]
            fasta:   [ meta, fasta ]
            fai:     [ meta, fai ]
            dict:    [ meta, dict ]
        }

    GATK4_VARIANTFILTRATION (
        ch_vf_all.vcf_tbi,
        ch_vf_all.fasta,
        ch_vf_all.fai,
        ch_vf_all.dict
    )
    ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions_gatk4.first())

    //
    // MODULE: Select filtered biallelic SNPs (exclude filtered)
    //
    ch_filt_vcf_tbi = GATK4_VARIANTFILTRATION.out.vcf
        .join(GATK4_VARIANTFILTRATION.out.tbi, by: [0])

    SELECT_FILTERED_BIALLELIC (
        ch_filt_vcf_tbi
    )
    ch_versions = ch_versions.mix(SELECT_FILTERED_BIALLELIC.out.versions_gatk4.first())

    //
    // MODULE: Generate IUPAC consensus FASTA
    // Join final VCF+TBI with per-sample reference bundle
    //
    ch_farm_all = SELECT_FILTERED_BIALLELIC.out.vcf
        .join(SELECT_FILTERED_BIALLELIC.out.tbi, by: [0])
        .join(ch_ref_bundle, by: [0])
        .multiMap { meta, vcf, tbi, fasta, fai, dict ->
            vcf_tbi: [ meta, vcf, tbi ]
            fasta:   [ meta, fasta ]
            fai:     [ meta, fai ]
            dict:    [ meta, dict ]
        }

    GATK4_FASTAALTERNATEREFERENCEMAKER (
        ch_farm_all.vcf_tbi,
        ch_farm_all.fasta,
        ch_farm_all.fai,
        ch_farm_all.dict
    )
    ch_versions = ch_versions.mix(GATK4_FASTAALTERNATEREFERENCEMAKER.out.versions_gatk4.first())

    //
    // ===== PER-LOCUS PARALLELIZATION =====
    //

    //
    // MODULE: Split BAM by reference contig (locus)
    // Use the processed BAM (after merge + optional dedup)
    //
    BAMTOOLS_SPLIT (
        ch_bam_for_calling
    )
    ch_versions = ch_versions.mix(BAMTOOLS_SPLIT.out.versions_bamtools.first())

    //
    // Channel operation: Flatten split BAMs into per-locus tuples
    // bamtools outputs files like sample.REF_locus1.bam, sample.REF_locus2.bam
    // We parse the filename to extract the locus name
    //
    ch_per_locus_bam = BAMTOOLS_SPLIT.out.bams
        .transpose()
        .map { meta, bam ->
            def locus = bam.name.replaceAll(/.*\.REF_/, '').replaceAll(/\.bam$/, '')
            def new_meta = meta + [locus: locus]
            [ new_meta, bam ]
        }

    //
    // MODULE: Index per-locus BAMs (required for H-PoPG random access)
    //
    SAMTOOLS_INDEX_PERLOCUS (
        ch_per_locus_bam
    )

    // Join per-locus BAM with its index
    ch_per_locus_bam_bai = ch_per_locus_bam
        .join(SAMTOOLS_INDEX_PERLOCUS.out.bai, by: [0])

    //
    // Channel operation: Pair per-locus BAMs with per-sample filtered VCF
    // Each sample has one filtered VCF but many per-locus BAMs.
    // combine(by:) creates one output per (sample, locus) since each sample has exactly one VCF.
    //
    ch_sample_vcf = SELECT_FILTERED_BIALLELIC.out.vcf
        .join(SELECT_FILTERED_BIALLELIC.out.tbi, by: [0])
        // [meta, vcf, tbi]

    ch_phase_input = ch_per_locus_bam_bai
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .combine(
            ch_sample_vcf.map { meta, vcf, tbi -> [ meta.id, vcf ] },
            by: [0]
        )
        .map { sample_id, meta, bam, bai, vcf -> [ meta, bam, bai, vcf ] }

    //
    // MODULE: H-PoPG haplotype phasing per sample x locus
    //
    HPOPG_PHASE (
        ch_phase_input
    )

    //
    // ===== ALLELE EXTRACTION =====
    //

    //
    // MODULE: Build phased consensus sequences per sample x locus
    // Each phase output needs: the filtered VCF, reference FASTA, and IUPAC FASTA
    // All three are per-sample, combined with per-locus phase output via sample ID
    //
    ch_consensus_input = HPOPG_PHASE.out.phase_out
        .map { meta, phase_out -> [ meta.id, meta, phase_out ] }
        .combine(
            ch_sample_vcf.map { meta, vcf, tbi -> [ meta.id, vcf ] },
            by: [0]
        )
        .combine(
            ch_sample_reference.map { meta, fasta -> [ meta.id, fasta ] },
            by: [0]
        )
        .combine(
            GATK4_FASTAALTERNATEREFERENCEMAKER.out.fasta.map { meta, fasta -> [ meta.id, fasta ] },
            by: [0]
        )
        .map { sample_id, meta, phase_out, vcf, ref_fasta, iupac_fasta ->
            [ meta, phase_out, vcf, ref_fasta, iupac_fasta ]
        }

    BUILD_PHASED_CONSENSUS (
        ch_consensus_input
    )
    ch_versions = ch_versions.mix(BUILD_PHASED_CONSENSUS.out.versions.first())

    //
    // MODULE: Collect phasing statistics across all samples
    //
    ch_all_phased = BUILD_PHASED_CONSENSUS.out.phased_fasta.collect { it[1] }
    ch_all_vcfs = ch_sample_vcf.map { meta, vcf, tbi -> vcf }.collect()
    ch_all_phases = HPOPG_PHASE.out.phase_out.collect { it[1] }
    ch_all_refs = ch_sample_reference.collect { it[1] }
    ch_ploidy_map = ch_samplesheet.map { meta, reads -> [id: meta.id, ploidy: meta.ploidy] }.collect()

    COLLECT_PHASING_STATS (
        ch_all_phased,
        ch_all_vcfs,
        ch_all_phases,
        ch_all_refs,
        ch_ploidy_map
    )
    ch_versions = ch_versions.mix(COLLECT_PHASING_STATS.out.versions.first())

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_versions_files = ch_versions.filter { it instanceof Path }
    softwareVersionsToYAML(ch_versions_files.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pate_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
