nextflow.enable.dsl=2

//---------------------------------------------------------------------//
// NANOSEQ NEXTFLOW PIPELINE - MAIN WORKFLOW (2023)                    //
// Author: Raul Alcantara (ra11), minor changes by A Baez-Ortega (ao7) //
// Pipeline parameters are hard-coded below                            //
//---------------------------------------------------------------------//

//---------------------------------------------------------------------//

def file_exists(x, name ) {
    if ( x == "" ) { return }
    assert file(x).exists() : "\n$name file $x was not found!\n\n"
}

include { CHECK_SIZE; DOWNLOAD; MERGE; SEQ_FILTER } from './modules/stageFiles.nf'
include { NANOSEQ } from './modules/NanoSeq_analysis.nf'
include { BWAMEM2_MAP; BWAMEM2_REMAP } from './modules/bwa.nf'
include { ADD_NANOSEQ_FASTQ_TAGS; MARKDUP; NANOSEQ_ADD_RB; NANOSEQ_DEDUP;
         VERIFY_BAMID; NANOSEQ_EFFI; NANOSEQ_VAF; FINALIZE; CLEANUP } from './modules/NanoSeq_aux.nf'

versions = Channel.empty() //accumulate all version files here


workflow {

  // Help message
  if (params.help) {
    println """
╔══════════════════════════════════════════════════════════════════════════════════════════════╗
║                              NANOSEQ NEXTFLOW PIPELINE - HELP                                ║
║                          Author: Raul Alcantara (ra11), A Baez-Ortega (ao7)                  ║
╚══════════════════════════════════════════════════════════════════════════════════════════════╝

DESCRIPTION:
    NanoSeq Nextflow pipeline for variant calling and analysis. The pipeline processes duplex 
    and normal BAM/CRAM files through quality control, preprocessing, and variant calling steps.

USAGE:
    nextflow run main.nf [OPTIONS]

REQUIRED PARAMETERS:
    --sample_sheet <file>       CSV file with 3 columns: id,duplex,normal
    --ref <file>                Reference FASTA file (must be named 'genome.fa', indexed)
    --post_triNuc <file>        Reference trinucleotide frequencies file

OPTIONAL INPUT FILES:
    --noise_bed <file>          SNP+NOISE BED mask file (bed.gz, indexed)
    --snp_bed <file>            SNP BED file (bed.gz, indexed) [auto-set for GRCh37/38]
    --panel <file>              Panel BED file for targeted analysis
    --part_excludeBED <file>    BED file with regions to exclude from partitioning

REFERENCE GENOME OPTIONS:
    --grch37                    Use predefined GRCh37 reference and annotations [false]
    --grch38                    Use predefined GRCh38 reference and annotations [false]

EXECUTION MODES:
    --short                     Short-run mode: only import, preprocessing & efficiency [false]
    --keep_temp                 Keep temporary work directory after completion [false]
    --remap                     Remap CRAMs using BWA-MEM2 [false]
    --location <string>         Data location: 'local' or 'irods' [irods]

PREPROCESSING OPTIONS:
    --seq <string>              Restriction site sequence for filtering [TGCA]
    --min_mapQ <int>            Minimum mapping quality for initial filtering [20]
    --cov_q <int>               Minimum mapping quality for duplex reads [0]
    --nanoseq_dedup_m <int>     Minimum reads in read bundle for deduplication [1]

COVERAGE ANALYSIS OPTIONS:
    --cov_exclude <string>      Comma-separated contigs to exclude [MT,MT%,hs37d5]
    --cov_include <string>      Comma-separated contigs to include (overrides exclude) ['']
    --cov_larger <int>          Only include contigs larger than this size [50000]

VARIANT CALLING PARAMETERS:
    --var_a <int>               Minimum AS-XS (alignment score difference) [50]
    --var_b <int>               Minimum matched normal reads per strand [0]
    --var_c <float>             Fraction of clips [0]
    --var_d <int>               Minimum duplex depth [2]
    --var_f <float>             Minimum fraction of reads for consensus [0.9]
    --var_i <float>             Maximum fraction of reads with indel [0.1]
    --var_m <int>               Minimum cycle number [8]
    --var_n <int>               Maximum number of mismatches [3]
    --var_p <float>             Minimum fraction of proper pairs [0]
    --var_q <int>               Minimum consensus base quality [60]
    --var_r <int>               Read length after 5' trimming [144]
    --var_v <float>             Maximum bulk VAF [0.01]
    --var_x <int>               Maximum cycle number [8]
    --var_z <int>               Minimum total normal reads [15]

INDEL CALLING PARAMETERS:
    --indel_rb <int>            Minimum reads in bundle [2]
    --indel_t3 <int>            Trim excess bases from 3' above this value [135]
    --indel_t5 <int>            Bases to trim from 5' [10]
    --indel_z <int>             Minimum total normal reads [15]
    --indel_v <float>           Maximum bulk VAF [0.01]
    --indel_a <int>             Minimum AS-XS for indels [50]
    --indel_c <float>           Fraction of clips for indels [0]

ANALYSIS PARAMETERS:
    --jobs <int>                Number of parallel jobs for NanoSeq [100]
    --dsa_d <int>               DSA minimum duplex depth [2]
    --dsa_q <int>               DSA minimum base quality for normal [30]
    --part_excludeCov <int>     Minimum coverage for partition exclusion [0]

VERIFYBAMID PARAMETERS:
    --vb_epsilon <string>       VerifyBAMid epsilon value [1e-12]
    --vb_ud <file>              VerifyBAMid UD file [auto-set for GRCh37/38]
    --vb_bed <file>             VerifyBAMid BED file [auto-set for GRCh37/38]
    --vb_mu <file>              VerifyBAMid MU file [auto-set for GRCh37/38]

OTHER OPTIONS:
    --study <string>            Study name for iRODS data access ['']
    --outDir <string>           Output directory [.]
    --help                      Show this help message and exit

SAMPLE SHEET FORMAT:
    The sample sheet must be a CSV file with headers: id,duplex,normal

    iRODS example:
    id,duplex,normal
    sample1,PD123a,PD123b
    sample2,PD456a,PD456b

    Local example:
    id,duplex,normal
    sample1,path/to/PD123a_duplex.bam,path/to/PD123b_normal.bam
    sample2,path/to/PD456a_duplex.bam,path/to/PD456b_normal.bam

    For iRODS data, prefix with study: study:sample_id
    Lines starting with # are ignored.

EXAMPLES:
    # Basic run with custom reference
    nextflow run main.nf --sample_sheet samples.csv --ref /path/to/genome.fa --post_triNuc trinuc.txt
    
    # Run with GRCh38 reference
    nextflow run main.nf --grch38 --sample_sheet samples.csv --post_triNuc trinuc.txt
    
    # Short run (preprocessing only)
    nextflow run main.nf --short --sample_sheet samples.csv --ref /path/to/genome.fa --post_triNuc trinuc.txt
    
    # Run with remapping
    nextflow run main.nf --remap --sample_sheet samples.csv --ref /path/to/genome.fa --post_triNuc trinuc.txt

For more information, see the pipeline documentation.
    """
    System.exit(0)
    }

    // checking parameters

    // input / output
    assert ( params.location == "local" || params.location == "irods" ) : "\nlocation parameter must be 'local' or 'irods'\n"
    assert ( params.ext == "cram" || params.ext == "bam" ) : "\next parameter must be 'cram' or 'bam'\n"
    // create named map with ext.file and ext.index, with cram and cram.crai if cram, bam and bam.bai if bam
    ext = params.ext == "cram" ? [file: "cram", index: "cram.crai", types: "uc"] : [file: "bam", index: "bam.bai", types: "m"]

    // remapping CRAMs
    assert ( params.remap == true || params.remap == false ) : "\nrealign parameter must be true or false\n"

    // use predefined parameter sets for GRCh37 & GRCh38
    assert ( params.ref != "" ) : "\nmust define a reference file genome.fa\n"
    assert ( params.ref.split("/")[-1] == "genome.fa" ) : "\nreference file must be named genome.fa\n"
    file_exists(params.ref, "ref")
    reference_path = params.ref.split("/")[0..-2].join('/')
    if ( params.remap ) {
        file_exists(reference_path + "/genome.fa.bwt.2bit.64", "bwa-mem2 index ")
        file_exists(reference_path + "/genome.fa.dict", "samtools dict ")
    }

    // check optional
    if ( params.panel == "" ) {
        panel_fh = file( 'NO_FILE_panel')
    } else {
        file_exists( params.panel, "panel BED" )
        panel_fh = file( params.panel )
    }
    if ( params.snp_bed != "" ) {
        file_exists( params.snp_bed, "snp_bed" )
        file_exists( params.snp_bed + ".tbi", "SNP BED index")
    }
    if ( params.noise_bed != "" ) {
        file_exists( params.noise_bed, "noise_bed" )
        file_exists( params.noise_bed + ".tbi", "noise BED index")
    }

    // part parameters
    file_exists(params.part_excludeBED, "part_excludeBED")
    file_exists(params.post_triNuc,"post_triNuc")

    // varify bam id params
    file_exists(params.vb_ud,"vb_ud")
    file_exists(params.vb_bed, "vb_bed")
    file_exists(params.vb_mu, "vb_mu")

    // process sample sheet
    file_exists(params.sample_sheet,"sample_sheet")

    list_ids = []
    samples = []
    sample_tags_a = []
    studies = []

    iline = 0
    new File(params.sample_sheet).splitEachLine(",") { fields -> 
        if ( iline == 0 ) {
            assert ( fields.contains("id")) : "\n Must specify a unique id column in the sample sheet\n"
            assert ( fields.contains("duplex")) : "\n Must specify a duplex column in the sample sheet\n"
            assert ( fields.contains("normal")) : "\n Must specify a normal column in the sample sheet\n"
            i_id = fields.findIndexOf{ it == "id"}
            i_duplex = fields.findIndexOf{ it == "duplex"}
            i_normal = fields.findIndexOf{ it == "normal"}
            iline = 1
        } else {
            row_id = fields[i_id]
            row_duplex_id = row_id
            row_normal_id = row_id

            if ( ! row_id.startsWith('#') ) {
                assert ( ! list_ids.contains(row_id) ) : "\nids in sample sheet must be unique\n\n"
                row_duplex = fields[i_duplex]
                row_normal = fields[i_normal]
                if ( row_duplex.contains(':') ) {
                    study_duplex = row_duplex.split(':')[0]
                    row_duplex_id = row_duplex.split(':')[1]
                } else {
                    assert ( params.study != "") : "\nMust define a study as an argument\n"
                    study_duplex = params.study
                }
                if ( row_normal.contains(':') ) {
                    study_normal = row_normal.split(':')[0]
                    row_normal_id = row_normal.split(':')[1]
                } else {
                    assert ( params.study != "" ) : "\nMust define a study as an argument\n"
                    study_normal = params.study
                }
                list_ids.add(row_id)
                tag_duplex = [ id: row_id, type : "duplex", name: row_duplex ]
                tag_normal = [ id: row_id, type : "normal", name: row_normal ]
                sample_tags_a.add([row_duplex, tag_duplex])
                sample_tags_a.add([row_normal, tag_normal])
                samples.add( [ [ study:study_duplex, id:row_duplex_id ], row_duplex] )
                samples.add( [ [ study:study_normal, id:row_normal_id ], row_normal] )
            }
        }
    }
    ch_samples = Channel.from(samples.unique())
    ch_samples.view()

    pout = params.toSorted{a,b -> a.key <=> b.key} // for ordered print
    if ( params.remap) {
        println("\nCarrying out remapping...\n")
    }
    println("\n")
    println("Run arguments: $pout")
    println("\n\n")
    println("Sample sheet contents:\n")
    ss = file( params.sample_sheet )
    println(ss.text)
    println("\n")

    // create sample files use the study name as the file name
    sam_file = ch_samples.map{ [it[0].study, it[1]] }.groupTuple().collectFile(cache: 'lenient', newLine: true) { 
        item ->[ "${item[0]}", item[1].join('\n') ]
    }
    
    if (params.location == "irods") {

      // *MODIFIED* (ao7): TEMPORARY - DISABLE CHECK_SIZE to avoid failure on farm22 
      // CHECK_SIZE( irods_data.getAbsolutePath(), sam_file, 1 )

      // download data from irods
      DOWNLOAD(params.outDir + "/irods_data", ch_samples, ext)
      versions = versions.concat(DOWNLOAD.out.versions.first())
      out_files = DOWNLOAD.out.files

    } else {
      
      // use local paths
      out_files = ch_samples
      
    }

    MERGE( params.outDir, out_files, params.remap ? 0 : params.min_mapQ, ext)
    versions = versions.concat(MERGE.out.versions.first())

    // *MODIFIED* (ao7): Commented out in order to place SEQ_FILTER after BWAMEM2_REMAP
    //Add metadata to samples
    //ch_cram_ss = MERGE.out.file.map{ [ [ name: it[0].id ], it[1], it[2] ] }

    if ( params.remap ) { //remapping

        MAP = BWAMEM2_REMAP( MERGE.out.file.map{ [ [ name: it[0].id ], it[1], it[2] ] }, reference_path, params.min_mapQ )

        MARKDUP( MAP.cram, reference_path )

        versions = versions.concat(MAP.versions.first())
        versions = versions.concat(MARKDUP.out.versions.first())

        ch_cram_ss = MARKDUP.out.cram
    
    } else {

        // *MODIFIED* (ao7): Moved line 309 to else clause
        ch_cram_ss = MERGE.out.file.map{ [ [ name: it[0].id ], it[1], it[2] ] }
        //
        
    }

    // *MODIFIED* (ao7): SEQ_FILTER step to filter restriction site motif
    SEQ_FILTER( ch_cram_ss, params.seq, reference_path )
    ch_cram = SEQ_FILTER.out.cram
    //

    NANOSEQ_ADD_RB( ch_cram , reference_path )

    versions = versions.concat(NANOSEQ_ADD_RB.out.versions.first())

    ch_tags = Channel.from(sample_tags_a)

    //* add labels for the pair analysis. This handles paired data using the same sample as duplex and normal
    ch_tag_n_add_rb = 
        NANOSEQ_ADD_RB.out.cram.map{[it[0].name,it[1],it[2]]}.cross(ch_tags).map{[it[1][1],it[0][1],it[0][2]]}

    NANOSEQ_DEDUP( NANOSEQ_ADD_RB.out.cram, reference_path, params.nanoseq_dedup_m )

    versions = versions.concat(NANOSEQ_DEDUP.out.versions.first())

    if ( params.vb_ud != "" &&  params.vb_bed != "" && params.vb_mu != "" ) {
        VERIFY_BAMID( NANOSEQ_DEDUP.out.cram, params.vb_epsilon, reference_path, params.vb_ud, params.vb_bed, params.vb_mu )

        versions = versions.concat(VERIFY_BAMID.out.versions.first())
    }
    
    //* add labels for the pair analysis. This handles paired data using the same sample as duplex and normal
    ch_tag_n_dedup = 
        NANOSEQ_DEDUP.out.cram.map{[it[0].name,it[1],it[2]]}.cross(ch_tags).map{[it[1][1],it[0][1],it[0][2]]}

    //collate channels to prepare input for efficiency calculation
    //must provide the CRAM from NANOSEQ_ADD_RB and NANOSEQ_DEDUP as arguments
    ch_add_rb_normal =ch_tag_n_add_rb.filter{ it[0].type == "normal"}.map{[it[0].id, it ]}
    ch_add_rb_duplex =ch_tag_n_add_rb.filter{ it[0].type == "duplex"}.map{[it[0].id, it ]}

    ch_dedup_normal = ch_tag_n_dedup.filter{ it[0].type == "normal"}.map{[it[0].id, it ]}
    ch_dedup_duplex = ch_tag_n_dedup.filter{ it[0].type == "duplex"}.map{[it[0].id, it ]}

    ch_normal_effi = ch_add_rb_normal.join( ch_dedup_normal ).map{it[1] + it[2][1..-1] }
    ch_duplex_effi = ch_add_rb_duplex.join( ch_dedup_duplex ).map{it[1] + it[2][1..-1] }

    ch_input_effi = ch_normal_effi.mix(ch_duplex_effi).unique{ it[0].name}.map{[[name: it[0].name]]+ it[1..-1]}

    NANOSEQ_EFFI( ch_input_effi, reference_path, panel_fh )

    versions = versions.concat(NANOSEQ_EFFI.out.versions.first())


    // *MODIFIED* (ao7): Added 'SHORT-RUN MODE' including only BAM import, preprocessing & efficiency calculation;
    //                   if params.short=true, the pipeline finishes here, otherwise it goes on to variant calling
    
    if ( params.short ) {
    
        println("\n(SHORT-RUN MODE selected: execution will finish after NANOSEQ_EFFI)\n")
    
    } else {

        //Collate CRAMs for NanoSeq call
        //Must provide :  NANOSEQ_ADD_RB.out.cram (duplex) & NANOSEQ_DEDUP.out.cram (normal) ; as input arguments

        ch_input_nanoseq = ch_add_rb_duplex.join(ch_dedup_normal).map{
            def meta = it[1][0].clone()
            meta.name = meta.id
            meta.type = "pair"
            [ meta ] + it[1][1..-1] + it[2][1..-1] }

        // *MODIFIED* (ao7): replaced params.min_mapQ (=20) with params.cov_q (=0)
        //NANOSEQ( ch_input_nanoseq, reference_path, params.jobs, params.min_mapQ, params.cov_exclude, params.cov_include,
        //    params.cov_larger, params.part_excludeBED, params.part_excludeCov, params.snp_bed, params.noise_bed, params.dsa_d, 
        //    params.dsa_q, params.var_a, params.var_b, params.var_c, params.var_d, params.var_f, params.var_i, params.var_m, 
        //    params.var_n, params.var_p, params.var_q, params.var_r, params.var_v, params.var_x, params.var_z, params.indel_rb,
        //    params.indel_t3, params.indel_t5, params.indel_z, params.indel_v, params.indel_a, params.indel_c, params.post_triNuc )
        NANOSEQ( ch_input_nanoseq, reference_path, params.jobs, params.cov_q, params.cov_exclude, params.cov_include,
            params.cov_larger, params.part_excludeBED, params.part_excludeCov, params.snp_bed, params.noise_bed, params.dsa_d, 
            params.dsa_q, params.var_a, params.var_b, params.var_c, params.var_d, params.var_f, params.var_i, params.var_m, 
            params.var_n, params.var_p, params.var_q, params.var_r, params.var_v, params.var_x, params.var_z, params.indel_rb,
            params.indel_t3, params.indel_t5, params.indel_z, params.indel_v, params.indel_a, params.indel_c, params.post_triNuc )
        //

        versions = versions.concat(NANOSEQ.out.versions.first())

        //NANOSEQ.out.results contains: muts.vcf, muts.vcf.tbi, indel.vcf, indel.vcf.tbi, cov.bed & cov.bed.tbi
        //VAF calculation also requires NANOSEQ_DEDUP.out.cram (duplex) in addition to the previous
        ch_input_vaf = NANOSEQ.out.results.map{[it[0].id, it ]}.join( ch_dedup_duplex ).map{it[1] + it[2][1..-1] }

        NANOSEQ_VAF( ch_input_vaf, reference_path )

        versions = versions.concat(NANOSEQ_VAF.out.versions.last())

        //combine the meta tags of duplex and normal
        ch_meta = ch_add_rb_duplex.join(ch_dedup_normal).map{[it[1][0],it[2][0]]}

        FINALIZE( versions.collect(), ch_meta )

        // *MODIFIED* (ao7): added new process for deletion of temp (work/) folder
        //                   (only if keep_temp == false)
        if ( ! params.keep_temp ) {
            versions = versions.concat(FINALIZE.out.versions.last())
            CLEANUP( versions.collect() )
        }

    }
}
