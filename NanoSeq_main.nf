nextflow.enable.dsl=2

//---------------------------------------------------------------------//
// NANOSEQ NEXTFLOW PIPELINE - MAIN WORKFLOW (2023)                    //
// Author: Raul Alcantara (ra11), minor changes by A Baez-Ortega (ao7) //
// Pipeline parameters are hard-coded below                            //
//---------------------------------------------------------------------//

// SHORT-RUN MODE: if SHORT=true, only sample import and preprocessing are done //
// (USE COMMAND-LINE OPTION: --short)
SHORT = false

// KEEP TEMP FILES: if KEEP_TMP=true, the `work` directory is not deleted after success //
// (USE COMMAND-LINE OPTION: --keep_temp)
KEEP_TEMP = false

// INPUT PATHS //
// SAMPLES: CSV file with 3 columns headed: id,duplex,normal        (--sample_sheet)
// REF:     Reference FASTA (must be named 'genome.fa', indexed)    (--ref)
// TRINUC:  Reference trinucleotide frequencies                     (--post_triNuc)
// SNP_NOISE: SNP+NOISE BED mask (bed.gz, indexed)                  (--noise_bed)
SAMPLES   = ""   // "sample_ids.csv"
REF       = ""   // "/lustre/scratch125/casm/team268im/ao7/nanoseq_pipe/SNPmask_pipeline/ref/mouse/genome.fa"                 
TRINUC    = ""   // "/lustre/scratch125/casm/team268im/ao7/nanoseq_pipe/SNPmask_pipeline/ref/mouse/trinuc_freq.txt"  
SNP_NOISE = ""   // "/lustre/scratch125/casm/team268im/ao7/nanoseq_pipe/SNPmask_pipeline/output/masks/SNP+NOISE.NV2.SeenNV1.mouse.bed.gz"

// OPTIONS FOR PREPROCESSING //
REMAP = false  // Remap CRAMs?                                      (--remap)
SEQ = "TGCA"   // Restriction site sequence (for filtering)         (--seq)

// OPTIONS FOR NANOSEQ COVERAGE ANALYSIS (runNanoSeq.py cov)
//EXCLUDE = "MT,MT%,GL%,JH%,NC_%,hs37d5" // Comma-separated list of contigs to exclude       (--cov_exclude)
EXCLUDE = "MT,MT%,hs37d5"               // Comma-separated list of contigs to exclude       (--cov_exclude)
LARGER = 50000                          // Only include contigs larger than this size (0)   (--cov_larger)

// OPTIONS FOR NANOSEQ VARIANT CALLING (runNanoSeq.py var) //
A = 50    // minimum AS-XS (diff between primary and secondary alignment scores) (50)       (--var_a)
B = 0     // min matched normal reads per strand (0 for nanoseq normal, 5 for WGS normal)   (--var_b)
C = 0     // fraction of clips (0.02)                                                       (--var_c)
D = 2     // minimum duplex depth (2)                                                       (--var_d)
F = 0.9   // minimum fraction of reads for consensus (0.9)                                  (--var_f)
I = 0.1   // maximum fraction of reads with an indel (1.0)                                  (--var_i)
M = 8     // minimum cycle number (8)                                                       (--var_m)
N = 3     // maximum number of mismatches (3)                                               (--var_n)
P = 0     // minimum fraction of reads that are proper-pairs (0)                            (--var_p)
Q = 60    // minimum consensus base quality (60)                                            (--var_q)
R = 144   // read length (after 5' trimming) (144)                                          (--var_r)
V = 0.01  // maximum bulk VAF (0.01)                                                        (--var_v)
X = 8     // maximum cycle number (8)                                                       (--var_x)
Z = 15    // minimum total number of normal reads (15)                                      (--var_z)

// OPTIONS FOR NANOSEQ INDEL CALLING (runNanoSeq.py indel) //
RB = 2    // minimum reads in a bundle (2)                                                  (--indel_rb)
T3 = 135  // excess bases above this value are trimmed from 3' (135)                        (--indel_t3)
T5 = 10   // bases to trim from 5' reads (10)                                               (--indel_t5)

//---------------------------------------------------------------------//

//Setting up parameters

// *MODIFIED* (ao7) - added short-run mode and keep-temp mode parameters
params.short = SHORT
params.keep_temp = KEEP_TEMP
//

params.bwa_image = "docker://quay.io/wtsicgp/pcap-core:5.7.0"
//params.nanoseq_image = "docker://quay.io/wtsicgp/nanoseq:3.5.5" // USING LOCAL VERSION: /nfs/dog_n_devil/adrian/software/NanoSeq_3.5.7/bin
//params.nanoseq_image = "docker://quay.io/wtsicgp/nanoseq:3.3.0" // previous version

// *** remapping CRAMs
params.remap = REMAP
assert ( params.remap == true || params.remap == false ) : "\nrealign parameter must be true or false\n"

//*use predefined parameter sets for GRCh37 & GRCh38
params.grch37 = false
assert ( params.grch37 == true || params.grch37 == false ) : "\ngrch37 parameter must be true or false\n"
params.grch38 = false
assert ( params.grch38 == true || params.grch38 == false ) : "\ngrch38 parameter must be true or false\n"

if ( params.grch37 ) {
    params.ref = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/genome.fa"
} else if ( params.grch38 ) {
    params.ref = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa"
} else {
    params.ref = REF
}
assert ( params.ref != "" ) : "\nmust define a reference file genome.fa\n"
assert ( params.ref.split("/")[-1] == "genome.fa" ) : "\nreference file must be named genome.fa\n"
file_exists(params.ref,"ref")
reference_path = params.ref.split("/")[0..-2].join('/')
if ( params.remap ) {
    file_exists(reference_path + "/genome.fa.bwt.2bit.64", "bwa-mem2 index ")
    file_exists(reference_path + "/genome.fa.dict", "samtools dict ")
}

// *MODIFIED* (ao7): new default directory
//params.outDir = projectDir
params.outDir = "."
//
// *** Preprocessing and mapping params
params.fastq_tags_m = 3
params.fastq_tags_s = 4
// *MODIFIED* (ao7) - corrected to m=1
//params.nanoseq_dedup_m = 2
params.nanoseq_dedup_m = 1  // minimum number of reads in RB (hard-coded)
//

// *** NanoSeq parameters
params.jobs = 100
params.panel = ""

// *MODIFIED* (ao7): restriction site sequence
params.seq = SEQ
//

// mapping and cov
params.min_mapQ = 20
// *MODIFIED* (ao7): separate min_mapQ and cov_q into separate parameters
params.cov_q = 0    // minimum mapQ to include a duplex read (hard-coded)
//
params.cov_larger = LARGER
params.cov_include = ""
if ( params.grch37 ) { //for GRCh37 reference
    params.cov_exclude = "MT,GL%,NC_%,hs37d5"
    params.snp_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/SNP.sorted.bed.gz"
    params.noise_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/botseq/NOISE.sorted.bed.gz"
} else if ( params.grch38 ) {//for GRCh38 reference
    params.cov_exclude = "chrM,chr%_random,chrUn_%,chr%_alt,HLA-%" 
    params.snp_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/botseq/SNP.sorted.GRCh38.bed.gz"
    params.noise_bed = "/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/botseq/NOISE.sorted.GRCh38.bed.gz"
} else {
    params.cov_exclude = EXCLUDE
    params.noise_bed = SNP_NOISE
    params.snp_bed = ""
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
params.part_excludeBED = ""
file_exists( params.part_excludeBED, "part_excludeBED" )
params.part_excludeCov = 0

// dsa parameters
params.dsa_d = D
params.dsa_q = 30            // minimum base quality for normal (hard-coded)
params.dsa_M = params.cov_q  // minimum mapping quality (hard-coded)

// variantcaller parameters
params.var_a = A
params.var_b = B
params.var_c = C
params.var_d = D
params.var_f = F
params.var_i = I
params.var_m = M
params.var_n = N
params.var_p = P
params.var_q = Q
params.var_r = R
params.var_v = V
params.var_x = X
params.var_z = Z
// indel parameters
params.indel_rb = RB
params.indel_t3 = T3
params.indel_t5 = T5
params.indel_z = params.var_z
params.indel_v = params.var_v
params.indel_a = params.var_a
params.indel_c = params.var_c
// post paramaters
params.post_triNuc = TRINUC
file_exists(params.post_triNuc,"post_triNuc")

// ** VerifyBAMid params
params.vb_epsilon ="1e-12"
if ( params.grch37 ) { //GRCh37
    params.vb_ud ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.UD"
    params.vb_bed ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.bed"
    params.vb_mu ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCH37d5/verifybamid/ALL_500K.strictmasked.ok.vcf.mu"
} else if ( params.grch38 ) { //GRCh38
    params.vb_ud ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.UD"
    params.vb_bed ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.bed"
    params.vb_mu ="/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/verifybamid/1000g.phase3.100k.b38.vcf.gz.dat.mu"
} else {
    params.vb_ud = ""
    params.vb_bed = ""
    params.vb_mu = ""
}
file_exists(params.vb_ud,"vb_ud")
file_exists(params.vb_bed, "vb_bed")
file_exists(params.vb_mu, "vb_mu")

// *** getting data from IRODs
params.study = ""

download_dir = params.outDir.toString() + "/irods_data" // download of raw data
File irods_data = new File( download_dir)
if ( ! irods_data.exists() ) {
    irods_data.mkdir()
}

// *** process sample sheet
params.sample_sheet = SAMPLES
file_exists(params.sample_sheet,"sample_sheet")

list_ids = []
samples = []
sample_tags_a = []
studies = []

iline = 0
new File(params.sample_sheet).splitEachLine(",") {fields -> 
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
        // *MODIFIED* (ao7): skip commented lines (#)
        if ( ! row_id.startsWith('#') ) {
        //
            assert ( ! list_ids.contains(row_id) ) : "\nids in sample sheet must be unique\n\n"
            row_duplex = fields[i_duplex]
            row_normal = fields[i_normal]
            if ( row_duplex.contains(':') ) {
                study_duplex = row_duplex.split(':')[0]
                row_duplex = row_duplex.split(':')[1]
            } else {
                assert ( params.study != "" ) : "\nMust define a study as an argument\n"
                study_duplex = params.study
            } 
            if ( row_normal.contains(':') ) {
                study_normal = row_normal.split(':')[0]
                row_normal = row_normal.split(':')[1]
            } else {
                assert ( params.study != "" ) : "\nMust define a study as an argument\n"
                study_normal = params.study
            } 
            list_ids.add(row_id)
            tag_duplex = [ id: row_id, type : "duplex", name: row_duplex ]
            tag_normal = [ id: row_id, type : "normal", name: row_normal ]
            sample_tags_a.add([row_duplex, tag_duplex ])
            sample_tags_a.add([row_normal, tag_normal ])
            samples.add( [ [study:study_duplex, id:row_duplex], row_duplex] )
            samples.add( [ [study:study_normal, id:row_normal], row_normal] )
        }
    }
}
ch_samples = Channel.from(samples.unique())

def file_exists(x, name ) {
    if ( x == "" ) { return }
    assert file(x).exists() : "\n$name file $x was not found!\n\n"
}

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

include { CHECK_SIZE; DOWNLOAD; MERGE; SEQ_FILTER } from './modules/stageFiles.nf'
include { NANOSEQ } from './modules/NanoSeq_analysis.nf'
include { BWAMEM2_MAP; BWAMEM2_REMAP } from './modules/bwa.nf'
include { ADD_NANOSEQ_FASTQ_TAGS; MARKDUP; NANOSEQ_ADD_RB; NANOSEQ_DEDUP;
         VERIFY_BAMID; NANOSEQ_EFFI; NANOSEQ_VAF; FINALIZE; CLEANUP } from './modules/NanoSeq_aux.nf'

versions = Channel.empty() //accumulate all version files here


workflow {
    //create sample files use the study name as the file name
    sam_file = ch_samples.map{ [it[0].study, it[1]] }.groupTuple().collectFile(cache: 'lenient', newLine: true) { 
        item ->[ "${item[0]}", item[1].join('\n') ]
    }
    // *MODIFIED* (ao7): TEMPORARY - DISABLE CHECK_SIZE to avoid failure on farm22 
    //CHECK_SIZE( irods_data.getAbsolutePath(), sam_file, 1 )
    // [line below was already commented out]
    //
    //CHECK_SIZE.out.view()
    DOWNLOAD( irods_data.getAbsolutePath(), ch_samples )
    versions = versions.concat(DOWNLOAD.out.versions.first())

    MERGE( params.outDir, DOWNLOAD.out.files, params.remap ? 0 : params.min_mapQ )
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
