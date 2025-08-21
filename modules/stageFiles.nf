nextflow.enable.dsl=2

/* export types can be comma separated. Valid values are:
#                       a - all data (default)
#                       m - mapped data (sample bam file and unmapped pairs files)
#                       r - results (raw, flagged and archives for results)
#                       s - statistics (mapping stats, pass/fail states)
#                       u - un-aligned bam files
#                     Alternatively, sub-set of results can be exported:
#                       cn - copynumber results
#                       re - rearrangement results
#                       su - substitution results
#                       in - indel results
#                       rna - rna results
#                       snp6 - snp6 results (a subset of cn)
#                       ms - methylation results
#                       ds - duplex results (botseq, nanoseq)
#                       tb - telbam
#                       jsc - jlShearwater counts results
#                       gps - GRIPSS results
#                       gds - GRIDSS results
*/
// *MODIFIED* (ao7): '-noconvert_cram_to_bam' replaced by '-t uc' in newer dataImportExport versions
//TYPE = "-t u"; EXT=".cram";EXTI=".cram.crai";ODIR="raw_lane_bams"
TYPE = "-t uc"; EXT=".cram";EXTI=".cram.crai";ODIR="raw_lane_bams"
//
//TYPE = "-t m"; EXT=".bam";ODIR="mapped_sample"
//TYPE = "-t jsc";EXT=".bed.gz";ODIR="jlSCounts"
//TYPE = "-t gps";EXT=".bam";ODIR="GRIDSS"
//TYPE = "-t gds";EXT=".vcf.gz";ODIR="GRIPSS"
//TYPE = "-t in";EXT=".vcf.gz";ODIR="pindel"
//TYPE = "-t tb";EXT=".bam";ODIR="telbam"
//TYPE = "-t su";EXT=".vcf.gz";ODIR="cavemanC"

//check to see if we have enough space to get all the files
process CHECK_SIZE {

    input :
        path dest_out
        path sample_file
        val study

    output :
        stdout

    errorStrategy 'terminate'
    executor = 'local'

    // *MODIFIED* (ao7): change to match farm22 module name
    //beforeScript 'module purge; module add dataImportExport'
    //beforeScript 'module purge; module add dataImportExport/1.56.0'
    beforeScript 'module purge; module add dataImportExport/1.60.4'
    //

    script:
        def opts1 = task.ext.args ?: '-noconvert_cram_to_bam'
        def opts2 = study == 1 ? "-study " + sample_file.getBaseName() : "-project " + sample_file.getBaseName()
        """
        touch ${task.process}
        NS1=`wc -l $sample_file | sed 's/ .*\$//'`

        ## *MODIFIED* (ao7): '-noconvert_cram_to_bam' ($opts1) replaced by '-t uc' in newer dataImportExport versions
        ##exportData.pl -n -f -st -l live $opts1 -s $sample_file -o $dest_out $opts2 $TYPE | grep -v bsub | tee log.txt
        exportData.pl -n -f -st -l live -s $sample_file -o $dest_out $opts2 $TYPE | grep -v bsub | tee log.txt
        ##

        NS2=`grep 'bam file sizes' log.txt | wc -l`

        if [ \$NS1 -ne \$NS2 ]; then
            echo "\nSample was not found in IRODs! Found \$NS2 out of \$NS1 samples"
            exit 1
        fi

        """

    stub:
        """
        echo 'done'
        """
}

process DOWNLOAD {

    tag "$sample"

    input :
        val dest_out
        tuple val(meta), val(sample)
    
    publishDir "$dest_out/${meta.containsKey("study") ? meta.study : meta.project}/$sample/$ODIR", mode: 'link', pattern : "*$EXT", overwrite: true

    output :
        tuple val(meta), path("*$EXT"), emit: files
        path ("versions.yml"), emit: versions

    cpus 1
    memory {  task.exitStatus == 130  ? 2.GB * task.attempt : 2.GB }

    maxForks 10 // DO NOT CHANGE
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' } //delayed retry
    maxRetries 5


    // *MODIFIED* (ao7): change to match farm22 module name
    //beforeScript 'module purge; module add dataImportExport'
    //beforeScript 'module purge; module add dataImportExport/1.56.0'
    beforeScript 'module purge; module add dataImportExport/1.60.4'
    //

    script :
        def study = meta.containsKey("study") ? meta.study : meta.project
        def opts1 = task.ext.args ?: '-noconvert_cram_to_bam'
        def opts2 = meta.containsKey("study") ? "-study $meta.study" : ''
        def opts3 = meta.containsKey("project") ? "-project $meta.project" : ''
        """
        touch ${task.process}_$sample

        ## *MODIFIED* (ao7): '-noconvert_cram_to_bam' ($opts1) replaced by '-t uc' in newer dataImportExport versions
        ##exportData.pl -n -st -l live $opts1 -s $sample -o $dest_out $opts2 $opts3 $TYPE
        exportData.pl -n -st -l live -s $sample -o $dest_out $opts2 $opts3 $TYPE
        ##
        mv $dest_out/$study/$sample/*/*$EXT .
        #mv $dest_out/${study}/$sample/*/*$EXTI .

        set +u #otherwise perl print causes error
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dataImportExport: \$echo \$(perl -MSanger::CGP::DataImportExport::Export::ExportData  -e 'print \$Sanger::CGP::DataImportExport::Export::ExportData::VERSION;')
        END_VERSIONS
        """

    stub:
        """
        sleep \$[ ( \$RANDOM % 20 )  + 1 ]s
        touch ${sample}$EXT

        set +u #otherwise perl print causes error
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dataImportExport: \$echo \$(perl -MSanger::CGP::DataImportExport::Export::ExportData  -e 'print \$Sanger::CGP::DataImportExport::Export::ExportData::VERSION;')
        END_VERSIONS
        """
}

process MERGE {

    tag "$meta.id"

    input :
        val dest_out
        tuple val(meta), path(files)
        val min_mapQ
    
    //publishDir "${dest_out}", mode: 'link', pattern: "merged/*", overwrite: true

    output :
        tuple val(meta),path("merged/${meta.id}$EXT"), path("merged/${meta.id}${EXTI}"), emit: file
        path ("versions.yml"), emit: versions

    cpus 3
    memory {  task.exitStatus == 130  ? 2.GB * task.attempt : 2.GB }

    maxRetries 3

    // *MODIFIED* (ao7): change to match farm22 module name
    //beforeScript 'module purge; module add samtools'
    //beforeScript 'module purge; module add samtools-1.19.2'
    //

    script:
        """
        set -o pipefail
        touch ${task.process}_${meta.id}
        mkdir -p merged
        samtools merge -@ $task.cpus -u - $files | \
            samtools view -@ $task.cpus -q $min_mapQ -C -o merged/${meta.id}$EXT
        
        samtools index -@ $task.cpus merged/${meta.id}$EXT

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        touch MERGE
        mkdir -p merged
        touch merged/${meta.id}${EXT}
        touch merged/${meta.id}${EXTI}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

// *MODIFIED* (ao7): Added process SEQ_FILTER
// Filters CRAM reads containing restriction site motif (params.seq)
process SEQ_FILTER {

    tag "${meta.name}"

    input:
        tuple val(meta), path("in.cram"), path("in.cram.crai")
        val seq
        path ref_path

    output:
        tuple val(meta), path("${meta.name}.cram"), path("${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    cpus 4
    memory { task.exitStatus == 137 ? 6.GB * task.attempt : 6.GB }

    script:
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
 
        samtools view -@ $task.cpus -T ${ref_path}/genome.fa -h in.cram | awk '\$10 !~ /$seq/ || \$0 ~/^@/' | \\
            samtools view -@ $task.cpus -C -T ${ref_path}/genome.fa  > ${meta.name}.cram
        samtools index ${meta.name}.cram

        cat <<-END_VERSIONS2 > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS2
        """

    stub:
        """
        touch ${meta.name}.cram
        touch ${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}
//

/* sample workflow
workflow {
    ch_samples = Channel.of("MD6258ab_ds0001","MD6258ab_ds0002")
    sam_file = ch_samples.collectFile(name: 'samples.txt', newLine: true)
    size_check = CHECK_SIZE("/lustre/scratch117/casm/team78/ra11/temp/data",sam_file,"","2508")
    //size_check.view()
    DOWNLOAD("/lustre/scratch117/casm/team78/ra11/temp/data",ch_samples,"","2508", size_check.collect() )
    MERGE("/lustre/scratch117/casm/team78/ra11/temp/merged",DOWNLOAD.out.files,"","2508",0)
    SEQ_FILTER(MERGE.out.file,"TGCA",reference_path)
}
*/
