nextflow.enable.dsl=2

if (! params.containsKey("outDir")) {
    params.outDir = projectDir
}

outDir = params.outDir

//Run this for fastq files prior to mapping
process ADD_NANOSEQ_FASTQ_TAGS {
    
    tag "${meta.name}"

    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(reads)
        val m
        val s

    output:
        tuple val(meta), path("out/*.fastq.gz"), emit: fastqs
        path  "versions.yml", emit: versions

    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { ( task.exitStatus == 130 ) ? 500.MB * task.attempt : 500.MB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 500.MB * task.attempt : 500.MB }
    //

    script:
        def read1 = reads[0]
        def read2 = reads[1]
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
        mkdir -p out
        L=`zcat $read1 | head -2 | tail -1 | awk '{ print length }'`
        extract_tags.py -a $read1 -b $read2 -c ./out/${meta.name}_R1.fastq.gz -d ./out/${meta.name}_R2.fastq.gz -m $m -s $s -l \$L
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            extract-tags.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """

    stub:
        def read1 = reads[0]
        def read2 = reads[1]
        """
        mkdir -p out
        touch ./out/${meta.name}_R1.fastq.gz
        touch ./out/${meta.name}_R2.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            extract-tags.py: \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

process MARKDUP {

    tag "${meta.name}"

    container params.bwa_image

    input:
        tuple val(meta), path(cram), path(index)
        path(ref_path)

    output:
        tuple val(meta), path("optdup/${meta.name}.cram"), path("optdup/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    cpus 3
    errorStrategy { sleep(Math.pow(2, task.attempt) * 10000 as long); return 'retry' } //delayed retry
    // *MODIFIED* (ao7): increased memory; added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    // memory { ( task.exitStatus == 130 || task.exitStatus == 140) ? 25.GB * task.attempt : 25.GB }
    // queue { task.exitStatus == 140 ? "long" : "normal" }
    memory { ( task.exitStatus == 130 || task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? 40.GB * task.attempt : 40.GB }
    queue { ( task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? "long" : "normal" }
    //

    script:
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
        mkdir -p nsorted
        mkdir -p optdup
        ln -sf ../$cram ./nsorted/${meta.name}.cram
        samtools view -H $cram | grep SO:queryname > /dev/null || \\
            ( rm ./nsorted/${meta.name}.cram; samtools sort -@ $task.cpus -O cram -m 2G -n -o ./nsorted/${meta.name}.cram $cram )
        samtools view -T ${ref_path}/genome.fa -@ $task.cpus --no-PG -h ./nsorted/${meta.name}.cram | bamsormadup inputformat=sam level=0 blocksortverbose=0 rcsupport=1 threads=$task.cpus fragmergepar=$task.cpus optminpixeldif=10 | \\
            bammarkduplicatesopt verbose=0 level=0 index=0 optminpixeldif=2500 | samtools view -T ${ref_path}/genome.fa --no-PG -@ $task.cpus -C -o ./optdup/${meta.name}.cram
        samtools index -@ $task.cpus ./optdup/${meta.name}.cram
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamsormadup: \$(bamsormadup -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            bammarkduplicatesopt: \$(bammarkduplicatesopt -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    stub:
        """
        mkdir -p nsorted
        mkdir -p optdup
        ln -sf ../$cram ./nsorted/${meta.name}.cram
        touch ./optdup/${meta.name}.cram
        touch ./optdup/${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamsormadup: \$(bamsormadup -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            bammarkduplicatesopt: \$(bammarkduplicatesopt -v 2>&1 | head -1 | sed 's/.*version //' | sed 's/\\.\$// ')
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process NANOSEQ_ADD_RB {

    tag "${meta.name}"
    
    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        path( ref_path )

    //publishDir "${outDir}", mode: 'link', pattern : "add_rb/*", overwrite: true

    output:
        tuple val(meta), path("add_rb/${meta.name}.cram"), path("add_rb/${meta.name}.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 2.GB * task.attempt : 2.GB }
    //queue { task.exitStatus == 140 ? "long" : "normal" }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 2.GB * task.attempt : 2.GB }
    queue { ( task.exitStatus == 140 || task.exitStatus == Integer.MAX_VALUE ) ? "long" : "normal" }
    //

    script:
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
        mkdir -p add_rb

        #determine if this a NanoSeq or WGS library
        NLINES=`samtools view $cram | head -10000 | grep rb: | grep rc: | grep mb: | grep mc: | wc -l` || true
        ## *MODIFIED* (ao7): fixed numeric comparison operator
        ## if [ \$NLINES != 0 ]; then
        if [ \$NLINES -ne 0 ]; then
        ##
            bamaddreadbundles -I $cram -O ./add_rb/${meta.name}.cram
            samtools index ./add_rb/${meta.name}.cram
        cat <<-END_VERSIONS1 > versions.yml
        "${task.process}":
            bamaddreadbundles: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS1
        else
            ln -sf ../$cram ./add_rb/${meta.name}.cram
            ln -sf ../$index ./add_rb/${meta.name}.cram.crai
        cat <<-END_VERSIONS2 > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS2
        fi
        """
    stub:
        """
        mkdir -p add_rb
        touch ./add_rb/${meta.name}.cram
        touch ./add_rb/${meta.name}.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bamaddreadbundles: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process NANOSEQ_DEDUP {

    tag "${meta.name}"
    
    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        path (ref_path)
        val m

    //publishDir "${outDir}", mode: 'link', pattern : "dedup/*", overwrite: true

    output:
        tuple val(meta), path("dedup/${meta.name}.neat.cram"), path("dedup/${meta.name}.neat.cram.crai"), emit: cram
        path  "versions.yml", emit: versions

    cpus 1
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 10.GB * task.attempt : 10.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 10.GB * task.attempt : 10.GB }
    //
    
    script:
        """
        set -o pipefail
        touch ${task.process}_${meta.name}
        mkdir -p dedup

        NLINES1=`samtools view -H $cram | grep ^@PG | grep ID:bamaddreadbundles | wc -l` || true
        NLINES2=`samtools view -H $cram | grep ^@PG | grep ID:randomreadinbundle | wc -l` || true
        ## *MODIFIED* (ao7): fixed numeric comparison operators
        ## if [ \$NLINES1 > 0 ] && [ \$NLINES2 == 0 ]; then
        if [ \$NLINES1 -gt 0 ] && [ \$NLINES2 -eq 0 ]; then
        ##
            ## *MODIFIED* (ao7): added -m parameter
            ##randomreadinbundle -I $cram -O ./dedup/${meta.name}.neat.cram
            randomreadinbundle -I $cram -O ./dedup/${meta.name}.neat.cram -m $m
            ##
            samtools index ./dedup/${meta.name}.neat.cram
        else
            ln -sf ../$cram ./dedup/${meta.name}.neat.cram
            ln -sf ../$index ./dedup/${meta.name}.neat.cram.crai
        fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            randomreadinbundle: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    stub:
        """
        mkdir -p dedup
        touch ./dedup/${meta.name}.neat.cram
        touch ./dedup/${meta.name}.neat.cram.crai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            randomreadinbundle: \$(runNanoSeq.py -v)
            samtools: \$(echo \$(samtools --version 2>&1) | head -1 | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

process VERIFY_BAMID {

    tag "${meta.name}"
    
    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index)
        val epsilon
        path ref_path
        path vb_ud
        path vb_bed
        path vb_mu

    publishDir "$params.outDir/outNextflow/QC/$meta.name", mode: 'copy', pattern: "verifyBAMid/*", overwrite: true

    output:
        path "verifyBAMid/${meta.name}.verifyBAMid.txt", emit: verifybamid
        path  "versions.yml", emit: versions

    maxRetries 1
    cpus 2
    memory 6.GB

    script:
        """
        touch ${task.process}_${meta.name}
        mkdir -p verifyBAMid
        #processing of CRAM is super slow so need to convert to BAM
        samtools view -b -@ $task.cpus  $cram > temp.bam 
        samtools index -@ $task.cpus temp.bam
        ( VerifyBamID --Epsilon $epsilon --UDPath $vb_ud \\
            --BedPath $vb_bed --MeanPath $vb_mu \\
            --Reference ${ref_path}/genome.fa --BamFile temp.bam > verifyBAMid/${meta.name}.verifyBAMid.txt 2> verifyBAMid/error ) || ( cp verifyBAMid/error verifyBAMid/${meta.name}.verifyBAMid.txt )
        rm temp.bam
        rm temp.bam.bai
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            VerifyBamID: \$(VerifyBamID -v 2>&1 | grep Version | sed 's/.*://')
        END_VERSIONS
        """
    stub:
        """
        mkdir -p verifyBAMid
        touch verifyBAMid/${meta.name}.verifyBAMid.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            VerifyBamID: \$(VerifyBamID -v 2>&1 | grep Version | sed 's/.*://')
        END_VERSIONS
        """
}

// *MODIFIED* (ao7): Increased number of retries
//MAXN = 1
MAXN = 3
//
process NANOSEQ_EFFI {
    
    tag "${meta.name}"

    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(cram), path(index), path(cram_neat), path(index_neat)
        path ref_path
        file panel

    publishDir "$params.outDir/outNextflow/QC/$meta.name", mode: 'copy', pattern : "effi/*" , overwrite: true

    output:
        path "effi/${meta.name}.effi.tsv", emit: effi
        path "versions.yml", emit: versions
        // *MODIFIED* (ao7): output coverage
        path "effi/${meta.name}.meancov.tsv", emit: meancov
        //


    maxRetries MAXN
    cpus 4
    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130  ? 25.GB * task.attempt : 25.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 25.GB * task.attempt : 25.GB }
    //

    errorStrategy { task.attempt == MAXN ? 'ignore' : 'retry' }

    script:
        def cpus = task.cpus - 1
        def panel_arg = panel.name != 'NO_FILE_panel' ? "-p $panel" : ''
        """
        touch ${task.process}_${meta.name}
        mkdir -p effi
        efficiency_nanoseq.pl -t $cpus -d $cram_neat -x $cram -o effi/${meta.name}.effi -r ${ref_path}/genome.fa $panel_arg

        # *MODIFIED* (ao7): calculate mean dedup coverage across sites with cov >1x
        cov=`samtools depth $cram_neat | awk '\$3 > 1 { sum += \$3; n++ } END { print sum / n; }'`
        echo -e "${meta.name}\\t\$cov" > effi/${meta.name}.meancov.tsv
        #
 
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            efficiency_nanoseq.pl : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
    stub:
        """
        mkdir -p effi
        touch effi/${meta.name}.effi.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            efficiency_nanoseq.pl : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

process NANOSEQ_VAF {

    tag "${meta.id}_${meta.type}"
    
    // *MODIFIED* (ao7): disabled Docker container; using local install on farm22
    // container params.nanoseq_image

    input:
        tuple val(meta), path(vcf_muts), path(index_muts), path(vcf_indel), path(index_indel),
             path(bed_cov), path(index_cov), path(cram_neat), path(index_neat)
        path ref_path

    publishDir "$params.outDir/outNextflow/NanoSeq/$meta.id", mode: 'copy', pattern: "*.vcf.*", overwrite: true

    output:
        tuple val(meta),path("${meta.id}.vcf.gz"), path("${meta.id}.vcf.gz.tbi"), emit: vcf
        path  "versions.yml", emit: versions

    cpus 2

    // *MODIFIED* (ao7): added task.exitStatus == Integer.MAX_VALUE;
    //                   the default value of task.exitStatus in the absence of an .exitcode file for the previous execution
    //memory { task.exitStatus == 130 ? 10.GB * task.attempt : 10.GB }
    memory { ( task.exitStatus == 130 || task.exitStatus == Integer.MAX_VALUE ) ? 10.GB * task.attempt : 10.GB }
    //

    script:
        """
        touch ${task.process}_${meta.id}_${meta.type}
        mkdir -p out
        #snv_merge uses bam2R which blows up memory when using CRAM
        samtools view -@ $task.cpus -b $cram_neat > tempBam.bam
        samtools index -@ $task.cpus tempBam.bam 
        snv_merge_and_vaf_calc.R $vcf_muts $vcf_indel tempBam.bam $bed_cov out/${meta.id}.vcf
        rm tempBam.bam
        bcftools sort -Oz out/${meta.id}.vcf -o ${meta.id}.vcf.gz
        bcftools index -t ${meta.id}.vcf.gz
        rm out/${meta.id}.vcf
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snv_merge_and_vaf_calc.R : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
    stub:
        """
        mkdir out
        mkdir sorted
        touch ${meta.id}.vcf.gz
        touch ${meta.id}.vcf.gz.tbi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snv_merge_and_vaf_calc.R : \$(runNanoSeq.py -v)
        END_VERSIONS
        """
}

process FINALIZE {

    input:
        val versions
        tuple val(meta_d), val(meta_n)

    output:
        // *MODIFIED* (ao7): emit versions to allow linking with CLEANUP process
        path  "versions.yml", emit: versions
        //path "versions.yml"
        //
        path "QC_*"

    publishDir "$params.outDir/outNextflow/NanoSeq/$meta_d.id", mode: 'copy', overwrite: true

    executor 'local'
    cpus 1
    memory 500.MB

    script:
        def allversions = versions.join(' ')
        """
        cat $allversions > versions.yml
        ln -sf $params.outDir/outNextflow/QC/$meta_d.name QC_duplex
        ln -sf $params.outDir/outNextflow/QC/$meta_n.name QC_normal
        """
    stub:
        def allversions = versions.join(' ')
        """
        cat $allversions > versions.yml
        ln -sf $params.outDir/outNextflow/QC/$meta_d.name QC_duplex
        ln -sf $params.outDir/outNextflow/QC/$meta_n.name QC_normal
        """
}

// *MODIFIED* (ao7): new process for deletion of `work` folder after success
//                   (only if keep_temp == false)
process CLEANUP {
    input:
        val versions
    executor 'local'
    cpus 1
    memory 500.MB
    script:
        """
        echo "keep_temp = ${params.keep_temp}"
        echo "Executing: find ${workDir} -type f,l ! -name '.command.*' -delete; find ${workDir}/../output/irods_data -type f,l -delete"
        find ${workDir} -type f,l ! -name '.command.*' -delete
        find ${workDir}/../output/irods_data -type f,l -delete
        #echo "Executing: rm -rf ${workDir}/*/*/* ${workDir}/../output/irods_data/*"
        #rm -rf ${workDir}/*/*/* ${workDir}/../output/irods_data/*  ## error: too many arguments
        """
    stub:
        """
        echo "keep_temp = ${params.keep_temp}"
        """
}
