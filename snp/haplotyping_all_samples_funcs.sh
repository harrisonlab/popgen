function convert2shapeit
{
    mapfile=$1
    genofile=$2 #assumed to be gzipped
    outdir=$3

    #rm -rf ${outdir}
    mkdir -p ${outdir}

    #split markers into one file per linkage group in map order ready for conversion to ped files
    #grid_run -Jsplit2lgs -M4 ${NODES}

    split_into_lgs.py ${mapfile} ${genofile} ${outdir}

    #grid_wait -Ljoblist

    #convert affy calls into ped format
    #rm -f joblist
    for fname in ${outdir}/*.tsv
    do
        outname=${fname/tsv/ped}
        #grid_run -Jaffy2ped ${NODES}
        affy2ped.py ${fname} > ${outname}
        #>> joblist

        pmname=${fname/tsv/pmap}
        recname=${fname/tsv/rec}
        #grid_run -Jmap2rec ${NODES}
        map2rec.py ${pmname} > ${recname}
        #>> joblist
    done
    #grid_wait -Ljoblist
}

function impute_haplotypes
{
    datadir=$1

    rm -f joblist
    for pedfile in ${datadir}/*.ped
    do
        pmfile=${pedfile/ped/pmap}
        recfile=${pedfile/ped/rec}
        phafile=${pedfile/ped/phased}
        grafile=${pedfile/ped/graph}
        logfile=${pedfile/ped/log}

        #impute haplotypes
        grid_run -Jshapeit ${MAXJOBS} ${NODES} -M4 -C4 \
            shapeit --input-ped  ${pedfile} ${pmfile} \
                    --input-map  ${recfile} \
                    --output-max ${phafile} \
                    --output-graph ${grafile} \
                    --output-log ${logfile} \
                    --thread 4 \
                    >> joblist
    done
    grid_wait -Ljoblist
}

function estimate_haploblocks
{
    datadir=$1

    for pedfile in ${datadir}/*.ped
    do
        prefix=${RANDOM}${RANDOM}
        pmfile=${pedfile/ped/pmap}
        blkfile=${pedfile/ped/blocks}

        #grid_run -Jplink -L14 ${NODES} -M4 -C4
        plink --ped  ${pedfile} --map ${pmfile} \
              --chr-set 95 no-xy \
              --blocks no-pheno-req no-small-max-span \
              --blocks-max-kb 1000 \
              --blocks-min-maf 0.03 \
              --out ${prefix}

        mv ${prefix}.blocks     ${blkfile}
        mv ${prefix}.blocks.det ${blkfile}.det
        mv ${prefix}.log        ${blkfile}.log

        #convert block boundaries into cM assuming 3.125 cM per megabase, 3.125e-6 cM per base
        echo $(basename ${pedfile} .ped) > ${blkfile}.cm
        tail -n +2 ${blkfile}.det \
            | awk '{print $2*3.125e-6} END{print $3*3.125e-6}' \
            >> ${blkfile}.cm
    done
    #grid_wait -Ljoblist
}

function extract_haplotypes2ped
{
    datadir=$1

    rm -f joblist

    for pedfile in ${datadir}/*.ped
    do
        extract_haplotypes_to_ped.py ${pedfile} ${pedfile/ped/phased.haps} ${pedfile/ped/ped.phased}
        extract_haplotypes_to_ped2.py ${pedfile} ${pedfile/ped/phased.haps} ${pedfile/ped/ped.phased2}
        cat ${pedfile/ped/pmap} | awk '{print $2,$4}' > ${pedfile/ped/info}
    done

    #haploview
}
