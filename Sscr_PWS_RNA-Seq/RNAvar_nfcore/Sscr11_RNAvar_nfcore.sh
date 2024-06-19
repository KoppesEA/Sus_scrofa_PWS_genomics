##did not use this script

nextflow run /ihome/crc/install/genomics_nextflow/nf-core-rnavar-1.0.0/workflow/ -r 1.0.0 \
	-name EAK_Sscr_RNAvar \
	-profile htc \
	-work-dir /bgfs/rnicholls/PWS_Sus_2017 \
	-resume \
	--input /bgfs/rnicholls/PWS_Sus_2017/RNAvar-nfcore/samplesheet.csv \
	--outdir /bgfs/rnicholls/PWS_Sus_2017/RNAvar-nfcore/ \
	--email eak37@pitt.edu \
	--multiqc_title EAK_Sscr_RNAvar_multiQC \
	--fasta /bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz \
	--gtf /bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/Sus_scrofa.Sscrofa11.1.104.gtf.gz \
	--gff /bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/Sus_scrofa.Sscrofa11.1.104.gff3.gz \
	--read_length 74 \
	--known_indels /bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/sus_scrofa.vcf.gz \
	--known_indels_tbi /bgfs/rnicholls/REFGenomes/Sscrofa_v11.1_v104/sus_scrofa.vcf.gz.tbi \
	--aligner star\
	--annotate_tools vep \
	
	/ihome/crc/install/genomics_nextflow/nf-core-cutandrun-3.0/workflow
