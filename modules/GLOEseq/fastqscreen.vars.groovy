FASTQSCREEN_OUTDIR=QC + "/fastqscreen"
FASTQSCREEN_THREADS=Integer.toString(ESSENTIAL_THREADS)
//fastqscreen additional param e.g. subset or bowtie /bowtie 2 parameters
FASTQSCREEN_PARAM="--aligner 'bowtie2' --nohits --subset 100000"
//the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
//this could be e.g.
FASTQSCREEN_CONF="human::" + ESSENTIAL_BOWTIE_Human + "," + "yeast::" + ESSENTIAL_BOWTIE_Yeast + "," + "bovine::" + ESSENTIAL_BOWTIE_Bovine + "," + "bacteria::" + ESSENTIAL_BOWTIE_Bacteria  + "," + "mycoplasma::" + ESSENTIAL_BOWTIE_Mycoplasma + "," + "phix::" + ESSENTIAL_BOWTIE_PhiX +  "," + "ERCC::" + ESSENTIAL_BOWTIE_ERCC +  ","  + "rRNA::" + ESSENTIAL_BOWTIE_rRNA