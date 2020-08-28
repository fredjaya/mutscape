#!/bin/bash

JS=/home/meep/Desktop/Biocomputing/jvarkit/dist/vcffilterjs.jar
DATA_PATH=/home/meep/Desktop/People/fred/Dropbox/meep/bee/01_data
VCF=$DATA_PATH/combined_nonVariants.vcf
OUT=$DATA_PATH/filteredCommons.vcf

# Run vcffilterjs to remove sites that share identical variants 
java -jar $JS -e 'function accept(v) {var g0= v.getGenotype(0);for(var i=1;i< v.getNSamples();i++) {if(!v.getGenotype(i).sameGenotype(g0)) return true;} return false;}accept(variant);' $VCF > $OUT
