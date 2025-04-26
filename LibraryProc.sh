#!/bin/bash -i
##################
# VG is required
# KhufuEnv is required
# KMC is required
##################
helpFunc()
{
   echo -e "Usage:       \033[46m$0 -gfa ref/asd.gfa -gams gams -id asd\033[0m
   \033[46m-t\033[0m        NumOfThread
   \033[46m-/--id\033[0m    SampleID
   \033[46m-/--gfa\033[0m   GFA graph file (fill path)
   \033[46m-/--gams\033[0m  gams (working dirctory)
   \033[46m-/--r1\033[0m    FastqPair1
   \033[46m-/--r2\033[0m    FastqPair2 [optional]
   \033[46m-/--mem\033[0m   memory to be used by kmer
   \033[46m-/--kmer\033[0m  kmer length, max 31, default is 25, kmer=0 stops the personalized graph.
   \033[46m-/--clean\033[0m flag to clean up gams folder after run"
   exit 1
}
##################
t=2
gams="gams"
clean=0
mem=12
kmer=25
clean=1
##################
##################
SOPT='t:q:h'
LOPT=('gams' 'id' 'gfa' 'r1' 'r2' 'mem' 'kmer' 'clean')
OPTS=$(getopt -q -a --options ${SOPT} --longoptions "$(printf "%s:," "${LOPT[@]}")" --name "$(basename "$0")" -- "$@")
eval set -- $OPTS
while [[ $# > 0 ]]; do
    case ${1} in
		-h) helpFunc ;;
		-t) t=$2 && shift ;;
		--id) id=$2 && shift ;;
		--gams) gams=$2 && shift ;;
		--gfa) gfa=$2 && shift ;;
		--r1) fq1=$2 && shift ;;
		--r2) fq2=$2 && shift ;;
      --mem) mem=$2 && shift ;;
      --kmer) kmer=$2 && shift ;;
      --clean) clean=$2 && shift ;;
		esac
    shift
done
##################
#echo $OPTS
echo "t=$t"
echo "id=$id"
echo "gams=$gams"
echo "gfa=$gfa"
echo "fq1=$fq1"
echo "fq2=$fq2"
echo "mem=$mem"GB
echo "kmer=$kmer"
if [[ $clean == 1 ]]; then echo "gams folder will be cleaned after run" ; elif [[ $clean == 0 ]]; then echo "" ; else echo "clean should be 1 or 0" ; fi
##################
if [[ $gfa == "" ]]; then echo "gfa should be provided" ; exit 0; fi
if [[ $id == "" ]]; then echo "sample ID should be provided" ; exit 0; fi
if [[ $fq1 == "" ]]; then echo "fastq1 should be provided" ; exit 0; fi
if [[ $fq2 == "" ]]; then echo "single-read end processing" ; fi
##################
#################
workDir=$(dirname $gfa)
prefix=$(echo $gfa | sed "s:.*/::g" | sed "s:.gfa.*$::g")
################
mkdir -p $gams
mkdir -p stds
###
dep=1; dep=$(cat "$workDir"/"$prefix".gfa.stats | awk '{if($1=="length") print $2}')
(
zcat $fq1 $fq2 | awk '{if(NR%4==2) print $0}' | tr -d '\n' | wc -c | sed "s:$: mapping to $dep:g" > "$gams"/"$id".len
cat "$gams"/"$id".len | cut -d' ' -f 1 | awk -v dep=$dep '{printf("%.2f\n", $0/(dep*1))}' | sed "s:^:$id\t:g"  > "$gams"/"$id".dep
sed -i "s:^:$id\t:g" "$gams"/"$id".len
) &
###
mkdir "$gams"/TMP_"$id"
fqR1=$fq1
if [[ -n $fq2 ]]; then fqR2=$fq2; fi
##########
#### mapping
if [[ -n $fq2  ]]
then
   echo "paired-end will be mapped"
   if [[ "$kmer" -eq 0  ]]
   then
      echo " >>> not personalized graph >>> "
      vg giraffe -p -t "$t" -Z $(echo $gfa | sed "s:.gfa$:.gbz:g") -d $(echo $gfa | sed "s:.gfa$:.dist:g") -m $(echo $gfa | sed "s:.gfa$:.min:g") -f "$fqR1" -f "$fqR2" > "$gams"/TMP_"$id"/"$id".gam
   else
      echo " >>> personalized graph >>> "
      ## merging/sorting fastq
      cat "$fqR1" "$fqR2" > "$gams"/TMP_"$id"/12.fq.gz
      ## get kff file
      kmc -k"$kmer" -okff -m"$mem" -t"$t" -hp "$gams"/TMP_"$id"/12.fq.gz "$gams"/TMP_"$id"/"$id" "$gams"/TMP_"$id"
      ## subgraphing
      vg haplotypes -v 2 -t "$t" --include-reference --diploid-sampling  -i "$workDir"/"$prefix".hapl -k "$gams"/TMP_"$id"/"$id".kff -g "$gams"/TMP_"$id"/"$id".gbz "$workDir"/"$prefix".gbz
      ## mapping
      vg giraffe -p -t "$t" -Z "$gams"/TMP_"$id"/"$id".gbz -f "$fqR1" -f "$fqR2" > "$gams"/TMP_"$id"/"$id".gam
   fi
else
   echo "single-end will be mapped"
   if [[ "$kmer" -eq 0  ]]
   then
      echo " >>> not personalized graph >>> "
      vg giraffe -p -t "$t" -Z $(echo $gfa | sed "s:.gfa$:.gbz:g") -d $(echo $gfa | sed "s:.gfa$:.dist:g") -m $(echo $gfa | sed "s:.gfa$:.min:g") -f "$fqR1" > "$gams"/TMP_"$id"/"$id".gam
   else
      echo " >>> personalized graph >>> "
      ## get kff file
      kmc -k"$kmer" -okff -m"$mem" -t"$t" -hp "$fqR1" "$gams"/TMP_"$id"/"$id" "$gams"/TMP_"$id"
      ## subgraphing
      vg haplotypes -v 2 -t "$t" --include-reference --diploid-sampling  -i "$workDir"/"$prefix".hapl -k "$gams"/TMP_"$id"/"$id".kff -g "$gams"/TMP_"$id"/"$id".gbz "$workDir"/"$prefix".gbz
      ## mapping
      vg giraffe -p -t "$t" -Z "$gams"/TMP_"$id"/"$id".gbz -f "$fqR1" > "$gams"/TMP_"$id"/"$id".gam
   fi
fi
##########
vg stats -a "$gams"/TMP_"$id"/"$id".gam > "$gams"/"$id".gam.stat
##################
#filtering
vg filter -t $t --min-mapq 60 "$gams"/TMP_"$id"/"$id".gam > "$gams"/TMP_"$id"/"$id".mq60.gam
vg stats -a "$gams"/TMP_"$id"/"$id".mq60.gam > "$gams"/"$id".mq60.gam.stat
##################
# packing & calling
vg pack -x "$workDir"/"$prefix".xg -g "$gams"/TMP_"$id"/"$id".gam -o "$gams"/TMP_"$id"/"$id".gam.pack -Q 5
vg call -a "$workDir"/"$prefix".xg -k "$gams"/TMP_"$id"/"$id".gam.pack -s "$id" -t "$t" -g "$workDir"/"$prefix".gbwt -r  "$workDir"/"$prefix".snarls > "$gams"/"$id".vcf
##################
# processing
indx=$(cat "$gams"/"$id".vcf | grep -v "^##" | head  -1 | tr '\t' '\n' | awk -v id=$id '{if($0==id) print NR }')
cat "$gams"/"$id".vcf | grep -v "^#" | cut -f 1,2,4,5,"$indx" | grep -v "\./\." |  tr ':' '\t'  | cut -f 1-7 | awk '{print $1"_"$2"\t"$3","$4"\t"$5"\t"$6"\t"$7 }' | sed "1ichr_pos\talleles\tGT\tDP\tAD" > "$gams"/"$id".vcf2
merge "$workDir"/"$prefix".vcf9  "$gams"/"$id".vcf2 | cut -f 4,6-  >  "$gams"/"$id".vcf3
cat "$gams"/"$id".vcf3  | cut -f 1,2,5 | sed 1d | awk '{nA=split($1,A,",");nB=split($2,B,",");split($3,C,","); S="-"; for(a=1;a<=nA;++a){s=0; for(b=1;b<=nB;++b){ if(A[a]==B[b]){s=C[b]} }; S=S","s  }; gsub("^-,","",S) ;  print$1"\t"S  }' > "$gams"/"$id".vcf4
cat "$gams"/"$id".vcf4  | awk '{nA=split($1,A,",");split($2,B,","); s="-"; dep=0; for(a=1;a<=nA;++a){if(B[a]>0){ s=s","a ;dep=dep+B[a]} }; gsub("^-,","",s); print s";"dep }' | sed "1i$id" > "$gams"/"$id".vcf5
cat "$gams"/"$id".vcf5 | sed "1d" | cut -d";" -f 2| sort -n | uniq -c | sed -E "s:^ +::g" | tr ' ' '\t' |  awk '{printf "%05d\t%s\n", $1,$2}' | sed "s:^:$(cat "$gams"/"$id".dep)\t:g" | sed "1i#id\tOverallDep\tcnt\tcoverage" > "$gams"/"$id".vcf5.stat
#################
# cleaning
if [[ $clean == 1 ]]; then
	echo "the following files will be deleted:"
	ls "$gams"/TMP_"$id"/*
	rm -r "$gams"/TMP_"$id"
fi
#################
echo ">>>>>>>>>   KhufuPAN: sample "$id" was processed"
#################