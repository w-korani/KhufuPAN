#!/bin/bash -i
##################
# VG is required
# KhufuEnv is required
##################
helpFunc()
{
   echo -e "
   Usage:   \033[46m$0 -gfa ref/asd.gfa -l sample.lists\033[0m
   \033[46m-t\033[0m       NumOfThread 
   \033[46m-o\033[0m       output 
   \033[46m-gfa\033[0m     gfa input file
   \033[46m-l\033[0m       ListOfSamples[Optional: all samples in the gams will processed if not used] 
   \033[46m-/--gams\033[0m gams
   \033[46m-/--min\033[0m  MinimumDepth[default=2]
   \033[46m-/--max\033[0m  MinimumDepth[default=2]
   "
   exit 1
}
##################
t=2
out=""
list=""
min_dep=2
max_dep=20
##################
SOPT='t:o:l:h'
LOPT=('gams' 'min' 'max' 'gfa')
OPTS=$(getopt -q -a --options ${SOPT} --longoptions "$(printf "%s:," "${LOPT[@]}")" --name "$(basename "$0")" -- "$@")
eval set -- $OPTS
while [[ $# > 0 ]]; do
    case ${1} in
		-h) helpFunc ;;
		-t) t=$2 && shift ;;
		-o) output=$2 && shift ;;
		--gams) gams=$2 && shift ;;
		-l) list=$2 && shift ;;
		--min) min_dep=$2 && shift ;;
      --max) max_dep=$2 && shift ;;
      --gfa) gfa=$2 && shift ;;
    esac
    shift
done
##################
tmpDir=$(mktemp -d "KhufuEnviron.XXXXXXXXX")
###
if [[ -z $gfa ]]; then echo "gfa should be provided"; exit 0 ;fi
workDir=$(dirname $gfa)
prefix=$(echo $gfa | sed "s:.*/::g" | sed "s:.gfa.*$::g")
###
if [[ ! -d $gams ]]; then echo "gams should be provided"; exit 0 ;fi
###
if [[ $list == "" ]]
then
   ls -l "$gams"/*.vcf5 | sed "s:.*/::g" | sed "s:.vcf5::g" > "$tmpDir"/geno.list
   list="$gams"
else
   ls -l "$bams"/*.vcf5 | sed "s:.*/::g" | sed "s:.vcf5::g" | grep -wf "$list" > "$tmpDir"/geno.list
fi
###
if [[ -z "$output" ]]
then
   output="khufuPAN_03_"$list.panmap
fi
if [[ $(echo $output | awk '{if($0 ~ ".panmap$" ) {print 1} else {print 0} }') == 0 ]]; then echo output suffix should be .panmap; exit 0 ;fi
#####################
echo "t=$t"
echo "gfa=$gfa"
echo "list=$list"
echo "min_dep=$min_dep"
echo "max_dep=$max_dep"
echo "gams=$gams"
echo "output=$output"
#####################
tcount=1
for id in $(cat "$tmpDir"/geno.list )
do
   (
   echo "$id"
   cat $gams/"$id".vcf5  | sed "1d" | tr ';' '\t' | awk -v min_dep=$min_dep  -v max_dep=$max_dep '{if($2 < min_dep || $2 > max_dep ) {print "-\t0"} else {print $0} }'  > "$tmpDir"/"$id".txt 
   cat "$tmpDir"/"$id".txt | cut -f 2 | sed "1i$id" > "$tmpDir"/"$id".dep
   cat "$tmpDir"/"$id".txt | cut -f 1 | sed "1i$id" | awk 'function sort(str){ nA=split(str,A,","); asort(A, B); S=B[1] ; for(a=2;a<=nA;++a){ S=S","B[a] } return S } {print sort($0) }' > "$tmpDir"/"$id".call
   ) &
   if [[ "$tcount" == $t ]]; then  wait; tcount=1; else tcount=$((tcount+1)); fi
done
wait
##
paste <(cat "$workDir"/"$prefix".vcf9| awk '{print $1"\t"$5"\t"$3}'  | tr '_' '\t') "$tmpDir"/*.call > "$output" &
paste <(cat "$workDir"/"$prefix".vcf9 | awk '{print $1"\t"$5"\t"$3}'  | tr '_' '\t') "$tmpDir"/*.dep > "$output".dep &
wait
cp "$workDir"/"$prefix".vcf9.fa "$output".fa
##################
rm -rf $tmpDir
trap "rm -rf $tmpDir" EXIT
##################