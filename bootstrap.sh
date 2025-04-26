#!/bin/bash -i
##################
# VG is required
# KhufuEnv is required
##################
helpFunc()
{
   echo -e "Usage:         \033[46m$0 -RefworkDir RefworkDir -PanID PanID -seqfile seqfile -refGen refGen -t 2\033[0m
   \033[46mt\033[0m           NumOfThread
   \033[46m-/--gfa\033[0m     GFA file"
   exit 1
}
##################
t=4
gfa=""
##################
##################
SOPT='t:h'
LOPT=('gfa' 'seqfile')
OPTS=$(getopt -q -a --options ${SOPT} --longoptions "$(printf "%s:," "${LOPT[@]}")" --name "$(basename "$0")" -- "$@")
eval set -- $OPTS
while [[ $# > 0 ]]; do
    case ${1} in
		-h) helpFunc ;;
		-t) t=$2 && shift ;;
		--gfa) gfa=$2 && shift ;;
		esac
    shift
done
##################
if [[ $t == "" ]]; then t=2; fi
if [[ $gfa == "" ]]; then echo "gfa file should be provided" ; exit 0; fi
workDir=$(dirname $gfa)
prefix=$(echo $gfa | sed "s:.*/::g" | sed "s:.gfa.*$::g")
##################
echo "t=$t"
echo "gfa=$gfa"
##################
check=1
if (file $gfa | grep -q compressed ) ; then check=1 ; fi
if [ -f "$workDir"/"$prefix".gfa ]; then check=0;fi
if [[ "$check" == 1 ]]
then
   gunzip -kf $gfa
fi
##################
check=1
if ! [ -f "$workDir"/"$prefix".gbz ]; then check=0;fi
if ! [ -f "$workDir"/"$prefix".min ]; then check=0;fi
if ! [ -f "$workDir"/"$prefix".dist ]; then check=0;fi
if [[ "$check" == 0 ]]
then
   echo "gbz, min & dist will be produced"
   vg autoindex -g "$workDir"/"$prefix".gfa --workflow giraffe -t $t -p "$workDir"/"$prefix"
   cp "$workDir"/"$prefix".giraffe.gbz "$workDir"/"$prefix".gbz
else
   echo "available gbz, min, & dist will be used"
fi
du -sh "$workDir"/"$prefix".gbz | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
du -sh "$workDir"/"$prefix".min | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
du -sh "$workDir"/"$prefix".dist | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
#########
check=1
if ! [ -f "$workDir"/"$prefix".gbwt ]; then check=0;fi
if ! [ -f "$workDir"/"$prefix".ri ]; then check=0;fi
if ! [ -f "$workDir"/"$prefix".gg ]; then check=0;fi
if [[ "$check" == 0 ]]
then
   vg gbwt -p --num-threads "$t" -r "$workDir"/"$prefix".ri -g "$workDir"/"$prefix".gg -Z "$workDir"/"$prefix".gbz -o "$workDir"/"$prefix".gbwt
   echo "ri & gbwt will be produced"
else
   echo "available gbwt, ri and gg will be used"
fi
du -sh "$workDir"/"$prefix".gbwt | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
du -sh "$workDir"/"$prefix".ri | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
du -sh "$workDir"/"$prefix".gg | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
#########
if ! [ -f "$workDir"/"$prefix".hapl ]
then
   echo "hapl will be produced"
   vg haplotypes -v 2 -t"$t" -d "$workDir"/"$prefix".dist -r "$workDir"/"$prefix".ri -H "$workDir"/"$prefix".hapl "$workDir"/"$prefix".gbz
else
   echo "available hapl will be used"
fi
du -sh "$workDir"/"$prefix".hapl | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
#########
if ! [ -f "$workDir"/"$prefix".xg ]
then
   echo "xg will be produced"
   vg convert "$workDir"/"$prefix".gbz -x -t $t > "$workDir"/"$prefix".xg
else
   echo "available xg will be used"
fi
du -sh "$workDir"/"$prefix".xg | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
#########
if ! [ -f "$workDir"/"$prefix".snarls ]
then
   echo "snarls will be produced"
   vg snarls -T  "$workDir"/"$prefix".gbz > "$workDir"/"$prefix".snarls
else
   echo "available snarls will be used"
fi
du -sh "$workDir"/"$prefix".snarls | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
#########
if [ -f "$workDir"/"$prefix".vcf.gz ]
then
   gunzip -k "$workDir"/"$prefix".vcf.gz
   cp "$workDir"/"$prefix".vcf.gz.tbi "$workDir"/"$prefix".vcf.tbi
fi
##
if ! [ -f "$workDir"/"$prefix".vcf ]
then
   echo "vcf will be produced"
   vg deconstruct "$workDir"/"$prefix".xg -g "$workDir"/"$prefix".gbwt -a -t $t > "$workDir"/"$prefix".vcf
else
   echo "available vcf will be used"
fi
du -sh "$workDir"/"$prefix".vcf | sed "s:^:>>>>>>>>>\tKhufuPAN\t:g"
VCF="$workDir"/"$prefix".vcf
#########
vg stats -p 2 -zl "$workDir"/"$prefix".gfa > "$workDir"/"$prefix".gfa.stats 2> /dev/null &
cat "$workDir"/"$prefix".gfa | awk '{if($1=="S") print $0}'> "$workDir"/"$prefix".gfa.S.txt &
cat "$workDir"/"$prefix".gfa | awk '{if($1=="W") print $0}'> "$workDir"/"$prefix".gfa.W.txt &
cat "$workDir"/"$prefix".gfa| awk '{if($1=="L") print $0}'> "$workDir"/"$prefix".gfa.L.txt &
wait
#################
head=$(cat $VCF | sed -n  "/CHROM/ {p;n; p; q}"  | transverseDS | awk '{nA=split($2,A,"|"); if(nA==1) {print $1} else {for(a=1;a<=nA;++a) {print $1"."a}}  }' | transverseDS)
cat <(cat $VCF | sed -n '/##/,/#CHROM/{/#CHROM/{q};p}') <(echo $head | tr ' ' '\t') <(cat $VCF | grep -v "^#" | sed "s:^.*#::g" | tr '|' '\t' ) > "$VCF".phased
################
cat "$workDir"/"$prefix".gfa |grep -E "^S"  | cut -f 4| sort | uniq | sed "s:^.*#::g" | grep -Ev "^$" > "$workDir"/"$prefix".gfa_CHR.list
mkdir "$workDir"/CHRs
tcount=0
for chr in $(cat "$workDir"/"$prefix".gfa_CHR.list)
do
   echo $chr
   (
   mkdir "$workDir"/CHRs/"$chr"
   cat "$VCF".phased | grep -v "^##" | awk -v chr=$chr '{if(NR==1 || $1==chr) print $0}' > "$workDir"/CHRs/"$chr"/vcf0

   cat \
      <(cat "$workDir"/CHRs/"$chr"/vcf0 | sed -n '1,2p'| transverseDS | awk '{if($2~"[|]") { nA=split($2,A,"|"); S=$1"."1; for (a=2;a<=nA;++a){S=S"\n"$1"."a};  print S } else {print $1} }' | tr '\n' '\t' | sed "s:\t$:\n:g" ) \
      <(cat "$workDir"/CHRs/"$chr"/vcf0 | sed "1d" | tr '|' '\t' ) > "$workDir"/CHRs/"$chr"/vcf00

   cat "$workDir"/CHRs/"$chr"/vcf00 | awk 'OFS="\t"{if($6 >= 60) {print $0} }'  > "$workDir"/CHRs/"$chr"/vcf1
   cat "$workDir"/CHRs/"$chr"/vcf1 | cut -f 1,2,4,5 | awk '{print $1"_"$2"\t"$3","$4}' | paste - <(cat "$workDir"/CHRs/"$chr"/vcf1 | cut -f 10- | tr '\t' ',' ) > "$workDir"/CHRs/"$chr"/vcf2
   cat "$workDir"/CHRs/"$chr"/vcf2 | cut -f 3  | sed "s:\.::g" |  tr ',' '\t' | awk '{ delete a; for (i=1; i<=NF; i++) a[$i]++; n=asorti(a, b); for (i=1; i<=n; i++) printf b[i]" "; print "" }' | tr ' ' '\t'  |  awk '{if (NF > 1) {print 1} else {print 0} }'  | paste - "$workDir"/CHRs/"$chr"/vcf2 | awk '{if($1==1) print $0}' | cut -f 2- > "$workDir"/CHRs/"$chr"/vcf3
   cat "$workDir"/CHRs/"$chr"/vcf3 |  sed "s:[.],:-2,:g" | sed "s:[.]$:-2:g" | awk '{split($2,A,","); nB=split($3,B,"," ); for (b=1;b<=nB;++b) {C=C","A[B[b]+1]}; print C; C="" }' | sed "s:^,::g" | sed '1d' | sed "1i$(cat "$workDir"/CHRs/"$chr"/vcf3 | head -1 | cut -f 3)" | paste "$workDir"/CHRs/"$chr"/vcf3 -  | cut -f 1,4 > "$workDir"/CHRs/"$chr"/vcf4
   # filtering out alleles with N or X base
   cat "$workDir"/CHRs/"$chr"/vcf4  | cut -f 2 | sed '1d' | awk '{if($0 ~ "[NX]") {print 0} else {print 1} }'  | sed "1i1" | paste - "$workDir"/CHRs/"$chr"/vcf4 | awk '{if($1==1) print $0}' | cut -f 2- > "$workDir"/CHRs/"$chr"/vcf5
   cat "$workDir"/CHRs/"$chr"/vcf5 | cut -f 2 | sed '1d' | tr ',' '\t'  | awk '{ delete a; for (i=1; i<=NF; i++) a[$i]++; n=asorti(a, b); for (i=1; i<=n; i++) printf b[i]" "; print "" }' | tr ' ' '\t' | awk '{print NF";"$0}' | tr '\t' ',' |  sed "s:,$::g"  | tr ';' '\t' | awk '{print $2"\t"$1}' | sed "1iunique_alleles\tnumOfuniqueAlleles" | paste "$workDir"/CHRs/"$chr"/vcf5 - | tr ' ' '\t' > "$workDir"/CHRs/"$chr"/vcf6
   # get single call variants, only unque per locus, & sorting
   cat "$workDir"/CHRs/"$chr"/vcf6  | cut -f 1 | sed "1d" | tr '_' '\t' | sort -k1,1 -k2,2n | uniq -c | sed -E "s:^ +::g" | tr ' ' '\t'  | awk '{if($1==1) print $0}' | cut -f 2- | tr '\t' '_' | sed '1i#CHROM_POS' > "$workDir"/CHRs/"$chr"/vcf6.unique.list
   ###########
   merge "$workDir"/CHRs/"$chr"/vcf6.unique.list "$workDir"/CHRs/"$chr"/vcf6 > "$workDir"/CHRs/"$chr"/vcf6.unique
   cat "$workDir"/CHRs/"$chr"/vcf6.unique | sed '1d' | awk '{nA=split($2,A,","); nB=split($3,B,","); for(a=1;a<=nA;++a) { if (A[a] == "") {S=S","0} else {for(b=1;b<=nB;++b)  { if(A[a] == B[b] ) S=S","b }}} print S; S="" }' | sed "s:^[,]::g"  | sed "1i$(cat "$workDir"/CHRs/"$chr"/vcf6.unique | head -1 | cut -f2)" | paste "$workDir"/CHRs/"$chr"/vcf6.unique - | awk 'OFS="\t"{print $1,$2,$5,$3,$4}' | sed "s:#CHROM_POS:chr_pos:g" > "$workDir"/CHRs/"$chr"/vcf7
   cat "$workDir"/CHRs/"$chr"/vcf7 | sed '1d' | awk '{nA=split($4,A,","); for (a=1;a<=nA;++a) {S=S","length(A[a])}; print substr(S,2); S="" }' | sed '1ilen'| paste <(cat "$workDir"/CHRs/"$chr"/vcf7 | cut -f 1-4) - > "$workDir"/CHRs/"$chr"/vcf9
   cat "$workDir"/CHRs/"$chr"/vcf9 | sed '1d' | awk '{nA=split ($4,A,","); for(a=1;a<=nA;++a) {print ">"$1"_"a"\n"A[a]  } }'  > "$workDir"/CHRs/"$chr"/vcf9.fa
   cat "$workDir"/CHRs/"$chr"/vcf9 | sed '1d' | awk '{nA=split ($4,A,","); for(a=1;a<=nA;++a) {print $1"\t"$1"_"a"\t"A[a]  } }'  > "$workDir"/CHRs/"$chr"/vcf9.fa.txt
   ) &
   if [[ "$tcount" == $t ]]; then  wait; tcount=1; else tcount=$((tcount+1)); fi
done
wait
cat "$workDir"/CHRs/*/vcf9 | awk '{if(NR==1 || $1!~"chr_pos") print $0}' | (sed -u 1q; sort -k1,1V ) > "$workDir"/"$prefix".vcf9
cat "$workDir"/CHRs/*/vcf9.fa > "$workDir"/"$prefix".vcf9.fa
cat "$workDir"/CHRs/*/vcf9.fa.txt > "$workDir"/"$prefix".vcf9.fa.txt
#################
