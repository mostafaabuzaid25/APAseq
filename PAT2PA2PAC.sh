#!/bin/bash -l

#-------------------------------------------------
## Ex.
## dos2unix PAT2PA2PAC.sh
## bash PAT2PA2PAC.sh "P_heterocycla.fa" "/xx/" "/xx/yy" 24 "1.sam,2.sam" "s1,s2"

## Requirements:
## MAP_parseSAM2PAT.pl PAT_setIP_big.pl bedtools bedmap awk sed

# Parameters
#1 chrfa="P_heterocycla.fa" # Chromosome sequence, used for setIP
#2 indir=/xx/ # Working directory
#3 odir=/xx/yy/ # Output result directory
#4 dist=24 # Distance between adjacent PA and PAC
#5 samfiles="1.sam,2.sam" # SAM files to be processed, note: SAM files should already be filtered based on different mapping tools' results and converted to PAT format using perl
#6 lbls="test1,test2" # Corresponding labels

## Output [odir] Note: start coordinates are 0-based or 1-based
#all.PAC.info           -- Info about PACs
#all.PAC.header         -- Header row of PAC matrix [chr UPA_start(1-base) UPA_end strand PAnum tot_tagnum coord(1-base) refPAnum]/s1...sN
#all.PAC.PAcount        -- Number of PAs for each sample under PAC [info]/s1...sN (sum of samples should be less than or equal to PAnum)
#all.PAC.PATcount       -- Number of PATs for each sample under PAC [info]/s1...sN (sum of samples should be equal to tot_tagnum)
#all.PA.uniq.bed             -- Coordinates of all PAs across samples (deduplicated, with total PAT count) chr/start(0-base)/end/./PATnum/strand
#----------- The following are intermediate files for individual samples -----------
# PAT01.txt.PA
#PAT01.txt.PA.bed
#PAT01.txt.PA.bed.PAinPAC.bed
#PAT01.txt.PA.bed.PATinPAC.bed
#PAT01.txt.PA.bed.real
#PAT01.txt.PAT
#PAT02.txt.PA
#PAT02.txt.PA.bed
#PAT02.txt.PA.bed.PAinPAC.bed
#PAT02.txt.PA.bed.PATinPAC.bed
#PAT02.txt.PA.bed.real
#PAT02.txt.PAT

#chr	UPA_start	UPA_end	strand	PAnum	tot_tagnum	coord	refPAnum	s1	s2
#chr1	99	103	+	5	25	100	12	5	5
#chr1	1099	1103	+	5	12	1100	6	0	5

#-------------------------------------------------


chrfa=$1
indir=$2
odir=$3
dist=$4
#samfiles=$5
#lbls=$5
temp=$5

## DEBUG Testing Pay attention to the DEBUG lines in the code ##
#chrfa="P_heterocycla.fa"
#indir="/data/gpfs01/szhou/WXH"
#odir="/data/gpfs01/szhou/WXH/test1/"
#dist=24
#samfiles="PAT01.txt,PAT02.txt"
#lbls="s1,s2"


# Printing parameter values
printf "chrfa=$chrfa\nindir=$indir\nodir=$odir\ndist=$dist\nsamfiles=$temp\nlbls=$temp\n"

## Parameter settings
cd $indir

# Loop through each file in 'temp'
for file in $temp
do
 echo "${file}"
 samfiles+=($file)
 lbls+=($file)
done

#echo "print $samfiles"

# The following lines seem to be attempts to modify 'samfiles' and 'lbls' arrays
# They are commented out in the original script

##samfiles=(${samfiles//[, ;]/ });
#lbls=(${lbls//[, ;]/ });

##samfiles=$(find $indir -maxdepth 1 -name "*.bam*" -type f)

## Checking if the chromosome fasta file exists
if [ ! -f $chrfa ]; then
    echo "$chrfa not found!"
    exit 1
fi

# Loop through each file in 'samfiles' array
for file in "${samfiles[@]}";
do
  echo "********check $file exist************"
  if [ ! -f $file ]; then
    echo "$file not found!"
    exit 1
fi
done

## Creating output directory recursively if not exists
mkdir -p $odir
echo


#################################################
# SAM2PAT: Convert SAM to PAT
#################################################
# SAM files have been pre-filtered (since different alignment tools have different SAM formats, manual filtering is recommended). Here we are calling a Perl program to convert them to PAT format.
echo "Number of SAM files: ${#samfiles[*]}"
echo ">>> $(date) - SAM2PAT (MAP_parseSAM2PAT.pl)"
patfiles=()
cd $indir
for file in "${samfiles[@]}";
do
  iname=$(basename $file .sam)
  #idir="$(dirname $file)/"
  oname=${iname}.PAT
  ## Make sure to modify the code location   seven
  perl MAP_parseSAM2PAT.pl -sam $file -poly T -m 20 -s 10 -more F
  #cp $file $oname ##DEBUG, remember to comment out this line during actual execution
  mv ${iname}.PAT $odir
  echo "$file >>> $oname"
  patfiles+=($oname)
done
cd $odir

###############################################################
# PAT2PA: Sort the same rows directly and calculate tag count
# PAT: chr, strand, coord				     
###############################################################

echo "Number of pat files: ${#patfiles[*]}"
echo ">>> $(date) - Converting PAT to PA (sorting and deduplication)"
pafiles=()
for file in "${patfiles[@]}"; do
  iname=$(basename $file .PAT)
  oname=${iname}.PA
  sort "$file" | uniq -c | awk '$1=$1' > "$oname"
  # PA: tagnum, chr, strand, coord
  echo "$file >>> $oname"
  pafiles+=("$oname")
done

cd "$odir"
# Convert PA to BED format (chr, start (end-1), end, ., count, strand; count is in the 5th column)
echo ">>> $(date) - Converting PA to BED (using awk)"
paBEDfiles=()
for file in "${pafiles[@]}"; do
  iname=$(basename "$file" .PA)
  oname=${iname}.PA.bed
  # Here, start = PAcoord-1
  awk -v OFS="\t" '{ print $2, $4-1, $4, ".", $1, $3 }' "$file" | sort-bed - > "$file.bed"
  paBEDfiles+=("$oname")
done

cd "$odir"

#################################################
# setIP
#################################################
echo "paBED files: ${#paBEDfiles[*]}"
echo ">>> $(date) - setIP (PAT_setIP_big.pl)"
realfiles=()
for file in "${paBEDfiles[@]}";
do
  oname=${file}.real
  #flds为chr,strand,coord的列下标（0为第1列）
  #0:5:2 对应bed中的相应列
  #注意添加对应的路径
  perl PAT_setIP_big.pl -in "$file" -skip 0 -suf "" -flds 0:5:2 -chr "$chrfa"
  ##cp $file $oname ##DEBUG，实际运行要注释掉这句
  echo "$file >>> $oname"
  realfiles+=($oname)
done

#################################################
# PA2PAC
#################################################
#合并所有样本的PA,用chr/strand/coord进行sort并取uniq
echo "real files: ${#realfiles[*]}"
echo ">>> $(date) - merge realPA files (cat>awk) >> all.PA.bed >> all.PA.uniq.bed"
cat ${realfiles[@]} > all.PA.bed
#cat *.PA.bed.real > all.PA.bed
sort-bed all.PA.bed > all.PA.bed.tmp
mv all.PA.bed.tmp  all.PA.bed

#计算所有样本PA的总tag数：处理成这种形式 chr1,100,101,+  1 , 再将第2列累加
awk -vOFS="\t" '{print $1","$2","$3","$6,$5}' all.PA.bed > all.PA.bed.tmp 
awk -F"\t" -vOFS="\t" '{a[$1]+=$2;}END{for(i in a)print i", "a[i];}' all.PA.bed.tmp > all.PA.bed.tmp2
#再还原成bed格式
awk -F'[,\t]' -vOFS="\t" '{print $1,$2,$3,".",$5,$4}' all.PA.bed.tmp2 |  sort-bed - > all.PA.uniq.bed
rm all.PA.bed all.PA.bed.tmp all.PA.bed.tmp2

#[szhou@login01 test1]$ cat all.PA.uniq.bed
#chr1    96      97      .        1      -  第5列是score

echo ">>> $(date) - PA2PAC (bedtools merge) >> all.PAC.bed"
## PA2PAC: 以dist距离合并
# 以下两个结果一样
#bedtools merge -i all.PA.bed -s -d $dist > all.PAC.bed
bedtools merge -i all.PA.uniq.bed -s -d $dist -c 6 -o distinct > all.PAC.bed
#bedtools merge -i all.PA.uniq.bed -s -d 24 -c 6 -o distinct > all.PAC.bed
#[szhou@login01 WXH]$ cat PAT01.txt.PA.bed (sorted)
#chr1    98      99      .       3       +
#chr1    98      99      .       3       -

#[szhou@login01 WXH]$ cat all.PAC.bed
#chr1    98      103     +
#chr1    98      103     -

#################################################
# countPA
# 比对每个样本的PA到PAC区域中
#################################################
#先将每个PA样本按方向划分
echo ">>> $(date) - split paBEDfiles by strand"
for file in "${realfiles[@]}";
do
  grep "+$" $file > "$file".for
  grep "\\-$" $file > "$file".rev
done



#PAC也按方向划分
echo ">>> $(date) - split all.PAC.bed by strand"
grep "+$" all.PAC.bed > all.PAC.bed.for
grep "\\-$" all.PAC.bed > all.PAC.bed.rev
rm all.PAC.bed

#总的uniq.PA也按方向划分
echo ">>> $(date) - split all.PA.uniq.bed by strand"
grep "+$" all.PA.uniq.bed > all.PA.uniq.bed.for
grep "\\-$" all.PA.uniq.bed > all.PA.uniq.bed.rev


## ----------------------------------------
#计算PAC的信息：总PA数、总PAT数、PAC坐标、refPA的PAT数
echo ">>> $(date) - PACinfo (bedmap --max-element) >> all.PAC.info"
bedmap --echo --count --delim  "\t" --sum --prec 0  --max-element all.PAC.bed.for all.PA.uniq.bed.for | awk -vOFS="\t" '{print $1,$2+1,$3,$4,$5,$6,$9,$11}' > all.PAC.info.for
bedmap --echo --count --delim  "\t" --sum --prec 0  --max-element all.PAC.bed.rev all.PA.uniq.bed.rev | awk -vOFS="\t" '{print $1,$2+1,$3,$4,$5,$6,$9,$11}' > all.PAC.info.rev
## 参数--count是uniqPA数；--sum是PAT总数
rm all.PA.uniq.bed.for all.PA.uniq.bed.rev
## chr1    98      103     +       5(PA数)       25(PAT数)      chr1    99      100(取得这1列就是PAC的坐标)     .       12(refPA的PAT数量)      +
#合并正反向
cat all.PAC.info.for all.PAC.info.rev > all.PAC.info
rm all.PAC.info.for all.PAC.info.rev

##[szhou@login01 test1]$ cat all.PAC.info.for
##chr1    98+1(首1)      103     +       5(uniqPA数)        25 (PAT数)     100(PAC坐标)     12(refPA的PAT数量)

##统计PAC在每个样本下的PA个数
#按方向比对PA到PAC （注意PA和PAC都已经sorted）
echo ">>> $(date) - countPAinPAC (bedmap --count)"
PAinPACfiles=()
for file in "${realfiles[@]}";
do
  bedmap --echo --count --delim "\t" all.PAC.bed.for "$file".for  > "$file".for.PAinPAC.bed
  bedmap --echo --count --delim "\t" all.PAC.bed.rev "$file".rev  > "$file".rev.PAinPAC.bed
  ##参数 --count  The number of overlapping elements in <map-file>.
  ##合并正反向
  cat "$file".for.PAinPAC.bed "$file".rev.PAinPAC.bed > "$file".PAinPAC.bed
  rm "$file".for.PAinPAC.bed "$file".rev.PAinPAC.bed  
  echo ">>> ""$file".PAinPAC.bed
  PAinPACfiles+=("$file".PAinPAC.bed)
done





#取得所有PAinPAC的第5列（计数列）
echo ">>> $(date) - mergeCount (awk) >>> all.PAC.PAcount"
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  ${PAinPACfiles[@]} | awk '{$1=$1}1' OFS="\t" > all.PAcount

#awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  *.PAinPAC.bed | awk '{$1=$1}1' OFS="\t" > all.PAcount



#连接PAC.info和all.PAcount
paste all.PAC.info all.PAcount > all.PAC.PAcount 
rm all.PAcount

#以下code不用
#取得前4列，并将start重新变成start+1（因为bed文件是0坐标，要改成1坐标）
#awk -vOFS="\t" '{ print $1,$2+1,$3,$4}' ${PAinPACfiles[0]}  > all.PAC.coord
#连接第1个文件的前4列和all.PAcount
##paste <(cut -f1,2,3,4 ${PAinPACfiles[0]}) <(cat all.PAcount) > all.PAC.PAcount 
#paste all.PAC.coord all.PAcount > all.PAC.PAcount 

#################################################
# countPAT
# 比对每个样本的PAT到PAC区域中
#################################################
##统计PAT的个数
#按方向比对PA到PAC （注意PA和PAC都要sorted）
echo ">>> $(date) - countPATinPAC (bedmap --sum)"
PAinPACfiles=()
for file in "${realfiles[@]}";
do
  bedmap --echo --sum --prec 0 --delim "\t" all.PAC.bed.for "$file".for  > "$file".for.PATinPAC.bed
  bedmap --echo --sum --prec 0 --delim "\t" all.PAC.bed.rev "$file".rev  > "$file".rev.PATinPAC.bed
  ##参数 --sum对score列求和，-prec保留0位小数
  ##合并正反向
  cat "$file".for.PATinPAC.bed "$file".rev.PATinPAC.bed > "$file".PATinPAC.bed
  rm "$file".for.PATinPAC.bed "$file".rev.PATinPAC.bed  "$file".for "$file".rev
  echo ">>> ""$file".PATinPAC.bed
  PAinPACfiles+=("$file".PATinPAC.bed)
done

rm all.PAC.bed.for all.PAC.bed.rev 


#取得所有PAinPAC的第5列计数列
echo ">>> $(date) - mergeCount (awk) >>> all.PAC.PATcount"
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  ${PAinPACfiles[@]} | awk '{$1=$1}1' OFS="\t" > all.PATcount

#awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  *.PATinPAC.bed | awk '{$1=$1}1' OFS="\t" > all.PATcount


#连接PAC.info和all.PAcount
paste all.PAC.info all.PATcount > all.PAC.PATcount 
rm all.PATcount

#以下code不用
#连接第1个文件的前4列和all.PAcount
##paste <(cut -f1,2,3,4 ${PAinPACfiles[0]}) <(cat all.PATcount) > all.PAC.PATcount 
#paste all.PAC.coord all.PATcount > all.PAC.PATcount 

#将统计时bedmap产生的NAN替换为0
sed "s/NAN/0/g" all.PAC.PATcount > all.PAC.PATcount1
mv all.PAC.PATcount1 all.PAC.PATcount

#################################################
# 输出展示
#################################################
#输出header文件：标题行，对应PAC矩阵的标题
echo ">>> $(date) - output header >>> all.PAC.header"
OLDIFS="$IFS"; IFS=$'\t'
header=(chr UPA_start UPA_end strand PAnum tot_tagnum coord refPAnum)
header+=("${lbls[@]}")
echo "${header[*]}" > all.PAC.header
IFS="$OLDIFS"

echo " ---------------------------------------------------- "
echo "[$odir]"
echo ">>> all.PAC.info     PAC坐标信息(1-base)"
echo ">>> all.PAC.PATcount PAC在所有样本下的PAT计数"
echo ">>> all.PAC.PAcount  PAC在所有样本下的PA计数(不常用)"
echo ">>> all.PAC.header   PAC矩阵的标题行"
echo ">>> all.PA.uniq.bed (0-base)   所有样本的PA总集bed格式坐标从0开始"
echo ">>> 单样本中间文件 --- *.PAT | .PA | PA.bed (0-base) | .real (0-base) | .IP (0-base) | .PAinPAC.bed (0-base)"

