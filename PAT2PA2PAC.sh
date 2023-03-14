#!/bin/bash -l

#------------------------------------------------- 
## Ex. 
## dos2unix PAT2PA2PAC.sh
## bash PAT2PA2PAC.sh "P_heterocycla.fa" "/xx/" "/xx/yy" 24 "1.sam,2.sam" "s1,s2"

## require:
## MAP_parseSAM2PAT.pl PAT_setIP_big.pl bedtools bedmap awk sed

#  ����
#1 chrfa="P_heterocycla.fa" #Ⱦɫ�����У�����setIP
#2 indir=/xx/ #����·��
#3 odir=/xx/yy/ #������Ŀ¼
#4 dist=24 #PA��PAC���ھ���
#5 samfiles="1.sam,2.sam" #Ҫ�����SAM�ļ���ע�⣺�Ѿ����ݲ�ͬmapping���ߵĽ�����˳���Ҫ��SAM������ֻ�ǵ���perl��ת��PAT��ʽ
#6 lbls="test1,test2" #��Ӧ�ı�ǩ

## ��� [odir] ע��start������0����1��ʼ
#all.PAC.info           -- PAC��info��Ϣ
#all.PAC.header         -- PAC��������� [chr UPA_start(1-base) UPA_end strand PAnum tot_tagnum coord(1-base) refPAnum]/s1...sN
#all.PAC.PAcount        -- PAC��ÿ�������µ�PA���� [info]/s1...sN �������ĺ�ҪС�ڵ���PAnum��
#all.PAC.PATcount       -- PAC��ÿ�������µ�PAT���� [info]/s1...sN  �������ĺ�Ҫ����tot_tagnum��
#all.PA.uniq.bed             -- ������������PA����(ȥ�غ��Ҽ����ܵ�PAT��) chr/start(0-base)/end/./PATnum/strand
#----------- ����Ϊ���������м��ļ� ----------- 
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

## DEBUG ���� ע������е� DEBUG �� ##
#chrfa="P_heterocycla.fa"
#indir="/data/gpfs01/szhou/WXH"
#odir="/data/gpfs01/szhou/WXH/test1/"
#dist=24
#samfiles="PAT01.txt,PAT02.txt"
#lbls="s1,s2"

printf  "chrfa=$chrfa\nindir=$indir\nodir=$odir\ndist=$dist\nsamfiles=$temp\nlbls=$temp\n"

## ��������
cd $indir

for file in $temp
do
 echo "${file}"
 samfiles+=($file)
 lbls+=($file)
done
#echo "print $samfiles"

#samfiles=(${samfiles//[, ;]/ });
#lbls=(${lbls//[, ;]/ });

##samfiles=$(find $indir -maxdepth 1 -name "*.bam*" -type f)

## �ж��ļ��Ƿ����
if [ ! -f $chrfa ]; then
    echo "$chrfa not found!"
	exit 1
fi

for file in "${samfiles[@]}";
do
  echo "********check $file exist************"
  if [ ! -f $file ]; then
    echo "$file not found!"
	exit 1
fi
done


## ���Ŀ¼,�����ݹ�Ŀ¼
mkdir -p $odir
echo "output results to $odir"

#################################################
# SAM2PAT: ��SAMתPAT
#################################################
# SAM���Ѿ����˹�����Ϊ��ͬ�ȶԹ���SAM��ʽ��ͬ�������ֶ�����Ϊ�ã�������ֻ�ǵ���perl����תΪPAT
echo "sam files: ${#samfiles[*]}"
echo ">>> $(date) - SAM2PAT (MAP_parseSAM2PAT.pl)"
patfiles=()
cd $indir
for file in "${samfiles[@]}";
do
  iname=$(basename $file .sam)
  #idir="$(dirname $file)/"
  oname=${iname}.PAT  
  ##ע���޸Ĵ�������·��   seven
  perl MAP_parseSAM2PAT.pl -sam $file -poly T -m 20 -s 10 -more F
  #cp $file $oname ##DEBUG��ʵ������Ҫע�͵����
  mv  ${iname}.PAT  $odir
  echo "$file >>> $oname"
  patfiles+=($oname)
done



cd $odir
#################################################
# PAT2PA��ֱ��sort��ͬ�в���tag count
# PAT: chr, strand, coord
#################################################
echo "pat files: ${#patfiles[*]}"
echo ">>> $(date) - PAT2PA (sort uniq)"
pafiles=()
for file in "${patfiles[@]}";
do
  iname=$(basename $file .PAT)
  oname=${iname}.PA
  sort $file | uniq -c | awk '$1=$1' > $oname
  # PA: tagnum, chr, strand, coord
  echo "$file >>> $oname"
  pafiles+=($oname)
done

cd $odir
#PAתbed��ʽ (chr, start(Ϊend-1), end, ., count, strand����5��Ϊcount��)
echo ">>> $(date) - PA2BED (awk)"
paBEDfiles=()
for file in "${pafiles[@]}";
do
  iname=$(basename $file .PA)
  oname=${iname}.PA.bed
  #����start=PAcood-1
  awk -vOFS="\t" '{ print $2, $4-1, $4, ".", $1, $3 }' $file | sort-bed - > $file.bed
  paBEDfiles+=($oname)
done

cd $odir
#################################################
# setIP
#################################################
echo "paBED files: ${#paBEDfiles[*]}"
echo ">>> $(date) - setIP (PAT_setIP_big.pl)"
realfiles=()
for file in "${paBEDfiles[@]}";
do
  oname=${file}.real
  #fldsΪchr,strand,coord�����±꣨0Ϊ��1�У�
  #0:5:2 ��Ӧbed�е���Ӧ��
  #ע����Ӷ�Ӧ��·��
  perl PAT_setIP_big.pl -in "$file" -skip 0 -suf "" -flds 0:5:2 -chr "$chrfa"
  ##cp $file $oname ##DEBUG��ʵ������Ҫע�͵����
  echo "$file >>> $oname"
  realfiles+=($oname)
done

#################################################
# PA2PAC
#################################################
#�ϲ�����������PA,��chr/strand/coord����sort��ȡuniq
echo "real files: ${#realfiles[*]}"
echo ">>> $(date) - merge realPA files (cat>awk) >> all.PA.bed >> all.PA.uniq.bed"
cat ${realfiles[@]} > all.PA.bed
#cat *.PA.bed.real > all.PA.bed
sort-bed all.PA.bed > all.PA.bed.tmp
mv all.PA.bed.tmp  all.PA.bed

#������������PA����tag���������������ʽ chr1,100,101,+  1 , �ٽ���2���ۼ�
awk -vOFS="\t" '{print $1","$2","$3","$6,$5}' all.PA.bed > all.PA.bed.tmp 
awk -F"\t" -vOFS="\t" '{a[$1]+=$2;}END{for(i in a)print i", "a[i];}' all.PA.bed.tmp > all.PA.bed.tmp2
#�ٻ�ԭ��bed��ʽ
awk -F'[,\t]' -vOFS="\t" '{print $1,$2,$3,".",$5,$4}' all.PA.bed.tmp2 |  sort-bed - > all.PA.uniq.bed
rm all.PA.bed all.PA.bed.tmp all.PA.bed.tmp2

#[szhou@login01 test1]$ cat all.PA.uniq.bed
#chr1    96      97      .        1      -  ��5����score

echo ">>> $(date) - PA2PAC (bedtools merge) >> all.PAC.bed"
## PA2PAC: ��dist����ϲ�
# �����������һ��
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
# �ȶ�ÿ��������PA��PAC������
#################################################
#�Ƚ�ÿ��PA���������򻮷�
echo ">>> $(date) - split paBEDfiles by strand"
for file in "${realfiles[@]}";
do
  grep "+$" $file > "$file".for
  grep "\\-$" $file > "$file".rev
done



#PACҲ�����򻮷�
echo ">>> $(date) - split all.PAC.bed by strand"
grep "+$" all.PAC.bed > all.PAC.bed.for
grep "\\-$" all.PAC.bed > all.PAC.bed.rev
rm all.PAC.bed

#�ܵ�uniq.PAҲ�����򻮷�
echo ">>> $(date) - split all.PA.uniq.bed by strand"
grep "+$" all.PA.uniq.bed > all.PA.uniq.bed.for
grep "\\-$" all.PA.uniq.bed > all.PA.uniq.bed.rev


## ----------------------------------------
#����PAC����Ϣ����PA������PAT����PAC���ꡢrefPA��PAT��
echo ">>> $(date) - PACinfo (bedmap --max-element) >> all.PAC.info"
bedmap --echo --count --delim  "\t" --sum --prec 0  --max-element all.PAC.bed.for all.PA.uniq.bed.for | awk -vOFS="\t" '{print $1,$2+1,$3,$4,$5,$6,$9,$11}' > all.PAC.info.for
bedmap --echo --count --delim  "\t" --sum --prec 0  --max-element all.PAC.bed.rev all.PA.uniq.bed.rev | awk -vOFS="\t" '{print $1,$2+1,$3,$4,$5,$6,$9,$11}' > all.PAC.info.rev
## ����--count��uniqPA����--sum��PAT����
rm all.PA.uniq.bed.for all.PA.uniq.bed.rev
## chr1    98      103     +       5(PA��)       25(PAT��)      chr1    99      100(ȡ����1�о���PAC������)     .       12(refPA��PAT����)      +
#�ϲ�������
cat all.PAC.info.for all.PAC.info.rev > all.PAC.info
rm all.PAC.info.for all.PAC.info.rev

##[szhou@login01 test1]$ cat all.PAC.info.for
##chr1    98+1(��1)      103     +       5(uniqPA��)        25 (PAT��)     100(PAC����)     12(refPA��PAT����)

##ͳ��PAC��ÿ�������µ�PA����
#������ȶ�PA��PAC ��ע��PA��PAC���Ѿ�sorted��
echo ">>> $(date) - countPAinPAC (bedmap --count)"
PAinPACfiles=()
for file in "${realfiles[@]}";
do
  bedmap --echo --count --delim "\t" all.PAC.bed.for "$file".for  > "$file".for.PAinPAC.bed
  bedmap --echo --count --delim "\t" all.PAC.bed.rev "$file".rev  > "$file".rev.PAinPAC.bed
  ##���� --count  The number of overlapping elements in <map-file>.
  ##�ϲ�������
  cat "$file".for.PAinPAC.bed "$file".rev.PAinPAC.bed > "$file".PAinPAC.bed
  rm "$file".for.PAinPAC.bed "$file".rev.PAinPAC.bed  
  echo ">>> ""$file".PAinPAC.bed
  PAinPACfiles+=("$file".PAinPAC.bed)
done





#ȡ������PAinPAC�ĵ�5�У������У�
echo ">>> $(date) - mergeCount (awk) >>> all.PAC.PAcount"
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  ${PAinPACfiles[@]} | awk '{$1=$1}1' OFS="\t" > all.PAcount

#awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  *.PAinPAC.bed | awk '{$1=$1}1' OFS="\t" > all.PAcount



#����PAC.info��all.PAcount
paste all.PAC.info all.PAcount > all.PAC.PAcount 
rm all.PAcount

#����code����
#ȡ��ǰ4�У�����start���±��start+1����Ϊbed�ļ���0���꣬Ҫ�ĳ�1���꣩
#awk -vOFS="\t" '{ print $1,$2+1,$3,$4}' ${PAinPACfiles[0]}  > all.PAC.coord
#���ӵ�1���ļ���ǰ4�к�all.PAcount
##paste <(cut -f1,2,3,4 ${PAinPACfiles[0]}) <(cat all.PAcount) > all.PAC.PAcount 
#paste all.PAC.coord all.PAcount > all.PAC.PAcount 

#################################################
# countPAT
# �ȶ�ÿ��������PAT��PAC������
#################################################
##ͳ��PAT�ĸ���
#������ȶ�PA��PAC ��ע��PA��PAC��Ҫsorted��
echo ">>> $(date) - countPATinPAC (bedmap --sum)"
PAinPACfiles=()
for file in "${realfiles[@]}";
do
  bedmap --echo --sum --prec 0 --delim "\t" all.PAC.bed.for "$file".for  > "$file".for.PATinPAC.bed
  bedmap --echo --sum --prec 0 --delim "\t" all.PAC.bed.rev "$file".rev  > "$file".rev.PATinPAC.bed
  ##���� --sum��score����ͣ�-prec����0λС��
  ##�ϲ�������
  cat "$file".for.PATinPAC.bed "$file".rev.PATinPAC.bed > "$file".PATinPAC.bed
  rm "$file".for.PATinPAC.bed "$file".rev.PATinPAC.bed  "$file".for "$file".rev
  echo ">>> ""$file".PATinPAC.bed
  PAinPACfiles+=("$file".PATinPAC.bed)
done

rm all.PAC.bed.for all.PAC.bed.rev 


#ȡ������PAinPAC�ĵ�5�м�����
echo ">>> $(date) - mergeCount (awk) >>> all.PAC.PATcount"
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  ${PAinPACfiles[@]} | awk '{$1=$1}1' OFS="\t" > all.PATcount

#awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $5 } END { for(i=1;i<=FNR;i++) print a[i]}'  *.PATinPAC.bed | awk '{$1=$1}1' OFS="\t" > all.PATcount


#����PAC.info��all.PAcount
paste all.PAC.info all.PATcount > all.PAC.PATcount 
rm all.PATcount

#����code����
#���ӵ�1���ļ���ǰ4�к�all.PAcount
##paste <(cut -f1,2,3,4 ${PAinPACfiles[0]}) <(cat all.PATcount) > all.PAC.PATcount 
#paste all.PAC.coord all.PATcount > all.PAC.PATcount 

#��ͳ��ʱbedmap������NAN�滻Ϊ0
sed "s/NAN/0/g" all.PAC.PATcount > all.PAC.PATcount1
mv all.PAC.PATcount1 all.PAC.PATcount

#################################################
# ���չʾ
#################################################
#���header�ļ��������У���ӦPAC����ı���
echo ">>> $(date) - output header >>> all.PAC.header"
OLDIFS="$IFS"; IFS=$'\t'
header=(chr UPA_start UPA_end strand PAnum tot_tagnum coord refPAnum)
header+=("${lbls[@]}")
echo "${header[*]}" > all.PAC.header
IFS="$OLDIFS"

echo " ---------------------------------------------------- "
echo "[$odir]"
echo ">>> all.PAC.info     PAC������Ϣ(1-base)"
echo ">>> all.PAC.PATcount PAC�����������µ�PAT����"
echo ">>> all.PAC.PAcount  PAC�����������µ�PA����(������)"
echo ">>> all.PAC.header   PAC����ı�����"
echo ">>> all.PA.uniq.bed (0-base)   ����������PA�ܼ�bed��ʽ�����0��ʼ"
echo ">>> �������м��ļ� --- *.PAT | .PA | PA.bed (0-base) | .real (0-base) | .IP (0-base) | .PAinPAC.bed (0-base)"

