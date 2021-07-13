#!/usr//bin/bash
# Mirror Project Data from Functional Genomics Center Zurich to NAS drive

#set -e
dir=/nfs/hardt_raw/FGCZ
log=${dir}/log

year=`date +%Y`
#date=`date +%Y%m%d`
date=`date`
project=p2379

cd ${dir}
[ -d ${log} ] || mkdir -p ${log} 
wget -m -e --robots=off --reject="index.html*" --no-parent --user=nguyenb --password=Salmonella101 http://fgcz-gstore.uzh.ch/projects/${project}/ 
echo -e "Script was executed on $date" >> ${log}/${year}.log
#echo -e "Script was executed on $date" >> /science/ssunagaw/test.log
