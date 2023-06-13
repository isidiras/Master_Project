#!/bin/bash
i=1
k=0
while [ $i -le 40 ]

do
        mkdir MONASH_BATCH_STATUS$i
	for z in {1..100}
	do
		 gcount=$(($z+$k))
       		 cp -r  /data/alice/isidiras/ccbar_MONASH_Hard_low/output_MONASH_STATUS1/Group$gcount /data/alice/isidiras/ccbar_MONASH_Hard_low/MONASH_BATCH_STATUS$i/
       	done
	k=$k+100
	 echo "Batch $i has been greated!"
        cd /data/alice/isidiras/ccbar_MONASH_Hard_low/MONASH_BATCH_STATUS$i/
	if [ $i -gt 1 ] 
	then
	j=1
       	 for file in *
       	        do
                mv -v $file Group$j
                let j++
                done
	fi
        cd ..
        let i++
done
