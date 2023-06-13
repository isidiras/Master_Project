#I use this small script to split my data into batches in order
#to make my analysis faster

#!/bin/bash
i=1
k=0
number_of_batches=40
while [ $i -le $number_of_batches ]

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
