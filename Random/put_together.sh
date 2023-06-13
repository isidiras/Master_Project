#I use this small script to combine the results
#from my batches.
 
 #!/bin/bash

 source /cvmfs/alice.cern.ch/etc/login.sh
 eval $(alienv printenv VO_ALICE@pythia::v8304-9,VO_ALICE@ROOT::v6-24-06-18)

 INPUT_FILE=/data/alice/isidiras/ccbar_MONASH_Hard_low

 cd ${INPUT_FILE}
 echo ${INPUT_FILE}
 mkdir -p complete_root
 echo "Here we go!"
 hadd complete_root/DplusDplus.root  MONASH_BATCH_STATUS*/DplusDplus.root
 hadd complete_root/DplusDminus.root  MONASH_BATCH_STATUS*/DplusDminus.root
 hadd complete_root/DplusLplus.root  MONASH_BATCH_STATUS*/DplusLplus.root
 hadd complete_root/DplusLminus.root  MONASH_BATCH_STATUS*/DplusLminus.root
 hadd complete_root/DminusDminus.root  MONASH_BATCH_STATUS*/DminusDminus.root
 hadd complete_root/DminusDplus.root  MONASH_BATCH_STATUS*/DminusDplus.root
 hadd complete_root/DminusLminus.root  MONASH_BATCH_STATUS*/DminusLminus.root
 hadd complete_root/LplusLplus.root  MONASH_BATCH_STATUS*/LplusLplus.root
 echo "Half way there!"
 hadd complete_root/LplusLminus.root  MONASH_BATCH_STATUS*/LplusLminus.root
 hadd complete_root/LplusDplus.root  MONASH_BATCH_STATUS*/LplusDplus.root
 hadd complete_root/LplusDminus.root  MONASH_BATCH_STATUS*/LplusDminus.root
 hadd complete_root/LminusLminus.root  MONASH_BATCH_STATUS*/LminusLminus.root
 hadd complete_root/LminusLplus.root  MONASH_BATCH_STATUS*/LminusLplus.root
 hadd complete_root/LminusDminus.root  MONASH_BATCH_STATUS*/LminusDminus.root
 hadd complete_root/LminusDplus.root  MONASH_BATCH_STATUS*/LminusDplus.root

 echo "Job Well Done?"
