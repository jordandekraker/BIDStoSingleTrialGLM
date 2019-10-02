#!/bin/bash
usage(){
echo "usage: LSS_submit <bids_dir> <exp_start volume>"
}
bids_dir=$1
expstart_vol=$2
if [ "$#" -eq 1 ]
then
	mkdir -p $bids_dir/derivatives/singletrial_GLM
	subs=$(cut -f 1 $bids_dir/participants.tsv|tail -n +2)
	for sub in $subs
	do
	sub_dir="$bids_dir/derivatives/singletrial_GLM/$sub"
	mkdir -p $sub_dir 
	neurogliaSubmit -t -j LongSkinny 'matlab -nosplash -nodisplay -r "addpath(genpath("'/project/6007967/hyang336/matlab/'")); singless_LSS("'$bids_dir,$sub,$expstart_vol'"); exit;"' 
	done
else
usage
exit 1
fi

