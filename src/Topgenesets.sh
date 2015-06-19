#! /bin/bash

topX=$2
topGenesets="../../topGenesets"$topX".txt"


output="../results/topGenesets"$topX".txt"
if [ -e $output ]
	then
	    rm $output
fi

raw_dir=`pwd`

while read -r line
do

	gseafile=$line
	# echo $gseafile
	cd $gseafile
	# echo "#####################################"


	# printf "top %d gene sets \n" $topX
	for dir in $(ls -d *Gsea*)
	do
		# echo "###" $dir
		key=${gseafile:11:3}${gseafile:20}${dir/_*/}
		# pwd
		# echo $key
		# echo $topGenesets
		# echo $(ls $dir/gsea*pos*.xls)
		awk 'BEGIN{FS=OFS="\t"; printf "'$key'"}NR!=1 && NR<=("'$topX'" + 1){printf "\t" $1, NR }END{printf "\n" }' $(ls $dir/gsea*pos*.xls)

	done
	# echo "#####################################"
	cd $raw_dir
	
done < $1


