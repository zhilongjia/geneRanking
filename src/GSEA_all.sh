#! /bin/bash

#This script aims to do gsea for different ranking methods for RNA-Seq and Microarray


#delete the old results if needed.
while read -r line
do
	oldfile=$line"/*"
	rm -rf $oldfile
done < $1



###########################################################
#use rank of other methods to do GSEA

gseapreranked () {
	for i in  "../results/AGR/" "../results/BGI/" "../results/NWU/" "../results/PSU/"
do
	if [ $1 = "../data/c2.cp.v4.0.symbols.gmt" ]
		then
		out=$i"gsea_c2cp"
    elif [ $1 = "../data/c2.cp.kegg.v4.0.symbols.gmt" ]
        then
        out=$i"gsea_c2kegg"
	elif [ $1 = "../data/c5.bp.v4.0.symbols.gmt" ]
		then
		out=$i"gsea_c5bp"
	fi
    
    for j in `ls $i*".rnk"`
    do
    	label=`basename $j`
        label=${label%.*}
    	#get the label
    	pwd
        echo $1 "," $j "," $label
    	#$1 is geneset 
    	java -cp ./gsea2-2.1.0.jar -Xmx1024m xtools.gsea.GseaPreranked -gmx $1 -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk $j -scoring_scheme classic -rpt_label $label -include_only_symbols true -make_sets false -plot_top_x 0 -rnd_seed 149 -set_max 500 -set_min 15 -zip_report false -out $out -gui false

    done
done
}

gseapreranked "../data/c2.cp.v4.0.symbols.gmt"
gseapreranked "../data/c2.cp.kegg.v4.0.symbols.gmt"
gseapreranked "../data/c5.bp.v4.0.symbols.gmt"

