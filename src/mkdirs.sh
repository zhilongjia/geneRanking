#! /bin/bash

for i in "AGR" "BGI" "NWU" "PSU"
do
	mkdir -p "../results/"$i"/gsea_c2cp"
	mkdir -p "../results/"$i"/gsea_c2kegg"
	mkdir -p "../results/"$i"/gsea_c5bp"
done
