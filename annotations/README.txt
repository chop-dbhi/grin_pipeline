cat nexterarapidcapture_exome_targetedregions_v1.2.bed | sed s/chr//g | sed s/M/MT/g > nexterarapidcapture_exome_targetedregions_v1.2.no_chr.MT.bed
cat NexteraRapidCapture_Exome_Probes_v1.2.txt | grep CEX | sed -e 's/chrM/chrMT/g;s/chr//g;' | cut -f2,3,4 > NexteraRapidCapture_Exome_Probes_v1.2.bed
picard BedToIntervalList I=annotations/NexteraRapidCapture_Exome_Probes_v1.2.bed O=annotations/NexteraRapidCapture_Exome_Probes_v1.2.interval_list SD=/mnt/isilon/cbmi/variome/reference/human/g1k_v37/human_g1k_v37.dict 
picard BedToIntervalList I=annotations/nexterarapidcapture_exome_targetedregions_v1.2.no_chr.MT.bed O=annotations/nexterarapidcapture_exome_targetedregions_v1.2.no_chr.MT.interval_list SD=/mnt/isilon/cbmi/variome/reference/human/g1k_v37/human_g1k_v37.dict 
