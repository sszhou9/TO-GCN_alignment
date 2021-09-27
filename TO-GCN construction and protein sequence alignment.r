#TO-GCN construction and protein sequence alignment
#UV-B_TO-GCN 
Cutoff 4 4 UV-B_TF_physiology.txt UV-B_physiology.txt
GCN 4 4 UV-B_TF_physiology.txt UV-B_physiology.txt 0.84 0.84 -1 -1
sed 's/ //g' C1+C2+.csv | awk -F "," '!(SEEN[$1,$3]++) && !(($3,$1) in SEEN)' | awk -F "," '{print $1"\t"$3"\t"$4}'> UV-B_0.84_C1+C2+.txt
head UV-B_0.84_C1+C2+.txt
Cutoff 4 4 UV-B_TF.txt UV-B.txt
GCN 4 4 UV-B_TF.txt UV-B.txt 0.95 0.95 -1 -1
sed 's/ //g' C1+C2+.csv | awk -F "," '!(SEEN[$1,$3]++) && !(($3,$1) in SEEN)' | awk -F "," '{print $1"\t"$3"\t"$4}'> UV-B_0.95_C1+C2+.txt
head UV-B_0.95_C1+C2+.txt
R
rm(list=ls())
library("dplyr")
UV-B_names <- c("Fv/Fm","Qy","Chl_a","Chl_b","Carotenoids","Cond","Photo","Phenolics","Flavonoids","Procyanidines","H2O2","GSH","GST","T-AOC")
UV-B_all <- read.table("UV-B_0.84_C1+C2+.txt",header = T)
UV-B_all_1 <- UV-B_all[(UV-B_all[,1] %in% UV-B_names)&(!UV-B_all[,2] %in% UV-B_names),]
c <- read.table("UV-B_TF.txt",row.names = 1)
UV-B_all_1 <- UV-B_all_1[!UV-B_all_1[,2] %in% rownames(c),]
write.table(UV-B_all_1,file="UV-B_physiology.txt",quote = F,row.names = F,col.names = F,sep="\t")
UV-B_all_2 <- read.table("UV-B_0.95_C1+C2+.txt",header = T)
UV-B_all <- rbind(UV-B_all_1,UV-B_all_2)
write.table(UV-B_all,file="UV-B_UV-B_all_gene_physiology.txt",quote = F,row.names = F,col.names = F,sep="\t")
head(UV-B_all)
q()
awk 'NR == FNR{UV-B_all[$1]=$2;next}{print $0"\t"UV-B_all[$2]}' UV-B_replacement_rule.txt UV-B_UV-B_all_gene_physiology.txt|awk 'NF == 3{print $0}NF == 4{print $1"\t"$4"\t"$3}' > replacement.txt
awk 'NR == FNR{UV-B_all[$1]=$2;next}{print $0"\t"UV-B_all[$1]}' UV-B_replacement_rule.txt replacement.txt | awk 'NF == 3{print $0}NF == 4{print $4"\t"$2"\t"$3}' > blast_UV-B_UV-B_all_gene_physiology.txt
awk '{print "UV-B_"$1"\tUV-B_"$2"\t"$3}' blast_UV-B_UV-B_all_gene_physiology.txt > UV-B_UV-B_all_gene_physiology.txt
rm replacement.txt

#UV-C_TO-GCN
Cutoff 4 4 UV-C_TF_physiology.txt UV-C_physiology.txt
GCN 4 4 UV-C_TF_physiology.txt UV-C_physiology.txt 0.84 0.84 -1.00 -1.00
sed 's/ //g' C1+C2+.csv | awk -F "," '!(SEEN[$1,$3]++) && !(($3,$1) in SEEN)' | awk -F "," '{print $1"\t"$3"\t"$4}'> UV-C_0.84_C1+C2+.txt
head UV-C_0.84_C1+C2+.txt
Cutoff 4 4 97UV-C_TF.txt 97UV-C.txt
GCN 4 4 97UV-C_TF.txt 97UV-C.txt 0.95 0.95 -1.00 -1.00
sed 's/ //g' C1+C2+.csv | awk -F "," '!(SEEN[$1,$3]++) && !(($3,$1) in SEEN)' | awk -F "," '{print $1"\t"$3"\t"$4}'> UV-C_0.95_C1+C2+.txt
head UV-C_0.95_C1+C2+.txt
R
rm(list=ls())
library("dplyr")
UV-C_names <- c("Fv/Fm","Qy","Chl_a","Chl_b","Carotenoids","Cond","Photo","Phenolics","Flavonoids","Procyanidines","H2O2","GSH","GST","T-AOC")
UV-C_all <- read.table("UV-C_0.84_C1+C2+.txt",header = T)
UV-C_all_1 <- UV-C_all[(UV-C_all[,1] %in% UV-C_names)&(!UV-C_all[,2] %in% UV-C_names),]
c <- read.table("UV-C_TF.txt",row.names = 1)
UV-C_all_1 <- UV-C_all_1[!UV-C_all_1[,2] %in% rownames(c),]
write.table(UV-C_all_1,file="C_physiology.txt",quote = F,row.names = F,col.names = F,sep="\t")
UV-C_all_2 <- read.table("UV-C_0.95_C1+C2+.txt",header = T)
UV-C_all <- rbind(UV-C_all_1,UV-C_all_2)
write.table(UV-C_all,file="C_UV-C_all_gene_physiology.txt",quote = F,row.names = F,col.names = F,sep="\t")
head(UV-C_all)
q()
awk 'NR == FNR{UV-C_all[$1]=$2;next}{print $0"\t"UV-C_all[$2]}' UV-C_replacement_rule.txt C_UV-C_all_gene_physiology.txt|awk 'NF == 3{print $0}NF == 4{print $1"\t"$4"\t"$3}' >replace.txt
awk 'NR == FNR{UV-C_all[$1]=$2;next}{print $0"\t"UV-C_all[$1]}' UV-C_replacement_rule.txt replace.txt|awk 'NF == 3{print $0}NF == 4{print $4"\t"$2"\t"$3}' > blast_C_UV-C_all_gene_physiology.txt
awk '{print "UV-C_"$1"\tUV-C_"$2"\t"$3}' blast_C_UV-C_all_gene_physiology.txt >C_UV-C_all_gene_physiology.txt
less blast_C_all_gene_physiology.txt |cut -f 1-2|sed  's/\t/\n/'|sort|uniq|grep -v -E "Fv/Fm|Qy|Chl_a|Chl_b|Carotenoids|Cond|Photo|Phenolics|Flavonoids|Procyanidines|H2O2|GSH|GST|T-AOC|Trmmol" >blast_C.txt
less blast_UV-B_all_gene_physiology.txt |cut -f 1-2|sed  's/\t/\n/'|sort|uniq|grep -v -E "Fv/Fm|Qy|Chl_a|Chl_b|Carotenoids|Cond|Photo|Phenolics|Flavonoids|Procyanidines|H2O2|GSH|GST|T-AOC|Trmmol" >blast_B.txt
faSomeRecords pita.HQLQgenes.proteins.fasta blast_C.txt C_FASTA.txt
faSomeRecords pita.HQLQgenes.proteins.fasta blast_B.txt B_FASTA.txt
makeblastdb -in C_FASTA.txt  -dbtype prot
blastp -query B_FASTA.txt -db C_FASTA.txt -evalue 1e-10 -num_threads 20 -max_target_seqs 20 -out B_C.evals.txt -outfmt "6 qseqid sseqid pident evalue"
awk '{print $1"\t"$2}' B_C.evals.txt >B_C.evals2.txt
awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)' B_C.evals2.txt > blast.txt
rm B_C.evals2.txt
awk '{print "UV-B_"$1"\tUV-C_"$2"\t"$3}' B_C.evals.txt >B_C.evalsreplace.txt
cat UV-B_all_gene_physiology.txt C_all_gene_physiology.txt>B_C_all_gene_physiology.txt
cat UV-B_all_gene_physiology.txt C_all_gene_physiology.txt B_C.evalsreplace.txt >B_C_evals_all_gene_physiology.txt