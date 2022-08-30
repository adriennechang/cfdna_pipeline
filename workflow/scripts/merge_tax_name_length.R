library(dplyr)

args=commandArgs(trailingOnly=TRUE)

taxid_name = read.delim(args[1], sep=",", header=F)
taxid_length = read.table(args[2], sep="\t")
tax_tax = read.table(args[3], sep=" ")

taxid_name$V1 = gsub(taxid_name$V1, pattern="-", replacement=NA)
taxid_name = taxid_name[!is.na(taxid_name$V1),]
taxid_name$V1 = trimws(taxid_name$V1, which="both")

file3= full_join(taxid_name, taxid_length, by="V1")
file3 = file3[!is.na(file3$V2.y),]
file3$V2.x = gsub(file3$V2.x, pattern="^", replacement="\"")
file3$V2.x = gsub(file3$V2.x, pattern="$", replacement="\"")
colnames(file3)=c("Taxid","Name","Length")
file3 = file3[!is.na(file3$Name),]

file4 = merge(file3, tax_tax, by.x="Taxid", by.y="V1", all=T)
colnames(file4) = c("Taxid","Name","Length","Genus", "Species")

file4 = file4[!is.na(file4$Name),]
write.table(file4, args[4], sep="\t", quote=F, row.names = F, col.names = T)
