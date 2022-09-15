args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=5) {
  stop("input table, primary entities, secondary entities, output_table, combo_length must be given", call.=FALSE)
}

library(glue)
library(RareComb)

input_table = args[1] #"/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/meta_rarecomb_in.csv"
primary_entities = args[2] #"/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/primary.csv"
secondary_entities = args[3] #"/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/secondary.csv"
output_table = args[4] #"/data5/deepro/ukbiobank/analysis/lifestyle_factors/data/meta_rarecomb_out.csv"
combo_length = as.numeric(args[5])

# input df
input_df = read.table(input_table, sep=",", header=TRUE)
# primary and secondary entities
pent = readLines(primary_entities)
sent = readLines(secondary_entities)

# cases_sample_list <- subset(input_df, input_df[["Output_pheno"]] == 1, select = c("Sample_Name"))

# gencode
filename = '/data5/bx_reference/hg38/annotations/gene_annotations/GENCODE39/gencode.v39.parsed.genes.tsv'
gencode = read.table(filename, sep='\t', header=T)
gencode[,'gene_id_stripped'] = unlist(lapply(gencode[,'gene_id'], function(x) strsplit(x, '.', fixed=T)[[1]][1]))
gencode = gencode[!duplicated(gencode$gene_id_stripped),]
gencode = gencode[,c('gene_id_stripped', 'Chrom', 'Start', 'End')]
colnames(gencode) = c('Gene', 'Chrom', 'Start', 'End')


max_freq_threshold=0.25
pval_filter_threshold=0.05
adj_pval_type='bonferroni'
min_power_threshold=0.7
sample_names_ind='Y'
quiet=FALSE
min_indv_threshold=5

result = compare_enrichment_modifiers(
    input_df, 
    gene_coordinates_df=gencode, 
    primary_input_entities=pent, 
    secondary_input_entities=sent, 
    combo_length=combo_length, 
    min_indv_threshold=min_indv_threshold, 
    max_freq_threshold=max_freq_threshold, 
    input_format='Input_', 
    output_format='Output_', 
    pval_filter_threshold=pval_filter_threshold, 
    adj_pval_type=adj_pval_type, 
    min_power_threshold=min_power_threshold, 
    sample_names_ind=sample_names_ind, 
    ld_block_size=0, 
    quiet=quiet
    )

write.csv(result, output_table, row.names=F)

rm(input_df)
rm(pent)
rm(sent)
rm(result)
