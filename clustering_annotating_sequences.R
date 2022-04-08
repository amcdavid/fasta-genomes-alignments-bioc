## ---- echo = TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, echo = TRUE)
suppressPackageStartupMessages({
  library(Biostrings)
library(dplyr)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(forcats)
library(ggdendro)
library(ggrepel)
})
# also need rtracklayer, msa
# BiocManager::install(c('Biostrings', 'dplyr', 'ggplot2', 
# 'stringr', 'GenomicRanges', 'rtracklayer', 'seqLogo', 'Mus.musculus', 'forcats', 'ggdendro', 'ggrepel'))

fasta = readDNAStringSet('data/ny_nt_2022.fasta')
genome = readDNAStringSet('data/sars-cov2-genome.fasta')


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fasta


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(width(fasta))
genome


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
anno = rtracklayer::import('data/sars-cov2.gff.gz')
anno_df = filter(as.data.frame(anno), type %in%  c('gene', 'five_prime_UTR', 'three_prime_UTR'))
anno_plot = ggplot(anno_df, aes(ymin = start, ymax = end, x= Name, color = type)) +
  geom_linerange(lwd = 3) + coord_flip() + theme_minimal()

anno_plot


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
anno_df = anno_df %>% mutate(Namef = fct_reorder(factor(Name), start, .fun = min))
(anno_plot %+% anno_df) + aes(x = Namef)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Mus.musculus)
my_entrezid = AnnotationDbi::select(org.Mm.eg.db, 'Itgax', columns = "ENTREZID", keytype = 'SYMBOL')
Itgax_region = exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, 'gene')[[my_entrezid$ENTREZID]]
Itgax_region


## usearch -usearch_global data/ny_nt_2022.fasta -db data/sars-cov2-genome.fasta \

## -nastout outputs/ny_nt.nast -strand plus -id .7

## 

## usearch -usearch_global data/ny_nt_2022.fasta -db data/sars-cov2-genome.fasta  \

## -nastout outputs/ny_nt.nast -strand plus -id .7


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
query = c('Txe rat can',
            'The rat can',
          'The ran')
reference = 'The cat ran'
pairwiseAlignment(query[1], reference, type = 'global')


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pairwiseAlignment(query[3], reference, type = 'global')



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pairwiseAlignment(query[1], reference, type = 'local')
pairwiseAlignment(query[2], reference, type = 'local')


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ralign = pairwiseAlignment(fasta[1:10], genome, type = 'overlap')
aligned(ralign,degap = FALSE, gapCode = '-', endgapCode = '+')


## kalign -i data/ny_nt_2022.fasta -o outputs/kalign_large.out


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aligned = readDNAStringSet('outputs/ny_nt.nast')
aligned


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
str(aligned)
slotNames(aligned)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
showMethods(classes = class(aligned), where = 'package:Biostrings')


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# First 20 bp of S gene
s_only = narrow(aligned, anno_df[which(anno_df$Name=='S'), 'start'], 
                end =  anno_df[which(anno_df$Name=='S'), 'start']+20)

s_only


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm = consensusMatrix(s_only,  baseOnly = TRUE)[1:4,]
normalized = t(t(cm)/colSums(cm))
seqLogo::seqLogo(normalized)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cm = consensusMatrix(aligned,  baseOnly = TRUE)[1:4,]
maf = tibble(maf = apply(cm, 2, function(x) 1-max(x)/sum(x)), pos = seq_along(maf))
maf = maf %>% mutate(gene = anno_df[findInterval(pos, anno_df$start), 'Namef']) %>% filter(!is.na(gene))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
maf_plot = ggplot(maf, aes(x = pos, xend = pos, yend= maf))+ theme_minimal() + facet_grid( ~ gene, scales = 'free', space = 'free') 
maf_plot +geom_segment(y =0) + scale_y_sqrt()
maf_plot + aes(y = stats::filter(maf > 0.005, rep(1/100, 100))) + 
  geom_path() + ylab('> 0.5% polymorphism  per base')

maf_plot + aes(y = stats::filter(maf > 0.01, rep(1/100, 100))) + 
  geom_path() + ylab('> 1% polymorphisms per base')



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
genic_only = narrow(aligned, start = filter(anno_df, type == 'gene') %>% pull(start) %>% min(),
                    end =  filter(anno_df, type == 'gene') %>% pull(end) %>% max())
genic_only_wild = chartr(old = '-', new="N", genic_only)


## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------
## genic_only_wild = lapply(genic_only, function(x){
##     sx = as.character(x)
##     sl = str_locate(sx, '^-*')
##     sr = str_locate(sx, '-*$')
##     if(sr[2] >= sr[1]) str_sub(sx, sr[1], sr[2]) = paste(rep('N', sr[2]-sr[1]+1), collapse = '')
##     if(sl[2] >= sl[1])
##     str_sub(sx, sl[1]+1, sl[2]) = paste(rep('N', sl[2]-sl[1]), collapse = '')
##     as(sx, 'DNAString')
## })


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fastanames = tibble(Accession = names(genic_only))
meta = readr::read_csv('data/sequence_meta_2022.csv')
fastanames = left_join(fastanames, meta)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
stringDist_fixed = function(x){
  res = lapply(x, function(y) neditStartingAt(y, x, fixed = FALSE))
  as.dist(matrix(unlist(res), nrow = length(x), dimnames = list(names(x), names(x))))
  }


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dists = stringDist_fixed(genic_only_wild)
summary(dists)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(fastanames, aes(x = Collection_Date)) + 
  geom_histogram() + theme_minimal()


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hc = hclust(dists)
ddata = dendro_data(hc)
ddata$labels = left_join(ddata$labels, fastanames, by = c(label = 'Accession'))

p = ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
   geom_jitter(data = ddata$labels, aes(x = x, y = y, color = Collection_Date), size = .5) +
  theme_minimal() +
  scale_color_viridis_c()
p



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pr = prcomp(as.matrix(dists), retx = TRUE)

pc_meta = cbind(pr$x[, 1:2], fastanames)
plt = ggplot(pc_meta, aes(x = `PC1`, y= `PC2`, color = Collection_Date)) + 
  geom_jitter(alpha = .5, width = 10, height = 5) + theme_minimal() +
  scale_color_viridis_c()

plt

plt + aes(color = Pangolin) + 
  scale_color_discrete() + theme(legend.pos = 'none')




## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------
## library(Biostrings)
## library(readr)
## library(dplyr)
## library(lubridate)
## 
## if(!dir.exists('data')) dir.create('data')
## if(!file.exists('data/ny_nt_2022.fasta')){
##     download.file('https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=Completeness_s%7DCompleteness_s:(%22complete%22)&fq=&fq=%7B!tag=USAState_s%7DUSAState_s:(%22NY%22)&fq=%7B!tag=SLen_i%7DSLen_i:(%5B29000%20TO%2030000%5D)&cmd=download&sort=CollectionDate_s%20asc,SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=fasta&fl=AccVer_s,Definition_s,Nucleotide_seq', destfile = 'data/ny_nt_2022.fasta')
## }
## 
## if(!file.exists(('data/sequence_meta_2022.csv'))){
##   download.file('https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=Completeness_s%7DCompleteness_s:(%22complete%22)&fq=&fq=%7B!tag=USAState_s%7DUSAState_s:(%22NY%22)&fq=%7B!tag=SLen_i%7DSLen_i:(%5B29000%20TO%2030000%5D)&cmd=download&sort=CollectionDate_s%20asc,SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=csv&fl=Accession:id,Submitters:Authors_csv,Release_Date:CreateDate_dt,Pangolin:Lineage_s,PangoVersions,Random_Sampling:BaselineSurveillance_s,Isolate:IsolateParsed_s,Species:VirusSpecies_s,Molecule_type:GenomicMoltype_s,Length:SLen_i,Nuc_Completeness:Completeness_s,Geo_Location:CountryFull_s,Country:Country_s,USA:USAState_s,Host:Host_s,Isolation_Source:Isolation_csv,Collection_Date:CollectionDate_s', destfile = 'data/sequence_meta_2022.csv')
## }
## 
## 
## if(!file.exists('data/sars-cov2-genome.fasta')){
##   download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz', destfile = 'data/sars-cov2-genome.fasta.gz')
##   system2('gunzip', args = 'data/sars-cov2-genome.fasta.gz')
## }
## 
## if(!file.exists('data/sars-cov2.gff.gz')){
##   download.file('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz', destfile = 'data/sars-cov2.gff.gz')
## }
## 
## meta = read_csv("data/sequence_meta_2022.csv")
## meta_sub = meta %>% mutate(year_week = {
##   w = lubridate::week(Collection_Date)
##   y = lubridate::year(Collection_Date)
##   sprintf("%4.0d, week %02d", y, w)
## }) %>% group_by(year_week) %>% slice_sample(n = 5)
## 
## fasta_all = Biostrings::readDNAStringSet('~/Downloads/sequences.fasta')
## names(fasta_all) = stringr::str_extract(names(fasta_all), "^[0-9A-Z]+(;ins.+)?")
## fasta = fasta_all[meta_sub$Accession]
## Biostrings::writeXStringSet(fasta, 'data/ny_nt_2022.fasta')
## Biostrings::writeXStringSet(fasta[1:20], 'data/small.fasta')

