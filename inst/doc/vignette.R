## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(vivaldi)
library(kableExtra)
library(tidyverse)

## -----------------------------------------------------------------------------
vardir = system.file("extdata", "vcfs", package="vivaldi") 

## -----------------------------------------------------------------------------
seg_sizes = system.file("extdata", "SegmentSize.csv", package="vivaldi")
sizes = read.csv(file=seg_sizes,header=T,sep=",",na.strings = c(''))

#select only the relevant segment sizes for H1N1
sizes = filter(sizes, STRAIN ==  "H1N1")

genome_size = 13133

## -----------------------------------------------------------------------------
rep_info = system.file("extdata", "reps.csv", package="vivaldi")
replicates = read.csv(file = rep_info, header = T, sep = ",", na.strings = c(""))
kable(head(replicates))
dim(replicates)

## ----message=FALSE, warning=FALSE---------------------------------------------
VCF_DF = arrange_data(vardir, ref = system.file("extdata", "H1N1.fa", package="vivaldi"), annotated = 'yes')
kable(head(VCF_DF))
VCF_DF %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  kable()
dim(VCF_DF)

## -----------------------------------------------------------------------------
SEGMENTS = c("H1N1_PB2","H1N1_PB1","H1N1_PA","H1N1_HA","H1N1_NP","H1N1_NA","H1N1_MP","H1N1_NS")
VCF_DF$CHROM = factor(VCF_DF$CHROM, levels = SEGMENTS)

## -----------------------------------------------------------------------------
cols = c("sample","CHROM","POS","REF","ALT","ANN","ALT_TYPE","major","minor")

DF_reps = merge_replicates(VCF_DF,replicates,"rep1","rep2",cols)

kable(head(DF_reps))

DF_reps %>%
  group_by(sample) %>%
  summarise(n = n()) %>%
  kable()

dim(DF_reps)

## -----------------------------------------------------------------------------
df = merge(replicates,VCF_DF, by.x = c("filename"), by.y = c("sample"))

df_rep1 = dplyr::filter(df, replicate == "rep1")
df_rep2 = dplyr::filter(df, replicate == "rep2")

df_merged_keep = merge(df_rep1, df_rep2, by = cols, all = TRUE)
df_merged_keep = df_merged_keep[!duplicated(df_merged_keep), ]

df_merged_keep$minorfreq.x[is.na(df_merged_keep$minorfreq.x)] = 0
df_merged_keep$minorfreq.y[is.na(df_merged_keep$minorfreq.y)] = 0

ggplot2::ggplot(df_merged_keep, ggplot2::aes(x = minorfreq.x, y = minorfreq.y)) + 
  ggplot2::geom_point()

## -----------------------------------------------------------------------------
ggplot2::ggplot(DF_reps, ggplot2::aes(x = minorfreq.x, y = minorfreq.y)) + 
  ggplot2::geom_point()

## -----------------------------------------------------------------------------
ggplot2::ggplot(DF_reps, ggplot2::aes(x = minorfreq, y = weighted_minorfreq)) + ggplot2::geom_point()

## -----------------------------------------------------------------------------
# Default coverage (200) and frequency (0.03) cutoffs 
#DF_filt = filter_variants(DF_reps)

# To run with custom values, specify these in the function
DF_filt = filter_variants(DF_reps, coverage_cutoff = 0, frequency_cutoff = 0.01 )

kable(head(DF_filt))

dim(DF_filt)

## -----------------------------------------------------------------------------
DF_filt = prepare_annotations(DF_filt)
kable(head(DF_filt))
dim(DF_filt)

## -----------------------------------------------------------------------------
DF_filt_ns = filter(DF_filt, feature_id != "H1N1_NS.1" & feature_id != "H1N1_NS.2" & 
                      feature_id != "H1N1_M1.1" & feature_id != "H1N1_M1.2")

DF_filt_s = filter(DF_filt, feature_id == "H1N1_NS.1" | feature_id == "H1N1_NS.2" | 
                          feature_id =="H1N1_M1.1" | feature_id =="H1N1_M1.2")

DF_filt_s_unique = DF_filt_s %>% dplyr::group_by(sample,CHROM,POS,REF,ALT) %>% 
  dplyr::mutate(count = 1, totalsamp = sum(count)) %>%
  filter(totalsamp == 1) %>%
  ungroup()

# if variants are duplicated, only take those from NS.1 or M.1
DF_filt_s_double = DF_filt_s %>% dplyr::group_by(sample,CHROM,POS,REF,ALT) %>% 
  dplyr::mutate(count = 1, totalsamp = sum(count)) %>%
  filter(totalsamp > 1) %>%
  filter(feature_id == "H1N1_NS.1" | feature_id =="H1N1_M1.1") %>%
  ungroup() %>%
  dplyr::select(!(count:totalsamp))
  
DF_filt_s_all = rbind(DF_filt_s_unique,DF_filt_s_double)
DF_filt_s_all = DF_filt_s_all[!duplicated(DF_filt_s_all), ] %>% droplevels()

DF_filt_SNVs = rbind(DF_filt_s_all,DF_filt_ns)

## -----------------------------------------------------------------------------
DF_filt_SNVs = add_metadata(DF_filt_SNVs, sizes, c('CHROM'), c('segment'))

kable(head(DF_filt_SNVs))
dim(DF_filt_SNVs)

## -----------------------------------------------------------------------------
af_distribution(DF_filt_SNVs)

## -----------------------------------------------------------------------------
group_list_seg = c('sample','CHROM', "SegmentSize")
seg_count = tally_it(DF_filt_SNVs, group_list_seg, "snv_count")

kable(seg_count)

## -----------------------------------------------------------------------------
group_list_gen = c('sample')
gen_count = tally_it(DF_filt_SNVs, group_list_gen, "snv_count")

kable(gen_count)

## -----------------------------------------------------------------------------
snv_location(DF_filt_SNVs)

## -----------------------------------------------------------------------------
snv_genome(DF_filt_SNVs)

## -----------------------------------------------------------------------------
snv_segment(DF_filt_SNVs)

## -----------------------------------------------------------------------------
DF_tstv = tstv_ratio(DF_filt_SNVs,genome_size)
kable(head(DF_tstv))

## -----------------------------------------------------------------------------
tstv_plot(DF_tstv)

## -----------------------------------------------------------------------------
DF_filt_SNVs = shannon_entropy(DF_filt_SNVs,genome_size)

kable(head(DF_filt_SNVs))
dim(DF_filt_SNVs)

## -----------------------------------------------------------------------------
plot_shannon(DF_filt_SNVs)

## -----------------------------------------------------------------------------
SPLICEFORMS = c("H1N1_PB2.1", "H1N1_PB1.1", "H1N1_PA.1", "H1N1_HA.1" ,"H1N1_NP.1", "H1N1_NA.1", "H1N1_M1.1", "H1N1_M1.2", "H1N1_NS.1", "H1N1_NS.2")
dNdS_segment(DF_filt)

## -----------------------------------------------------------------------------
shared_snv_plot(DF_filt_SNVs)

## -----------------------------------------------------------------------------
shared_snv_plot(DF_filt_SNVs, samples = c("a_1_fb","a_1_iv"))

## -----------------------------------------------------------------------------
shared_snv_table(DF_filt_SNVs) %>%
  head() %>%
  kable()

## -----------------------------------------------------------------------------
position_allele_freq(DF_filt_SNVs,"H1N1_NP", "1247")

