setwd("~/jolab001/backup_workspace/20200109-flavo_pyparanoid/")
set.seed(8675309)

library(dplyr)
library(reshape2)
library(ggplot2)
library(ggridges)
library(cowplot)
library(colorspace)
fj.ani <- read.table("ANI_fj.tsv", header=T, sep="\t")

ggplot(fj.ani %>% filter(ANI>0), aes(x=QueryAligned))+
  geom_histogram(binwidth = 200000, size=.25, color='white', fill='gray30')+
  theme_classic()+scale_x_continuous(breaks=seq(0,10000000,500000))

ggplot(fj.ani %>% filter(ANI>0), aes(x=ANI))+
  geom_histogram(size=.25, color='white', fill='gray30')+
  theme_classic()

leaveout <- fj.ani %>% filter(ANI==0 | QueryAligned <= aligncut)

aligncut <- 1 * 1000000 ## How much should be aligned
fj <- data.frame(Reference=89329, ANI=100)
fj.76 <- fj.ani %>% filter(ANI>76 & QueryAligned > aligncut) %>% select(Reference, ANI) %>% rbind(fj)
fj.80 <- fj.ani %>% filter(ANI>80 & QueryAligned > aligncut) %>% select(Reference, ANI) %>% rbind(fj)
fj.95 <- fj.ani %>% filter(ANI>95 & QueryAligned > aligncut) %>% select(Reference, ANI) %>% rbind(fj)

i <- rbind(fj, fj.ani %>% select(Reference, ANI) %>% filter(!Reference %in% leaveout$Reference))
i$in100 <- 'No'
i[i$Reference %in% fj$Reference,]$in100 <- 'Yes'
i[i$Reference %in% fj$Reference,]$in100 <- 'FJ'
i$in95 <- 'No'
i[i$Reference %in% fj.95$Reference,]$in95 <- 'Yes'
i[i$Reference %in% fj$Reference,]$in95 <- 'FJ'
i$in80 <- 'No'
i[i$Reference %in% fj.80$Reference,]$in80 <- 'Yes'
i[i$Reference %in% fj$Reference,]$in80 <- 'FJ'
i$in76 <- 'No'
i[i$Reference %in% fj.76$Reference,]$in76 <- 'Yes'
i[i$Reference %in% fj$Reference,]$in76 <- 'FJ'
write.table(i, 'in_groups.tsv', sep="\t", quote=F, row.names = F)

straintax <- read.table("flavo_stats_cut.tsv", header=F, col.names = c('Reference', 'contigs', 'len', 'lenmb', 'n50', 'l50', 'n90', 'gc', 'sp'), sep="\t") %>% select(Reference, sp)
ani.st <- left_join(fj.ani, straintax, by='Reference')
#spcount <- ani.st %>% group_by(sp) %>% summarize(count=n(), avgani=mean(ANI)) %>% arrange(-count) %>% filter(count > 2)
#f <- ani.st %>% filter(grepl('^Flavobacterium', sp)) %>% filter(sp %in% spcount$sp)
f <- ani.st %>% filter(grepl('^Flavobacterium johnsoniae', sp))
of <- ani.st %>% filter(!(sp %in% f$sp))
of$sp <- 'Other Flavobacteriaceae spp.'
ani.st.ss <- rbind(f, of)

ggplot(ani.st.ss %>% filter(!(Reference %in% leaveout$Reference)), aes(x=reorder(sp, ANI), y=ANI))+
  geom_hline(yintercept=80, color='#5a86c4', linetype='dashed')+
  geom_hline(yintercept=95, color='#af5ac4', linetype='dashed')+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=.3)+
  coord_flip()+
  theme_classic()+
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_text(face='italic')
  )+
  ylab("ANI to Fj UW101")
ggsave("fjani_by_sp.pdf", dpi=600, height=3, width=7)


ggplot(ani.st.ss %>% filter(!(Reference %in% leaveout$Reference)), aes(y=reorder(sp, ANI), x=ANI))+
  geom_vline(xintercept=80, color='#5a86c4', linetype='dashed')+
  geom_vline(xintercept=95, color='#af5ac4', linetype='dashed')+
  geom_density_ridges(quantile_lines=T, quantiles=2, alpha=.7, point_alpha=.3, jittered_points=T, point_size=.4, position='raincloud', scale=.7)+
  theme_classic()+theme(
    legend.title = element_blank(),
    legend.position = 'top',
    axis.text.y = element_text(face='italic'),
    axis.title.y = element_blank()
  )+
  xlab("ANI to Fj UW101")+
  scale_x_continuous(breaks=seq(75,105,5))
ggsave("fjani_by_sp_rain.pdf", dpi=600, height=3, width=7)


ggplot(fj.ani %>% filter(ANI>0), aes(x=ANI))+
  geom_vline(xintercept=80, color='#5a86c4', linetype='dashed')+
  geom_vline(xintercept=95, color='#af5ac4', linetype='dashed')+
  geom_histogram(binwidth = .75, size=.25, color='white', fill='gray30')+
  theme_classic()+scale_x_continuous(breaks=seq(75,100,2.5))+
  ylab("Number of Genomes")+xlab("ANI to Fj UW101")
ggsave("fjani_hist.pdf", dpi=600, height=3, width=5)


### Break

mibig.groups <- read.table("mibig_groups.tsv", header=T, sep="\t", quote="\"", comment.char = "")
hm.t <- read.table("ortho_tall.tsv", header=T, sep="\t", quote="\"", comment.char = "")
hm.t$in100 <- 'No'
hm.t[hm.t$strain_id %in% fj$Reference,]$in100 <- 'Yes'
hm.t$in76 <- 'No'
hm.t[hm.t$strain_id %in% fj.76$Reference,]$in76 <- 'Yes'
hm.t$in80 <- 'No'
hm.t[hm.t$strain_id %in% fj.80$Reference,]$in80 <- 'Yes'
hm.t$in95 <- 'No'
hm.t[hm.t$strain_id %in% fj.95$Reference,]$in95 <- 'Yes'

## 100
by100 <- hm.t %>% group_by(ortho_group, in100) %>% summarize(n=n(), avg_ortho=mean(n_ortho), sd_ortho=sd(n_ortho), med_ortho=median(n_ortho))
n <- nrow(by100[is.na(by100$sd_ortho),])
if(n>0){
 by100[is.na(by100$sd_ortho),]$sd_ortho <- 0
}
by100n <- by100 %>% filter(in100=='No')
colnames(by100n) <- c('ortho_group', 'in100', 'no_n', 'no_avg', 'no_sd', 'no_med')
by100n$in100 <- NULL
by100y <- by100 %>% filter(in100=='Yes')
colnames(by100y) <- c('ortho_group', 'in100', 'yes_n', 'yes_avg', 'yes_sd', 'yes_med')
by100y$in100 <- NULL
by100 <- left_join(by100y, by100n, by='ortho_group')
by100$avg_expansion <- by100$yes_avg-by100$no_avg
by100$med_expansion <- by100$yes_med-by100$no_med
by100$med_log2fc <- log2(by100$yes_med) - log2(by100$no_med)
by100$avg_log2fc <- log2(by100$yes_avg) - log2(by100$no_avg)

## 80
by80 <- hm.t %>% group_by(ortho_group, in80) %>% summarize(n=n(), avg_ortho=mean(n_ortho), sd_ortho=sd(n_ortho), med_ortho=median(n_ortho))
n <- nrow(by80[is.na(by80$sd_ortho),])
if(n>0){
  by80[is.na(by80$sd_ortho),]$sd_ortho <- 0
}
by80n <- by80 %>% filter(in80=='No')
colnames(by80n) <- c('ortho_group', 'in80', 'no_n', 'no_avg', 'no_sd', 'no_med')
by80n$in80 <- NULL
by80y <- by80 %>% filter(in80=='Yes')
colnames(by80y) <- c('ortho_group', 'in80', 'yes_n', 'yes_avg', 'yes_sd', 'yes_med')
by80y$in80 <- NULL
by80 <- left_join(by80y, by80n, by='ortho_group')
by80$avg_expansion <- by80$yes_avg-by80$no_avg
by80$med_expansion <- by80$yes_med-by80$no_med
by80$med_log2fc <- log2(by80$yes_med) - log2(by80$no_med)
by80$avg_log2fc <- log2(by80$yes_avg) - log2(by80$no_avg)

## 95
by95 <- hm.t %>% group_by(ortho_group, in95) %>% summarize(n=n(), avg_ortho=mean(n_ortho), sd_ortho=sd(n_ortho), med_ortho=median(n_ortho))
n <- nrow(by95[is.na(by95$sd_ortho),])
if(n>0){
  by95[is.na(by95$sd_ortho),]$sd_ortho <- 0
}
by95n <- by95 %>% filter(in95=='No')
colnames(by95n) <- c('ortho_group', 'in95', 'no_n', 'no_avg', 'no_sd', 'no_med')
by95n$in95 <- NULL
by95y <- by95 %>% filter(in95=='Yes')
colnames(by95y) <- c('ortho_group', 'in95', 'yes_n', 'yes_avg', 'yes_sd', 'yes_med')
by95y$in95 <- NULL
by95 <- left_join(by95y, by95n, by='ortho_group')
by95$avg_expansion <- by95$yes_avg-by95$no_avg
by95$med_expansion <- by95$yes_med-by95$no_med
by95$med_log2fc <- log2(by95$yes_med) - log2(by95$no_med)
by95$avg_log2fc <- log2(by95$yes_avg) - log2(by95$no_avg)

## Function stuff
fj.func <- read.table("89329.pep.fa.kofam_parsed.tsv", header=T, sep="\t")
#fj.func$Gene <- sub("^89329\\|", "", fj.func$Gene)
keggmap <- read.table("kegg_heir.tsv", sep="\t", header=T, quote="")
fj.func <- fj.func %>% left_join(keggmap, by='KO.D')
g2g <- read.table("group_list.tsv", header=T, sep="\t") %>% filter(strain_id %in% fj)
fj.func <- fj.func %>% left_join(g2g %>% select(Gene, ortho_group), by='Gene')
write.table(fj.func, "fjfunc.tsv", sep="\t", row.names = F)

met <- fj.func %>% filter(DescA=='Metabolism') %>% group_by(Gene) %>% distinct(Gene, .keep_all = T)
aa <- fj.func %>% filter(DescB %in% c('Amino acid metabolism', 'Metabolism of other amino acids')) %>% group_by(Gene) %>% distinct(Gene, .keep_all = T)
vit <- fj.func %>% filter(DescB %in% c('Metabolism of cofactors and vitamins')) %>% group_by(Gene) %>% distinct(Gene, .keep_all = T)
aavit.keep <- rbind(aa, vit) %>% group_by(ortho_group) %>% distinct(ortho_group, .keep_all = T)
met.keep <- met %>% group_by(ortho_group) %>% distinct(ortho_group)

instrain <- by100 %>% filter(yes_med>0)

by100$group <- 'All ortholog groups'
by100m <- by100 %>% filter(ortho_group %in% met.keep$ortho_group)
by100m$group <- 'Metabolism'
by100a <- by100 %>% filter(ortho_group %in% aavit.keep$ortho_group)
by100a$group <- 'Metabolism - AA & Vit. Biosynthesis'
by100mi <- by100 %>% filter(ortho_group %in% mibig.groups$ortho_group) %>% filter(ortho_group %in% instrain$ortho_group)
by100mi$group <- 'MIBiG'
all100 <- rbind(by100m, by100a, by100mi)

by95$group <- 'All ortholog groups'
by95m <- by95 %>% filter(ortho_group %in% met.keep$ortho_group)
by95m$group <- 'Metabolism'
by95a <- by95 %>% filter(ortho_group %in% aavit.keep$ortho_group)
by95a$group <- 'Metabolism - AA & Vit. Biosynthesis'
by95mi <- by95 %>% filter(ortho_group %in% mibig.groups$ortho_group) %>% filter(ortho_group %in% instrain$ortho_group)
by95mi$group <- 'MIBiG'
all95 <- rbind(by95m, by95a, by95mi)

by80$group <- 'All ortholog groups'
by80m <- by80 %>% filter(ortho_group %in% met.keep$ortho_group)
by80m$group <- 'Metabolism'
by80a <- by80 %>% filter(ortho_group %in% aavit.keep$ortho_group)
by80a$group <- 'Metabolism - AA & Vit. Biosynthesis'
by80mi <- by80 %>% filter(ortho_group %in% mibig.groups$ortho_group) %>% filter(ortho_group %in% instrain$ortho_group)
by80mi$group <- 'MIBiG'
all80 <- rbind(by80m, by80a, by80mi)

a100m <- ggplot(all100, aes(x=avg_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  theme_classic()+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (mean")+ylab("Orthologs")
a95m <- ggplot(all95, aes(x=avg_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  theme_classic()+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (mean)")+ylab("Orthologs")
a80m <- ggplot(all80, aes(x=avg_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  theme_classic()+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (mean)")+ylab("Orthologs")

plot_grid(
  a100m, a95m, a80m,
  ncol=1, labels=LETTERS[1:3]
)
ggsave("exp_mean.pdf", dpi=600, height=8, width=8)

a100me <- ggplot(all100, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")
a95me <- ggplot(all95, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")
a80me <- ggplot(all80, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, aes(fill=group), size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,1000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")

plot_grid(
  a100me, a95me, a80me,
  ncol=1, labels=LETTERS[1:4]
)
ggsave("exp_median.pdf", dpi=600, height=8, width=8)

med100 <- ggplot(by100, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, fill='gray30', size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,16000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")
med95 <- ggplot(by95, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, fill='gray30', size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,16000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")
med80 <- ggplot(by80, aes(x=med_expansion))+
  geom_histogram(binwidth = 1, fill='gray30', size=.25, color='white')+
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-.25)+
  scale_y_continuous(limits=c(0,16000))+
  scale_x_continuous(limits=c(-1.5,9.5), breaks=seq(-1,9,1))+
  theme_classic()+
  facet_wrap(~group, nrow = 1)+theme(
    panel.grid.major = element_line(color='gray80', linetype = 'dotted'),
    legend.position = 'none'
  )+scale_fill_discrete_qualitative(palette='Dynamic')+
  xlab("Expansion (median)")+ylab("Orthologs")

plot_grid(
  med100, med95, med80,
  ncol=1, labels=LETTERS[1:3]
)
ggsave("exp_median_all.pdf", dpi=600, height=8, width=8)

huge_100_expansions <- by100 %>% filter(med_expansion >6) %>% select(ortho_group, med_expansion) %>% left_join(fj.func, by='ortho_group')

#write.table(all100, "expansions_100_met.tsv", row.names = F, sep="\t")
#write.table(all95, "expansions_95_met.tsv", row.names = F, sep="\t")
# write.table(all80, "expansions_80_met.tsv", row.names = F, sep="\t")

func1 <- fj.func %>% group_by(ortho_group) %>% distinct(ortho_group, .keep_all = T)

e <- rbind(
  by100 %>% mutate(cut='med100'),
  by95 %>% mutate(cut='med095'),
  by80 %>% mutate(cut='med080')
) %>% dcast(ortho_group ~ cut, value.var = 'med_expansion') %>% filter(med100>0)
e$Metabolism <- 'No'
e[e$ortho_group %in% met.keep$ortho_group,'Metabolism'] <- 'Yes'
e$AAVit <- 'No'
e[e$ortho_group %in% aavit.keep$ortho_group,'AAVit'] <- 'Yes'
e$MIBiG <- 'No'
e[e$ortho_group %in% mibig.groups$ortho_group,'MIBiG'] <- 'Yes'
e <- left_join(e, func1 %>% select(ortho_group, DescD), by='ortho_group')
e$DescD <- as.character(e$DescD)
e[is.na(e$DescD),'DescD'] <- 'No KO annotation'
colnames(e)[8] <- 'Annotation'
write.table(e, "med_expansions.tsv", row.names = F, sep="\t")
  
