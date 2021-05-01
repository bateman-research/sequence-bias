# Analysis of DNA composition for a subset of genes or genomes
# Report the nucleotide fraction, GC content and asymmetries
# 
# Aleix Lafita - October 2019

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(gridExtra))

theme_set(theme_bw() + theme(panel.grid = element_blank()))


###################### Argparse #############################

input = "examples/periscopes_dna.fa"
prefix = "examples/periscopes_dna"

# create parser object
parser = ArgumentParser(description = 'Analysis of DNA composition for a subset of genes or genomes')

# specify our desired options 
parser$add_argument("-d", "--dna", default=input,
                    help="DNA sequences in FASTA format [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

dna = args$dna
prefix = args$prefix

###################### Parsing #############################

# Parse the gene sequences and convert to DF
genes.fa = read.fasta(dna)
genes.gc = data.frame(
  acc=names(genes.fa), 
  dna=unlist(getSequence(genes.fa, as.string=T)),
  stringsAsFactors = F
) %>% rowwise() %>%mutate(
  dna = toupper(dna),
  codon1 = paste0(unlist(strsplit(dna, ""))[seq(1,nchar(dna),3)], collapse = ""),
  codon2 = paste0(unlist(strsplit(dna, ""))[seq(2,nchar(dna),3)], collapse = ""),
  codon3 = paste0(unlist(strsplit(dna, ""))[seq(3,nchar(dna),3)], collapse = ""),
)
genes.fa = NULL

###################### Calculate #############################

nts = "ACTG"
nt = strsplit(nts, "")[[1]]

genes.gc = genes.gc %>%
  mutate(bias = 0)

# Routine to calculate the amino acid frequencies for each sequence
for (n in nt) {
  for (f in 1:3) {
    nf = paste0("n", n, f, "fr")
    nc = paste0("n", n, f)
    c = paste0("codon", f)
    genes.gc = genes.gc %>%
      mutate(
        !!nc := str_count(!!sym(c), n),
        !!nf := !!sym(nc) / str_count(!!sym(c))
        #bias = bias - ifelse(eval(as.symbol(nf)) > 0, eval(as.symbol(nf)) * log2(eval(as.symbol(nf))), 0)
      )
  }
}

# Calculate the total nucleotide frequencies
genes.gc = genes.gc %>% 
  mutate(
    nA = nA1 + nA2 + nA3,
    nT = nT1 + nT2 + nT3,
    nG = nG1 + nG2 + nG3,
    nC = nC1 + nC2 + nC3,
    nAfr = nA / nchar(dna),
    nTfr = nT / nchar(dna),
    nGfr = nG / nchar(dna),
    nCfr = nC / nchar(dna)
  )

# Calculate the global nucleotide frequencies (null model)
nt_freq = colSums(genes.gc %>% select(paste0("n", nt)) / sum(genes.gc %>% select(paste0("n", nt))))

# Calculate the nucleotide bias
for (n in nt) {
  nf = paste(c("n", n, "fr"), collapse = "")
  nb = paste(c("bias", n), collapse = "")
  genes.gc = genes.gc %>%
    rowwise() %>%
    mutate(
      !!nb := eval(as.symbol(nf)) * log2(eval(as.symbol(nf)) / nt_freq[paste0("n", n)]),
      bias = bias + ifelse(eval(as.symbol(nf)) > 0, eval(as.symbol(nb)), 0)
    )
}

# Calculate some other properties derived from nt fractions
genes.gc = genes.gc %>%
  rowwise() %>%
  mutate(
    nTotal = nA + nT + nC + nG,
    nGCfr = nGfr + nCfr,
    aAT = nAfr - nTfr,
    sAT = aAT / max((nAfr + nTfr), 0.00000001),
    aGC = nGfr - nCfr,
    sGC = aGC / max((nGfr + nCfr), 0.00000001),
    aPurPyr = (nAfr + nGfr) - (nCfr + nTfr)
  )


###################### Save results #############################

write.table(
  genes.gc %>% select(-dna),
  sprintf("%s_composition.tsv", prefix),
  quote = F,
  sep = "\t",
  row.names = F
)

###################### Plots #############################

# Plot the distribution of GC content
p = ggplot(genes.gc, aes(x = nGCfr)) +
  geom_step(stat = "bin", boundary = 0, binwidth = 0.01, color = "black") +
  xlab("Fraction of GC") +
  xlim(0,1)

pdf(sprintf("%s_gc-fraction.pdf", prefix), 4, 3)
plot(p)
log = dev.off()


# Plot the distribution of sequence bias
p = ggplot(genes.gc, aes(x = bias)) +
  geom_step(stat = "bin", boundary = 0, binwidth = 0.01, color = "black") +
  #xlim(0,max(genes.gc$bias)) +
  xlab("Relative entropy")

pdf(sprintf("%s_nt-bias.pdf", prefix), 4, 3)
plot(p)
log = dev.off()

# Gather the amino acid compositions into a column
genes.ntcomp = gather(
  genes.gc %>% select(nAfr, nTfr, nCfr, nGfr, acc),
  key, fraction, -acc
) %>% mutate(nt = substring(key, 2, 2))

# Plot the distribution of each nucleotide fraction
p = ggplot(genes.ntcomp, aes(x = nt, y = fraction)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0, scale = "width") + 
  geom_hline(yintercept = 0.25, alpha = 0.3) +
  ylab("Nucleotide fraction") +
  ylim(0, max(genes.ntcomp$fraction) + 0.05) +
  xlab("")

pdf(sprintf("%s_nt-fraction.pdf", prefix), 6, 6)
plot(p)
log = dev.off()

# Gather compositon for each codon 
genes.codons = genes.gc %>%
  select(acc, nA1fr, nA2fr, nA3fr,
         nT1fr, nT2fr, nT3fr,
         nG1fr, nG2fr, nG3fr,
         nC1fr, nC2fr, nC3fr) %>%
  gather(key, fraction, -acc) %>%
  mutate(
    nt = substring(key, 2, 2),
    frame = substring(key, 3, 3),
  )

# Plot the distribution of each nucleotide fraction
p = ggplot(genes.codons, aes(x = frame, y = fraction)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0, scale = "width") + 
  facet_grid(~nt) +
  geom_hline(yintercept = 0.25, alpha = 0.3) +
  ylab("Nucleotide fraction") +
  ylim(0, max(genes.codons$fraction) + 0.05) +
  xlab("Codon position")

pdf(sprintf("%s_nt-codons.pdf", prefix), 10, 5)
plot(p)
log = dev.off()


# Gather the other properties into a column
genes.ntasymm = gather(
  genes.gc %>% select(aAT, sAT, aGC, sGC, aPurPyr, acc),
  nt, fraction, -acc
)

# Plot the distribution of nucleotide asymmetry and skew
p = ggplot(genes.ntasymm, aes(x = nt, y = fraction)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0, scale = "width") + 
  geom_hline(yintercept = 0, alpha = 0.5) +
  xlab("") +
  ylab("Nucleotide asymmetry (X-Y) and skew (X-Y)/(X+Y)") +
  ylim(-1, 1)
  

pdf(sprintf("%s_nt-asymmetry.pdf", prefix), 6, 6)
plot(p)
log = dev.off()


###################### PCA #############################

genes.pca = genes.gc %>% select(nAfr, nTfr, nGfr, nCfr)
names(genes.pca) = c("A", "T", "G", "C")

genes.pca.calc = prcomp(genes.pca, center = TRUE, scale. = TRUE)

#pdf(sprintf("%s_pca_var.pdf", prefix), 5, 4)
#plot(genes.pca.calc, type = 'l')
#dev.off()

p = autoplot(
    genes.pca.calc,
    data = genes.pca.calc,
    #colour = 'type',
    loadings = TRUE,
    loadings.colour = 'black',
    loadings.label = TRUE,
    loadings.label.size = 4,
    loadings.label.colour = "black",
    loadings.label.vjust = 2,
    size = 1,
    alpha = 0.3
  )

pdf(sprintf("%s_nt-pca.pdf", prefix), 5, 5)
print(p)
log = dev.off()


###################### Dendrogram #############################

genes.matrix = genes.pca
rownames(genes.matrix) = genes.gc$acc
genes.dendro = as.dendrogram(hclust(d = dist(x = genes.matrix)))

gene.order = order.dendrogram(genes.dendro)

p1 = ggplot(data = genes.ntcomp, aes(x = gsub("n", "", gsub("fr", "", nt)), y = factor(acc, levels = genes.gc$acc[gene.order]))) +
  geom_tile(aes(fill = fraction)) +
  scale_fill_gradient(name = "Fraction", low = "white", high = "black") +
  theme(
    panel.border = element_blank(),
    #axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_blank(),
    legend.position = "left"
  ) + scale_x_discrete(expand = c(0, 0))
  

p2 = ggdendrogram(genes.dendro, rotate = FALSE, size = 1) +
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  ) + coord_flip() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.005, 0.005))

# Get the ggplot grobs
gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)

# Disable for now
#pdf(sprintf("%s_nt-dendro.pdf", prefix), 10, 10 + nrow(genes.gc) / 50)
# Combine the two plots
#grid.arrange(gp1, gp2, ncol = 2)
#log = dev.off()


###################### Expected aa composition #############################

# Create all possible codons
nt.double = merge(nt, nt)
nt.triple = merge(nt.double, nt, by =c()) 
names(nt.triple) = c("x", "y", "z")

codons = nt.triple %>%
  mutate(codon = paste0(x, y, z)) %>%
  gather(pos, nt, -codon)

# Create the matrix of all codons per gene and calculate theoretical fraction
genes.codons = genes.matrix %>%
  ungroup() %>%
  mutate(acc = rownames(genes.matrix)) %>%
  gather(nt, frac, -acc) %>%
  merge(codons) %>%
  group_by(codon, acc) %>%
  summarise(
    prob = prod(frac)
  ) %>%
  mutate(aa = translate(unlist(strsplit(codon, "")))) %>%
  group_by(acc, aa) %>%
  summarise(prob = sum(prob))

genes.codons.table = genes.codons %>%
  filter(aa != "*") %>%
  spread(aa, prob)

write.table(
  genes.codons.table,
  sprintf("%s_expected-aa.tsv", prefix),
  quote = F,
  sep = "\t",
  row.names = F
)

# Plot the distribution of each nucleotide fraction
p = ggplot(genes.codons %>% filter(aa != "*"), aes(x = aa, y = prob)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0, scale = "width") + 
  ylab("Expected AA fraction") +
  xlab("") +
  ylim(-0.01, max(genes.codons$prob) + 0.01)

pdf(sprintf("%s_expected-aa.pdf", prefix), 10, 6)
plot(p)
log = dev.off()





