# Amino acid composition analysis for a dataset of proteins
# Aleix Lafita - January 2020

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

input = "examples/periscopes_prot.fa"
prefix = "examples/periscopes_prot.fa"
minlen = 30

# create parser object
parser = ArgumentParser(description = 'Amino acid composition analysis for a subset of proteins')

# specify our desired options 
parser$add_argument("-s", "--seqs", default=input,
                    help="Protein sequences in FASTA format [default \"%(default)s\"]")
parser$add_argument("-l", "--minlen", default=minlen,
                    help="Minimum length of protein sequence; filter out shorter sequences [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

proteins = args$seqs
prefix = args$prefix
minlen = as.integer(args$minlen)

###################### Parsing #############################

# Parse the gene sequences and convert to DF
proteins.fa = read.fasta(proteins)
proteins.df = data.frame(
  acc=names(proteins.fa), 
  seq=unlist(getSequence(proteins.fa, as.string=T)),
  stringsAsFactors = F
) %>% mutate(
  seq = toupper(seq),
  seqlen = nchar(seq)
) %>% filter(seqlen > minlen)
proteins.fa = NULL

###################### Sequence bias #############################

aas = "ACDEFGHIKLMNPQRSTVWY"
aa = strsplit(aas, "")[[1]]

prot.comp = proteins.df %>%
  mutate(bias = 0, bias.u = 0)

# Matrix of amino acids along the protein
for (a in aa) {
  af = paste(c(a, "fr"), collapse = "")
  prot.comp = prot.comp %>%
    mutate(
      !!a := str_count(seq, a),
      !!af := str_count(seq, a) / str_count(seq)
    )
}

# Calculate the global amino acid frequencies (null model)
aa_freq = colSums(prot.comp %>% select(aa)) / sum(prot.comp %>% select(aa))

# Calculate sequence bias as a relative entropy
for (a in aa) {
  af = paste(c(a, "fr"), collapse = "")
  ab = paste(c("bias", a), collapse = "")
  au = paste(c("ubias", a), collapse = "") # This is the uniform background bias
  prot.comp = prot.comp %>%
    rowwise() %>%
    mutate(
      !!ab := eval(as.symbol(af)) * log2(eval(as.symbol(af)) / aa_freq[a]),
      !!au := eval(as.symbol(af)) * log2(eval(as.symbol(af)) / 0.05),
      bias.u = bias.u + ifelse(!is.na(eval(as.symbol(ab))), eval(as.symbol(ab)), 0),
      bias = bias + ifelse(!is.na(eval(as.symbol(au))), eval(as.symbol(au)), 0),
    )
}

# Save table
write.table(
  prot.comp %>% select(-seq),
  sprintf("%s_composition.tsv", prefix),
  quote = F,
  sep = "\t",
  row.names = F
)

# Plot the distribution of sequence bias
p = ggplot(prot.comp, aes(x = bias)) +
  geom_step(stat = "bin", boundary = 0, bins = 20, color = "black") +
  xlim(0, max(prot.comp$bias) + 0.05) +
  xlab("Sequence bias (Relative entropy)")

pdf(sprintf("%s_aa-bias.pdf", prefix), 5, 4)
plot(p)
log = dev.off()

# Gather the amino acid compositions into a column
proteins.aacomp = gather(
  prot.comp %>% select(paste0(aa, "fr"), acc),
  aa, fraction, -acc
) %>% mutate(
  aa = gsub("fr", "", aa)
)

# Plot the distribution of each nucleotide fraction
p = ggplot(proteins.aacomp, aes(x = aa, y = fraction)) +
  geom_jitter(alpha = 0.1) +
  geom_violin(alpha = 0, scale = "width") + 
  ylab("Amino acid fraction") +
  xlab("") +
  ylim(-0.01, max(proteins.aacomp$fraction) + 0.01)

pdf(sprintf("%s_aa-fraction.pdf", prefix), 10, 6)
plot(p)
log = dev.off()


###################### PCA #############################

proteins.pca = prot.comp %>% select(paste0(aa, "fr"))
names(proteins.pca) = aa

proteins.pca.calc = prcomp(proteins.pca, center = TRUE, scale. = TRUE)

#pdf(sprintf("%s_pca_var.pdf", prefix), 5, 4)
#plot(proteins.pca.calc, type = 'l')
#dev.off()

p = autoplot(
  proteins.pca.calc,
  data = proteins.pca.calc,
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

pdf(sprintf("%s_aa-pca.pdf", prefix), 5, 5)
print(p)
log = dev.off()


###################### Dendrogram #############################

proteins.matrix = proteins.pca
rownames(proteins.matrix) = prot.comp$acc
proteins.dendro = as.dendrogram(hclust(d = dist(x = proteins.matrix)))

protein.order = order.dendrogram(proteins.dendro)

p1 = ggplot(data = proteins.aacomp, aes(x = aa, y = factor(acc, levels = prot.comp$acc[protein.order]))) +
  geom_tile(aes(fill = fraction)) +
  scale_fill_gradient(name = "Fraction", low = "white", high = "black") +
  theme(
    panel.border = element_blank(),
    #axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #axis.ticks.y = element_blank(),
    legend.position = "left"
  ) + scale_x_discrete(name = "", expand = c(0, 0))


p2 = ggdendrogram(proteins.dendro, rotate = FALSE, size = 1) +
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  ) + coord_flip() + 
  scale_y_continuous(name = "", expand = c(0, 0)) +
  scale_x_continuous(name = "", expand = c(0.005, 0.005))

# Get the ggplot grobs
gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)

pdf(sprintf("%s_aa-dendro.pdf", prefix), 15, 10 + nrow(prot.comp) / 50)
# Combine the two plots
grid.arrange(gp1, gp2, ncol = 2)
log = dev.off()

# Could do a dendogram but with enrichment compared to average

