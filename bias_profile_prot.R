# Amino acid bias and composition profile along a protein sequence
# Aleix Lafita - November 2019

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(zoo))

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))

suppressPackageStartupMessages(library(ggplot2))

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))

###################### Argparse #############################

input = "SasG_prot.fa"
prefix = "SasG_prot"
N = 100

# create parser object
parser = ArgumentParser(description = 'Amino acid bias and composition profile along a protein sequence')

# specify our desired options 
parser$add_argument("-s", "--seq", default=input,
                    help="Protein sequence in FASTA format [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")
parser$add_argument("-N", "--rollN", default=N,
                    help="Number of amino acids to use for rolling mean calculation [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

input = args$seq
prefix = args$prefix
N = as.integer(args$rollN)

############################# Parsing ####################################
message(sprintf("# Parsing sequence '%s'", input))

# Parse the DNA sequences and convert to DF
prot.fasta = read.fasta(input)[1]

seqlen = length(getSequence(prot.fasta)[[1]])

prot.df = data.frame(
  id = names(prot.fasta),
  pos = seq(1, seqlen, 1),
  seq = unlist(getSequence(prot.fasta)),
  stringsAsFactors = F
) %>% mutate(seq = toupper(seq))

prot.all = prot.df[order(prot.df$pos),]

########################### Calculate ###################################
message("# Calculating rolling sequence properties...")

aas = "ACDEFGHIKLMNPQRSTVWY"
aa = strsplit(aas, "")[[1]]

prot.seqbias = prot.all %>%
  mutate(bias = 0)

# Matrix of amino acids along the protein
for (a in aa) {
  prot.seqbias = prot.seqbias %>%
    mutate(!!a := str_count(seq, a))
}

# Calculate the rolling means of each amino acid
for (a in aa) {
  rolla = paste0("roll", a)
  prot.seqbias[rolla] = 
    c(rollmean(prot.seqbias[,a], N), rep(mean(prot.seqbias[prot.seqbias$pos > (seqlen-N),][,a]), N-1))
  prot.seqbias = prot.seqbias %>%
    mutate(bias = bias + ifelse(eval(as.symbol(rolla)) == 0, 0, eval(as.symbol(rolla)) *log2(eval(as.symbol(rolla)) / 0.05)))
}

# Gather the nucleotide composition into a key-value
prot.rollcomp = prot.seqbias %>% 
  select(pos, paste0("roll", aa)) %>%
  gather(key, value, -pos) %>%
  mutate(
    key = gsub("roll", "", key)
  )

# Hydrophobicity and charge
prot.seqbias = prot.seqbias %>%
  mutate(
    charge.aaf = rollK + rollR + rollD + rollE,
    charge.fr = rollK + rollR - rollD - rollE,
    hphos.aaf = rollV+rollL+rollI+rollM+rollF+rollW
  )

# Correlation of sequence bias with amino acid fraction
prot.rollbiasaa = prot.seqbias %>% 
  select(pos, bias, paste0("roll", aa)) %>%
  gather(aa, fraction, -pos, -bias) %>%
  mutate(
    aa = gsub("roll", "", aa)
  )

########################### Plots ###################################
message("# Generating figures...")

# Plot width
width = 6

# Plot the rolling composition
p = ggplot(prot.rollcomp, aes(x = pos, y = key, fill = value)) +
  geom_raster() +
  theme(legend.position = "none") +
  scale_fill_gradient2(name = "Rolling frequency", low = "white", high = "black") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, seqlen)) +
  xlab("Sequence position") +
  ylab("")

pdf(sprintf("%s_roll-aa.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the hydrophobics fraction - disabled for now
p = ggplot(prot.seqbias, aes(x = pos, y = hphos.aaf)) +
  geom_step() +
  ylim(0,1) +
  ylab("Fraction of hydrophobic amino acids")

#pdf(sprintf("%s_hydrophobics.pdf", prefix), width, 3)
#plot(p)
#log = dev.off()

# Plot the charged residues fraction - disabled for now
p = ggplot(prot.seqbias, aes(x = pos, y = charge.aaf)) +
  geom_step() +
  ylim(0,1) + ylab("Fraction of charged amino acids")

#pdf(sprintf("%s_charge.pdf", prefix), width, 3)
#plot(p)
#log = dev.off()

# Plot the rolling sequence bias
p = ggplot(prot.seqbias, aes(x = pos, y = bias)) +
  geom_step() +
  scale_x_continuous(expand = c(0, 0), limits = c(0, seqlen)) +
  xlab("Sequence position") +
  ylab("Relative entropy")

pdf(sprintf("%s_roll-bias.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the correlation between sequence bias and AA fraction
p = ggplot(prot.rollbiasaa, aes(x = bias, y = fraction)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_smooth(method = "lm") +
  facet_wrap(~aa) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, seqlen)) +
  theme(panel.grid = element_blank()) +
  xlab("Relative entropy (sequence bias)") +
  ylab("Fraction")

pdf(sprintf("%s_roll-aabias.pdf", prefix), 10, 8)
plot(p)
log = dev.off()

message("# Done!")

