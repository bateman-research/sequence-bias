# Analyze the nucleotide bias along the DNA sequence of a gene
# Aleix Lafita - October 2019

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(zoo))

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(stringr))

suppressPackageStartupMessages(library(ggplot2))

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
options(stringsAsFactors = F)

###################### Argparse #############################

input = "example/SasG_dna.fa"
prefix = "example/SasG_dna"
N = 100
frame = "1,2,3"

# create parser object
parser = ArgumentParser(description = 'Nucleotide bias and composition along the DNA sequence of a gene')

# specify our desired options 
parser$add_argument("-s", "--seq", default=input,
                    help="DNA sequence in FASTA format [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output files, results, plots and tables [default \"%(default)s\"]")
parser$add_argument("-N", "--rollN", default=N,
                    help="Number of nucleotides to use for rolling mean calculation [default \"%(default)s\"]")
parser$add_argument("-f", "--frame", default=frame,
                    help="Codon positions (1,2,3) to use, comma separated [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

input = args$seq
prefix = args$prefix
N = as.integer(args$rollN)

frame = as.integer(unlist(strsplit(args$frame, ",")))
if (any(frame > 3 | frame < 1)) {
  stop(sprintf("Frame numbers should be in {1,2,3} only was: %s", args$frame))
}

############################# Parsing ####################################
message(sprintf("# Parsing sequence '%s'", input))

# Parse the DNA sequences and convert to DF
dna.fasta = read.fasta(input)[1]
dnalen = length(getSequence(dna.fasta)[[1]])

dna.df = data.frame(
  id = names(dna.fasta),
  pos = seq(1, dnalen, 1),
  seq = unlist(getSequence(dna.fasta))
) %>% mutate(seq = toupper(seq)) %>%
  filter(is.element(pos %% 3, frame %% 3))

# Extract the protein sequence too, to connect DNA and protein
prot.df = data.frame(
  pos = c(seq(1, dnalen, 3),seq(2, dnalen, 3),seq(3, dnalen, 3)),
  aaseq = translate(unlist(getSequence(dna.fasta)))
) %>% mutate(aapos = 1 + floor((pos-1) / 3))

dna.df = merge(dna.df, prot.df)

seqlen = nrow(dna.df)

########################### Calculate ###################################
message("# Calculating rolling sequence properties...")

# Create matrix of nucleotides
dna.df = dna.df %>%
  mutate(
    nA = str_count(seq, 'A'),
    nC = str_count(seq, 'C'),
    nT = str_count(seq, 'T'),
    nG = str_count(seq, 'G'),
  )

# Calculate the rolling averages of each nucleotide
dna.df$rollA = c(rollmean(dna.df$nA, N), rep(mean(dna.df$nA[(seqlen-N):seqlen]), N-1))
dna.df$rollC = c(rollmean(dna.df$nC, N), rep(mean(dna.df$nC[(seqlen-N):seqlen]), N-1))
dna.df$rollT = c(rollmean(dna.df$nT, N), rep(mean(dna.df$nT[(seqlen-N):seqlen]), N-1))
dna.df$rollG = c(rollmean(dna.df$nG, N), rep(mean(dna.df$nG[(seqlen-N):seqlen]), N-1))

# Calculate cumsums of each nucleotide
dna.df$cumA = cumsum(dna.df$nA)
dna.df$cumC = cumsum(dna.df$nC)
dna.df$cumT = cumsum(dna.df$nT)
dna.df$cumG = cumsum(dna.df$nG)

# Create matrix of amino acids
aas = "ACDEFGHIKLMNPQRSTVWY"
aa = strsplit(aas, "")[[1]]

dna.df = dna.df %>%
  mutate(bias = 0)

# Matrix of amino acids along the protein
for (a in aa) {
  dna.df = dna.df %>%
    mutate(!!a := str_count(aaseq, a))
}

# Calculate the rolling means of each amino acid and sequence bias
for (a in aa) {
  rolla = paste0("r", a)
  dna.df[rolla] = 
    c(rollmean(dna.df[,a], N), rep(mean(dna.df[dna.df$pos > (seqlen-N),][,a]), N-1))
  dna.df = dna.df %>%
    mutate(bias = bias + ifelse(eval(as.symbol(rolla)) == 0, 0, eval(as.symbol(rolla)) *log2(eval(as.symbol(rolla)) / 0.05)))
}

# Calculate the sequence bias
dna.bias = dna.df %>%
  #rowwise() %>%
  #filter(!is.na(rollA)) %>%
  mutate(
    bias = ifelse(rollA == 0, 0, rollA *log2(rollA / 0.25)) +
      ifelse(rollC == 0, 0, rollC *log2(rollC / 0.25)) + 
      ifelse(rollT == 0, 0, rollT *log2(rollT / 0.25)) + 
      ifelse(rollG == 0, 0, rollG *log2(rollG / 0.25)),
    rollGC = rollG + rollC,
    rollAT = rollA + rollT,
    rollATasym = rollA - rollT,
    rollGCasym = rollG - rollC,
    rollATskew = ifelse(rollAT == 0, 0, rollATasym / rollAT),
    rollGCskew = ifelse(rollGC == 0, 0, rollGCasym / rollGC),
    cumATskew.tot = (cumA - cumT),
    cumGCskew.tot = (cumG - cumC),
    cumATskew = cumATskew.tot / max(abs(cumATskew.tot)),
    cumGCskew = cumGCskew.tot / max(abs(cumGCskew.tot))
  )

step = ceiling(N/100)

dna.subsample = dna.bias[seq(1, nrow(dna.bias), step),]

# Gather the nucleotide composition into a key-value
dna.rollcomp = dna.subsample %>% 
  select(pos, rollA, rollC, rollT, rollG) %>%
  gather(key, value, -pos) %>%
  mutate(
    key = gsub("roll", "", key)
  )

dna.rollasym = dna.subsample %>%
  select(pos, rollATasym, rollGCasym) %>%
  gather(key, value, -pos)

dna.rollskew = dna.subsample %>%
  select(pos, rollATskew, rollGCskew, cumATskew, cumGCskew) %>%
  gather(key, value, -pos)

# Correlation of sequence bias with amino acid fraction
dna.rollskewnt = dna.bias %>% 
  select(pos, seq, rollATskew, rollA, rollT, rollG, rollC) %>%
  gather(nt, fraction, -pos, -seq, -rollATskew) %>%
  mutate(
    nt = gsub("roll", "", nt)
  )

# Correlation of AT skew with amino acid fraction
dna.rollskewaa = dna.bias %>% 
  select(pos, rollATasym, paste0("r", aa)) %>%
  gather(aa, fraction, -pos, -rollATasym) %>%
  mutate(
    aa = gsub("r", "", aa)
  )

  
########################### Plots ###################################
message("# Generating figures...")

# Plot width
width = max(min((dnalen)/(5*N), 8), 5)

# Plot the rolling composition of amino acids
p = ggplot(dna.rollcomp, aes(x = pos, y = value, color = key)) +
  geom_hline(yintercept = 0.25, alpha = 0.5) +
  geom_step() +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  #facet_grid(~frame) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, dnalen)) +
  xlab("Sequence position") +
  ylab("Nucleotide frequency")

pdf(sprintf("%s_roll-nt.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the rolling sequence bias
p = ggplot(dna.subsample, aes(x = pos, y = bias)) +
  geom_step() +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, dnalen)) +
  xlab("Sequence position") +
  ylab("Relative entropy")

pdf(sprintf("%s_roll-bias.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the nuclotide asymmetry
p = ggplot(dna.rollasym, aes(x = pos, y = value, color = key)) +
  geom_step() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, dnalen)) +
  xlab("Sequence position") +
  ylab("Nucleotide asymmetry (X-Y)") +
  ylim(-1,1)

pdf(sprintf("%s_roll-asymm.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the nuclotide skew
p = ggplot(dna.rollskew, aes(x = pos, y = value, color = key)) +
  geom_step() +
  geom_hline(yintercept = 0, alpha = 0.5) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, dnalen)) +
  xlab("Sequence position") +
  ylab("Nucleotide skew (X-Y)/(X+Y)") +
  ylim(-1,1)

pdf(sprintf("%s_roll-skew.pdf", prefix), width, 3)
plot(p)
log = dev.off()

# Plot the correlation between sequence bias and AA fraction
p = ggplot(dna.rollskewaa, aes(x = rollATasym, y = fraction)) +
  geom_point(alpha = 0.5, size = 1) +
  #geom_histogram(aes(fill = aa), alpha = 1, color = "black", position = "fill", binwidth = 0.1) +
  #geom_histogram(alpha = 0.1, color = "black", binwidth = 0.1) +
  geom_smooth(method = "lm") +
  facet_wrap(~aa) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, seqlen)) +
  theme(panel.grid = element_blank()) +
  xlab("AT skew") +
  ylab("Fraction")

pdf(sprintf("%s_roll-aaskew.pdf", prefix), 10, 8)
plot(p)
log = dev.off()

# Plot the correlation between sequence bias and NT fraction
p = ggplot(dna.rollskewnt, aes(x = rollATskew, y = fraction)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_smooth(method = "lm") +
  facet_wrap(~nt) +
  #scale_x_continuous(expand = c(0, 0), limits = c(0, seqlen)) +
  theme(panel.grid = element_blank()) +
  xlab("AT skew") +
  ylab("Fraction")

pdf(sprintf("%s_roll-ntskew.pdf", prefix), 5, 5)
plot(p)
log = dev.off()

message("# Done!")
