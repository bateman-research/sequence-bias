# Generate a dot-plot of two sequences
# Not optimized: using 2Gb RAM for 3K sequences - don't go over 5K
# 
# Aleix Lafita - November 2019

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(ggplot2))

theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))


###################### Argparse #############################

seq1 = "examples/SasG_prot.fa"
seq2 = "examples/SasG_prot.fa"
window = 5
size = 1.0
prefix = "examples/SasG_prot"
type = "prot"

# create parser object
parser = ArgumentParser(description = 'Generate a dot-plot of two sequences')

# specify our desired options 
parser$add_argument("-x", "--seq1", default=seq1,
                    help="First sequence (horizontal) in FASTA format [default \"%(default)s\"]")
parser$add_argument("-y", "--seq2", default=seq2,
                    help="Second sequence (vertical) in FASTA format [default \"%(default)s\"]")
parser$add_argument("-w", "--window", default=window,
                    help="Window size of the comparison [default \"%(default)s\"]")
parser$add_argument("-s", "--size", default=size,
                    help="Scaling size for the plot dimensions [default \"%(default)s\"]")
parser$add_argument("-t", "--type", default=type,
                    help="Sequence type, either 'prot' or 'nucl' [default \"%(default)s\"]")
parser$add_argument("-p", "--prefix", default=prefix,
                    help="Prefix for the output dot-plot file [default \"%(default)s\"]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args = parser$parse_args()

seq1 = args$seq1
seq2 = args$seq2
window = as.integer(args$window)
size = as.integer(args$size)
prefix = args$prefix
type = args$type

if (!(type == "prot" | type == "nucl"))
    stop(sprintf("Option 'type' must be either 'prot' or 'nucl', but it is: '%s'", type))

###################### Parsing #############################
message("# Parsing sequences...")

# Parse the sequences and convert to DF
seq1.fasta = read.fasta(seq1)[1]
seq1.len = length(getSequence(seq1.fasta)[[1]])
seq1.name = names(seq1.fasta)

seq1.df = data.frame(
  id = seq1.name,
  pos = seq(1, seq1.len, 1),
  seq = unlist(getSequence(seq1.fasta)),
  stringsAsFactors = T
) %>% rowwise() %>%
  mutate(
    #seq = toupper(seq),
    seqw = paste0(unlist(getSequence(seq1.fasta))[pos:min((pos+window)-1, seq1.len)], collapse = ""),
    # If it is not DNA the complement gets assigned the same sequence
    seqwc = ifelse(type != "nucl", seqw, paste0(comp(unlist(getSequence(seq1.fasta))[pos:min((pos+window)-1, seq1.len)]), collapse = "")),
  ) %>% filter(nchar(seqw) == window)

seq2.fasta = read.fasta(seq2)[1]
seq2.len = length(getSequence(seq2.fasta)[[1]])
seq2.name = names(seq2.fasta)

seq2.df = data.frame(
  id = seq2.name,
  pos = seq(1, seq2.len, 1),
  seq = unlist(getSequence(seq2.fasta)),
  stringsAsFactors = T
) %>% rowwise() %>%
  mutate(
    #seq = toupper(seq),
    seqw = paste0(unlist(getSequence(seq2.fasta))[pos:min((pos+window)-1, seq2.len)], collapse = ""),
    seqwc = ifelse(type != "nucl", seqw, paste0(comp(unlist(getSequence(seq2.fasta))[pos:min((pos+window)-1, seq1.len)]), collapse = "")),
  ) %>% filter(nchar(seqw) == window)


###################### Calculate #############################
message(sprintf("# Generating dotplot matrix (%ix%i), this might take a while...", nrow(seq1.df), nrow(seq2.df)))

dotplot = merge(seq1.df, seq2.df, by =c())

message("# Extracting dotplot points, also takes a while...")

# First calculate the values, not memory efficient but cleaner
dotplot = dotplot %>%
  mutate(
    ident = seqw.x == seqw.y,
    rev = seqw.x == stri_reverse(seqw.y),
    comp = seqw.x == seqwc.y,
    revcomp = seqw.x == stri_reverse(seqwc.y)
  )

# Now filter the points that do not match any
dotplot = dotplot %>%
  filter(ident | rev | comp | revcomp)

# The seqinr package also has a dotplot function, but takes equally long or longer than this script
# dotplot(as.character(seq1.df$seq), as.character(seq2.df$seq))

###################### Plot #############################

# Plot the distribution of GC content
p = ggplot(dotplot, aes(x = pos.x, y = pos.y)) +
  #geom_tile(fill = "black") +
  geom_point(shape = "square", size = 1) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab(seq1.name) +
  ylab(seq2.name) +
  #scale_fill_gradient(low = "white", high = "black") +
  scale_x_continuous(limits = c(1, seq1.len), expand = c(0, 0), breaks = c(1, seq(100, seq1.len, 100))) +
  scale_y_reverse(limits = c(seq2.len, 1), expand = c(0, 0), breaks = c(1, seq(100, seq2.len, 100))) +
  coord_fixed()

pdf(sprintf("%s_dotplot.pdf", prefix), size * max(8, seq1.len/200), size * max(8, seq2.len/200))
plot(p)
log = dev.off()

message(sprintf("# Dotplot saved to %s_dotplot.pdf", prefix))
