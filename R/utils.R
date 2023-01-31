melting_temp <- function(x) TmCalculator::Tm_GC(
  ntseq = as.character(x),
  variant = "Primer3Plus", Na = 50, outlist = FALSE)

gc_fract <- function(x) letterFrequency(x, letters = "CG", as.prob = TRUE)

design_probes <- function(mutation_name, st, sp, s, tmin = 59, tmax = 61, len_min = 20, len_max = 25,
                          origin_mut_start = 400 + 1) {
  probes <- GRanges()
  for (len in len_min:len_max) { # for each length
    for (i in st:(sp-len)) { # for each bases
      si <- s[i:(i+len)]
      tm <- melting_temp(si)
      if (tm < tmin | tm > tmax) next
      probes <- c(probes,
                  GRanges(seqnames = mutation_name,
                          ranges = IRanges(start = i, width = len),
                          strand = "+",
                          original = as.character(si),
                          shift = i - origin_mut_start,
                          length = len,
                          Tm = tm,
                          GC = gc_fract(si)))
    }
  }
  return(probes)
}
