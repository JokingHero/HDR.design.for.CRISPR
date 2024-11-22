rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

library(crisprScore)
library(Biostrings)

# sg1 ttccaagtggattctgctgg
# tcatgcagttcagatttgctcacccaactccccgtccatcagaaaaatgttccaagtggattctgctggaggaCtaCAgAaagcgggtgcagaacgtcactgagtttgatgacaggtgagtagtagttcagaaagcacatgtcccagg
# tcatgcagttcagatttgctcacccaactccCCGTCCATCAGAAAAATGTTCCAAGTGGATTCTGCTGGAGGATTATCAGAAGCGGGTGCAGAACGTCACTGAGTTTGATGACAGgtgagtagtagttcagaaagcacatgtcccagg
#                                                  ttccaagtggattctgctgg
# sg2 gctggaggattatcggaagc WT
# tcatgcagttcagatttgctcacccaactccccgtccatcagaaaaatgttccaagtggattcTGCTGGAGGATTATCGGAAGCGGgtgcagaacgtcactgagtttgatgacaggtgagtagtagttcagaaagcacatgtcccagg
#                                                                 gctggaggattatcggaagc
# sg3 tgctggaggattatcggaag WT
# tcatgcagttcagatttgctcacccaactccccgtccatcagaaaaatgttccaagtggattcTGCTGGAGGATTATCGGAAGCGGgtgcagaacgtcactgagtttgatgacaggtgagtagtagttcagaaagcacatgtcccagg
#                                                                tgctggaggattatcggaag
# sg4 ggattctgctggaggattat - removed because cant be tested in patient cells...
# tcatgcagttcagatttgctcacccaactccccgtccatcagaaaaatgttccaagtggattctgctggaggaCtaCAgAaagcgggtgcagaacgtcactgagtttgatgacaggtgagtagtagttcagaaagcacatgtcccagg
# tcatgcagttcagatttgctcacccaactccCCGTCCATCAGAAAAATGTTCCAAGTGGATTCTGCTGGAGGATTATCAGAAGCGGGTGCAGAACGTCACTGAGTTTGATGACAGgtgagtagtagttcagaaagcacatgtcccagg
#                                                          ggattctgctggaggattat
# sg5 atgttccaagtggattctgc
# tcatgcagttcagatttgctcacccaactccCCGTCCATCAGAAAAATGTTCCAAGTGGATTCTGCTGGAGGATTATCAGAAGCGGGTGCAGAACGTCACTGAGTTTGATGACAGgtgagtagtagttcagaaagcacatgtcccagg
#                                               atgttccaagtggattctgc
# sg6 catcagaaaaatgttccaag
# tcatgcagttcagatttgctcacccaactccCCGTCCATCAGAAAAATGTTCCAAGTGGATTCTGCTGGAGGATTATCAGAAGCGGGTGCAGAACGTCACTGAGTTTGATGACAGgtgagtagtagttcagaaagcacatgtcccagg
#                                     catcagaaaaatgttccaag
# sg7 atcctccagcagaatccact reverse AGTGGATTCTGCTGGAGGAT
# tcatgcagttcagatttgctcacccaactccCCGTCCATCAGAAAAATGTTCCAAGTGGATTCTGCTGGAGGATTATCAGAAGCGGGTGCAGAACGTCACTGAGTTTGATGACAGgtgagtagtagttcagaaagcacatgtcccagg
#                                                       AGTGGATTCTGCTGGAGGAT
spacer <- c("ttccaagtggattctgctgg", "gctggaggattatcggaagc", "tgctggaggattatcggaag",
            "atgttccaagtggattctgc", "catcagaaaaatgttccaag", "atcctccagcagaatccact") # "ggattctgctggaggattat",
pam <- c("AGG", "GGg", "CGG", "TGG", "TGG", "TGG") # "CGG", sg4 has CAG real PAM, but we input CGG because alghoritms
flank3_3bp <- c("ATT", "tgc", "gtg", "AGG", "AAC") # "AAG",
flank3_6bp <- c("ATTATC", "tgcaga", "gtgcag", "AGGATT", "AACATT") # "AAGCGG",
flank5_4bp <- c("AATG", "ttcT", "attc", "AAAA", "CGTC", "GATA") # "AAGT",
flank5_6bp <- c("AAAATG", "gattcT", "ggattc", "AGAAAA", "AATCCA") # "CCAAGT",

doench_2014 <- sapply(toupper(paste0(flank5_4bp, spacer, pam, flank3_3bp)),
                      function(x) crisprScore::getRuleSet1Scores(x)$score)
moreno_mateos_2015 <- sapply(toupper(paste0(flank5_6bp, spacer, pam, flank3_6bp)),
                             function(x) crisprScore::getCRISPRscanScores(x)$score)

kim_2019 <- sapply(toupper(paste0(flank5_4bp, spacer, pam, flank3_3bp)),
                           function(x) crisprScore::getDeepSpCas9Scores(x)$score)
labuhn_2018 <- sapply(toupper(spacer),
                      function(x) crisprScore::getCRISPRaterScores(x)$score)

# RANKS dont include sg4
dt <- data.frame(sg = paste0("sg", c(1:3, 5:7)),
                 PBMC_Kata = c(3, 2, 1, NA, NA, NA),
                 Fibroblasts_Frida = c(4, 2, 1, NA, NA, 3),
                 Stem_cells_Pavel = c(2, 1, 3, 4, 5, 6),
                 Stem_cells_Carolina = c(4, 3, 1, NA, NA, 2),
                 cell_free_screening = c(2, 1, 6, 4, 3, 5),
                 doench_2014 = rank(-1 * doench_2014),
                 moreno_mateos_2015 = rank(-1 * moreno_mateos_2015),
                 kim_2019 = rank(-1 * kim_2019),
                 labuhn_2018 = rank(-1 * labuhn_2018))
data.table::fwrite(dt, "ADA2_guides_scores.csv")
