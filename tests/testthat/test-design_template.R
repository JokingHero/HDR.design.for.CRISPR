# mutation loci on that gene is ENST00000349496.11: 110 in CDS coord, 25118 in transcript coordinates
# which is C > T, chr3:4122 4622, exon 3
# with exteionsion 4222-5022
# part of the genome used for CDS that overlaps our mutation loci + extension
# [1]     chr3 4526-4753      +
# [2]     chr3 4954-5022      +
# in CDS coordinates is
# [1]        14       241       228
# [2]       242       310        69


# genomic site: chr3 41224222-41225022
#

# annot <- system.file("data", "gencode.v42.annotation_mini.gff3", package = "HDR.design.for.CRISPR")
# txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(annot))
# withr::local_tempdir()
# test_that("get_cds", {
#   expect_equal(
#     sort(comb_along("NN", m = 1, letters = c("A", "B"))),
#     c("AN", "BN", "NA", "NB"))
# })

# design_one_template_for_each_guide(
#   ensemble_transcript_id = "ENST00000349496.11",
#   mutation_loci = 25118, # 25118
#   mutation_original = "C",
#   mutation_replacement = "T",
#   mutation_name = "C25118T",
#   output_dir = tmp,
#   annotation = "../gencode.v42.annotation.gff3",
#   genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#   snps = SNPlocs.Hsapiens.dbSNP155.GRCh38,
#   clinvar = "../clinvar.vcf.gz")


# ensemble_transcript_id = "ENST00000264657.10"
# mutation_loci = 1175
# mutation_original = "A"
# mutation_replacement = "G"
# mutation_name = "STAT3"
