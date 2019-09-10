library(VariantAnnotation)
library(vcfR)
wgs.less100 <- readVcfAsVRanges('../WGS/R19026469LD01-T0502-CAR_vs_R19026468LD01-T0501-MOCK_TN.PASS.vcf-anno.for_crispr.vcf')
wgs.more100 <- read.vcfR('../WGS/CAR_vs_Normal_somaticSV.vcf.PASS.vcf-anno.vcf')
wgs.less100.gr <- GRanges(seqnames = wgs.less100@seqnames[1:37],
                          ranges = IRanges(start = start(wgs.less100@ranges)[1:37]-200,
                                           width = 400),
                          strand = wgs.less100@strand[1:37])

wgs.more100.fix <- as.data.frame(wgs.more100@fix)

wgs.more100.gr <- GRanges(seqnames = wgs.more100.fix$CHROM,
                          ranges = IRanges(start = as.integer(as.character(wgs.more100.fix$POS)) -200,
                                           width = 400),
                          strand = rep('*', 18))

