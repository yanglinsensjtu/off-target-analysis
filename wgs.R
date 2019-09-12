library(VariantAnnotation)
library(vcfR)
wgs.less100 <- readVcfAsVRanges('../somatic-WGS/R19026469LD01-T0502-CAR_vs_R19026468LD01-T0501-MOCK_TN.PASS.vcf-anno.for_crispr.vcf')
wgs.more100 <- read.vcfR('../somatic-WGS/CAR_vs_Normal_somaticSV.vcf.PASS.vcf-anno.vcf')
wgs.less100.gr <- GRanges(seqnames = wgs.less100@seqnames[1:37],
                          ranges = IRanges(start = start(wgs.less100@ranges)[1:37]-200,
                                           width = 400),
                          strand = wgs.less100@strand[1:37])

wgs.more100.fix <- as.data.frame(wgs.more100@fix)

wgs.more100.gr <- GRanges(seqnames = wgs.more100.fix$CHROM,
                          ranges = IRanges(start = as.integer(as.character(wgs.more100.fix$POS)) -200,
                                           width = 400),
                          strand = rep('*', 18))

wgs.less100.jc <- readVcfAsVRanges('../jointcall-WGS/JC_SNP_INDEL_recal.PASS.vcf-anno.for_crispr.vcf')
wgs.more100.jc <- read.vcfR('../jointcall-WGS/R19026469LD01-T0502-CAR_diploidSV.vcf.PASS.vcf-anno.rm_MOCK.vcf')

wgs.less100.jc.gr <- GRanges(seqnames = wgs.less100.jc@seqnames,
                             ranges = IRanges(start = start(wgs.less100.jc@ranges)-200,
                                              width = 400),
                             strand = wgs.less100.jc@strand)
wgs.more100.jc.fix <- as.data.frame(wgs.more100.jc@fix)
wgs.more100.jc.gr <- GRanges(seqnames = wgs.more100.jc.fix$CHROM,
                             ranges = IRanges(start = as.integer(as.character(wgs.more100.jc.fix$POS)) - 200,
                                              width = 400),
                             strand = rep('*',nrow(wgs.more100.jc.fix)))
wgs.jc <- c(wgs.more100.jc.gr, wgs.less100.jc.gr)
