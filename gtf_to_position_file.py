from RiboCode import prepare_transcripts
gtfFile = "refdata-cellranger-mm10-3.0.0/genes/genes.gtf" #download from 10X genomic website
gene_dict,_ = prepare_transcripts.readGTF(gtfFile)
chrs = map(str,range(1,20))
with open("mm10_gene_pos.txt","w") as fout:
  gnames = set()
  for c in chrs:
    genes = [g for g in gene_dict.values() if g.chrom == c]
    sorted_genes = sorted(genes,key=lambda x: x.genomic_iv.start)
    for g in sorted_genes:
      if g.gene_name not in gnames:
        tmp = [g.gene_name,c,g.genomic_iv.start+1,g.genomic_iv.end]
        fout.write("\t".join(map(str,tmp)) + "\n")
        gnames.add(g.gene_name)
