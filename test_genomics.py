import os, tempfile, unittest
from pathlib import Path
import numpy as np
from popgenml.genomics import sequence_statistics, extract_interval, load_vcf_dosage, variant_statistics, filter_variants
from popgenml.genomics import read_gff_genes, annotate_variant, window_coordinates, genes_to_bed, extract_gene_sequences

class TestGenomics(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.d=Path(self.tmp.name)
        (self.d/"g.fa").write_text(">chr1\nACGTACGTNN\n")
        (self.d/"g.gff3").write_text("##gff-version 3\nchr1\ttest\tgene\t2\t8\t.\t+\t.\tID=gene1;Name=G1\n")
        (self.d/"v.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\nchr1\t3\trs1\tG\tA\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\nchr1\t5\trs2\tA\tC\t40\tPASS\t.\tGT\t0/1\t./.\t0/0\n")
    def tearDown(self): self.tmp.cleanup()
    def test_fasta(self):
        s=sequence_statistics("ACGTNN")
        self.assertEqual(s["length"],6); self.assertEqual(s["gc_percent"],50.0)
        self.assertEqual(extract_interval(self.d/"g.fa","chr1",2,4),"CGT")
    def test_vcf(self):
        m,x=load_vcf_dosage(self.d/"v.vcf",snps_only=True)
        self.assertEqual(x.shape,(2,3)); self.assertEqual(x.iloc[0].tolist(),[0,1,2])
        st=variant_statistics(x.to_numpy()); self.assertAlmostEqual(st[0]["maf"],0.5)
        self.assertEqual(filter_variants(x.to_numpy(),0.5,0.3).tolist(),[True,False])
    def test_annotation_windows(self):
        genes=read_gff_genes(self.d/"g.gff3")
        self.assertEqual(annotate_variant("chr1",4,genes)["gene_id"],"gene1")
        self.assertEqual(window_coordinates(250,100,100),[(1,100),(101,200),(201,250)])

    def test_gwas_qtl_grm(self):
        from popgenml.genomics import single_marker_gwas, marker_qtl_scan, vanraden_relationship, benjamini_hochberg
        rng=np.random.default_rng(11)
        G=rng.binomial(2,0.3,size=(80,12)).astype(float)
        y=1.5*G[:,2]+rng.normal(0,1,size=80)
        gwas=single_marker_gwas(G,y,min_samples=20)
        self.assertEqual(len(gwas),12)
        self.assertLess(gwas.loc[2,"p_value"],0.01)
        qtl=marker_qtl_scan(G,y,[{"chromosome":"1","position":i*1000} for i in range(12)],min_samples=20)
        self.assertEqual(len(qtl),12)
        K=vanraden_relationship(G)
        self.assertEqual(K.shape,(80,80))
        self.assertTrue(np.allclose(K,K.T))
        q=benjamini_hochberg([0.01,0.04,0.2])
        self.assertTrue(np.all(q>=np.array([0.01,0.04,0.2])))

    def test_gene_outputs(self):
        self.assertEqual(genes_to_bed(self.d/"g.gff3",self.d/"g.bed"),1)
        self.assertEqual(extract_gene_sequences(self.d/"g.gff3",self.d/"g.fa",self.d/"genes.fa"),1)

if __name__=="__main__": unittest.main()
