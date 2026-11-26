import unittest
from popgenml.molecular_genetics.sequence import reverse_complement, gc_content
from popgenml.molecular_genetics.translation import translate, transcribe
from popgenml.molecular_genetics.codon_usage import codon_usage
from popgenml.molecular_genetics.motif_scanning import scan_motif
from popgenml.molecular_genetics.primer_design import primer_metrics
from popgenml.molecular_genetics.gene_structure import summarize_gene
from popgenml.molecular_genetics.protein_properties import protein_properties

class TestMolecularGenetics(unittest.TestCase):
 def test_sequence(self):
  self.assertEqual(reverse_complement("ATGC"),"GCAT")
  self.assertAlmostEqual(gc_content("GGCCAT"),66.6666667,places=5)
 def test_translation(self):
  self.assertEqual(transcribe("ATG"),"AUG")
  self.assertEqual(translate("ATGGCTTAA"),"MA*")
 def test_codon(self):
  self.assertEqual(codon_usage("ATGGCT")["complete_codons"],2)
 def test_motif(self):
  self.assertEqual([x["start"] for x in scan_motif("ATATAT","ATA")],[1,3])
 def test_primer(self):
  self.assertEqual(primer_metrics("ATGCATGC")["length"],8)
 def test_gene(self):
  r=summarize_gene([(1,10),(21,30)])
  self.assertEqual(r["intron_intervals"],[(11,20)])
 def test_protein(self):
  self.assertEqual(protein_properties("MA")["length_aa"],2)
if __name__=="__main__": unittest.main()
