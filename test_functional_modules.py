import unittest
from popgenml.molecular_genetics.go_enrichment import go_enrichment
from popgenml.molecular_genetics.kegg_enrichment import kegg_enrichment
from popgenml.molecular_genetics.regulatory_network import build_regulatory_network, network_summary, reachable_targets
class TestFunctionalModules(unittest.TestCase):
 def test_go(self):
  m={"g1":["GO:1"],"g2":["GO:1"],"g3":["GO:2"],"g4":["GO:2"]}
  r=go_enrichment(["g1","g2"],m,background_genes=m,min_term_size=1)
  self.assertEqual(r[0]["overlap"],2)
  self.assertIn("p_adjusted",r[0])
 def test_kegg(self):
  m={"g1":["p1"],"g2":["p1"],"g3":["p2"]}
  self.assertTrue(kegg_enrichment(["g1"],m,background_genes=m,min_pathway_size=1))
 def test_network(self):
  n=build_regulatory_network([{"regulator":"TF1","target":"g1"},{"regulator":"g1","target":"g2"}])
  self.assertEqual(network_summary(n)["edge_count"],2)
  self.assertEqual(reachable_targets(n,"TF1"),["g1","g2"])
if __name__=="__main__":unittest.main()
