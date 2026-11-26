import unittest
from popgenml.quantitative_genetics.variance_components import one_way_random_effects
from popgenml.quantitative_genetics.heritability import broad_sense_heritability, entry_mean_heritability
from popgenml.quantitative_genetics.genetic_advance import genetic_advance
from popgenml.quantitative_genetics.combining_ability import gca_sca
from popgenml.quantitative_genetics.gxe import stability_parameters
class TestQuantitative(unittest.TestCase):
 def test_variance(self):
  r=one_way_random_effects({"A":[3,4,5],"B":[7,8,9]})
  self.assertGreater(r["genetic_variance"],0)
 def test_h2(self):
  self.assertAlmostEqual(broad_sense_heritability(2,4),.5)
  self.assertGreater(entry_mean_heritability(2,1,3,3,2),0)
 def test_ga(self):
  self.assertAlmostEqual(genetic_advance(10,2,.5,1.4)["genetic_advance"],1.4)
 def test_gca(self):
  self.assertIn("gca",gca_sca({("A","B"):10,("A","C"):12,("B","C"):8}))
 def test_gxe(self):
  rec=[]
  vals={"A":[2,4,6],"B":[3,5,7]}
  for g,vs in vals.items():
   for e,v in enumerate(vs):rec.append({"genotype":g,"environment":str(e),"trait":v})
  self.assertAlmostEqual(stability_parameters(rec)["A"]["environmental_slope"],1)
if __name__=="__main__":unittest.main()
