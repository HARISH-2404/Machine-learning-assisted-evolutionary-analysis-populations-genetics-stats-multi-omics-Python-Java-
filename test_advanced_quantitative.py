import unittest
import numpy as np
from popgenml.quantitative_genetics.gblup import gblup
from popgenml.quantitative_genetics.rrblup import rrblup
from popgenml.quantitative_genetics.bayesian_models import bayesian_ridge_gibbs
from popgenml.quantitative_genetics.multi_trait import multi_trait_rrblup
class TestAdvanced(unittest.TestCase):
 def setUp(self):
  self.X=np.array([[0,1],[1,0],[2,1],[1,2],[0,2]],float)
  self.y=np.array([1,2,4,3,2],float)
 def test_rrblup(self):
  self.assertEqual(len(rrblup(self.X,self.y)["predictions"]),5)
 def test_gblup(self):
  G=self.X@self.X.T/2+np.eye(5)*.01
  self.assertEqual(len(gblup(G,self.y)["predictions"]),5)
 def test_bayes(self):
  r=bayesian_ridge_gibbs(self.X,self.y,iterations=80,burn_in=30,thin=5,seed=7)
  self.assertGreater(r["posterior_samples"],0)
 def test_multi(self):
  Y=np.column_stack([self.y,self.y*2])
  self.assertEqual(multi_trait_rrblup(self.X,Y)["predictions"].shape,(5,2))
if __name__=="__main__":unittest.main()
