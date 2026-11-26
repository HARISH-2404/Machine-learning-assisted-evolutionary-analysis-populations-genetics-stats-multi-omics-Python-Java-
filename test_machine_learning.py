import unittest
import numpy as np
from popgenml.machine_learning.regression import ridge_regression, lasso_coordinate_descent
from popgenml.machine_learning.classification import logistic_regression, classification_metrics
from popgenml.machine_learning.cross_validation import grouped_kfold_indices
from popgenml.machine_learning.feature_selection import select_top_features
from popgenml.machine_learning.model_evaluation import regression_metrics
from popgenml.machine_learning.hyperparameter_tuning import grid_search_ridge
from popgenml.machine_learning.explainability import permutation_importance
class TestML(unittest.TestCase):
 def setUp(self):
  self.X=np.array([[0,1],[1,0],[2,1],[3,2],[4,1],[5,3]],float);self.y=np.array([1,2,3,4,5,6],float)
 def test_regression(self):
  self.assertEqual(len(ridge_regression(self.X,self.y)["predictions"]),6)
  self.assertEqual(len(lasso_coordinate_descent(self.X,self.y)["coefficients"]),2)
 def test_classification(self):
  r=logistic_regression(self.X,np.array([0,0,0,1,1,1]),epochs=100)
  self.assertEqual(len(r["predictions"]),6)
  self.assertIn("f1",classification_metrics([0,1],[0,1]))
 def test_group_cv(self):
  folds=list(grouped_kfold_indices(["a","a","b","b","c","c"],3))
  self.assertEqual(len(folds),3)
  for tr,te in folds:self.assertFalse(set(np.array(["a","a","b","b","c","c"])[tr])&set(np.array(["a","a","b","b","c","c"])[te]))
 def test_features_metrics(self):
  self.assertEqual(len(select_top_features(self.X,self.y,1)["indices"]),1)
  self.assertGreater(regression_metrics(self.y,self.y)["r_squared"],.99)
 def test_grid_perm(self):
  g=grid_search_ridge(self.X,self.y,alphas=[.1,1],n_splits=3)
  self.assertIn("best_alpha",g)
  imp=permutation_importance(lambda x: ridge_regression(self.X,self.y)["intercept"]+((x-self.X.mean(0))/np.where(self.X.std(0)==0,1,self.X.std(0)))@ridge_regression(self.X,self.y)["coefficients"],self.X,self.y,n_repeats=2)
  self.assertEqual(len(imp),2)
if __name__=="__main__":unittest.main()
