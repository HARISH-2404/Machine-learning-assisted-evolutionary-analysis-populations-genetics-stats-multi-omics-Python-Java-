import unittest
from popgenml.machine_learning.deep_learning_models import DEEP_MODELS
class TestAdvancedCatalogue(unittest.TestCase):
 def test_catalogue(self): self.assertEqual(len(DEEP_MODELS),10)
if __name__=="__main__":unittest.main()
