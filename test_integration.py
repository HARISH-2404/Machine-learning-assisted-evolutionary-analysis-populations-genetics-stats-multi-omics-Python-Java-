import pandas as pd
from popgenml.integration import merge_omics_blocks, scale_feature_block, late_integration_scores

def test_merge_prefix_and_intersection():
    a=pd.DataFrame({"x":[1,2]},index=["s1","s2"])
    b=pd.DataFrame({"x":[3,4]},index=["s2","s3"])
    out=merge_omics_blocks({"a":a,"b":b})
    assert list(out.index)==["s2"] and list(out.columns)==["a::x","b::x"]

def test_scale_uses_training_stats():
    tr=pd.DataFrame({"x":[0.,2.]}); te=pd.DataFrame({"x":[4.]})
    a,b=scale_feature_block(tr,te)
    assert abs(a.x.mean())<1e-9 and b.x.iloc[0]>1

def test_late_score():
    a=pd.Series([0.,1.],index=["s1","s2"]); b=pd.Series([1.,1.],index=["s1","s2"])
    out=late_integration_scores({"a":a,"b":b})
    assert list(out)==[0.5,1.0]
