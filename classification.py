"""Binary logistic regression and threshold-based classification metrics."""
import numpy as np

def _sigmoid(z): return 1/(1+np.exp(-np.clip(z,-35,35)))
def logistic_regression(X,y,learning_rate=0.05,epochs=3000,l2=0.01,threshold=0.5):
    X=np.asarray(X,float); y=np.asarray(y,float)
    if X.ndim!=2 or y.ndim!=1 or len(y)!=len(X) or set(np.unique(y))-{0.0,1.0}: raise ValueError("Binary y must contain 0 and 1")
    mu=X.mean(0); sd=X.std(0); sd[sd==0]=1; Z=(X-mu)/sd
    A=np.column_stack([np.ones(len(Z)),Z]); w=np.zeros(A.shape[1])
    for _ in range(epochs):
        p=_sigmoid(A@w); grad=A.T@(p-y)/len(y); grad[1:]+=l2*w[1:]
        w-=learning_rate*grad
    prob=_sigmoid(A@w); labels=(prob>=threshold).astype(int)
    return {"coefficients":w[1:],"intercept":w[0],"probabilities":prob,"predictions":labels,
            "feature_mean":mu,"feature_scale":sd}

def classification_metrics(y_true,y_pred):
    y=np.asarray(y_true,int); p=np.asarray(y_pred,int)
    if y.shape!=p.shape or y.ndim!=1: raise ValueError("Shapes must match")
    tp=int(((y==1)&(p==1)).sum()); tn=int(((y==0)&(p==0)).sum())
    fp=int(((y==0)&(p==1)).sum()); fn=int(((y==1)&(p==0)).sum())
    precision=tp/(tp+fp) if tp+fp else 0.0; recall=tp/(tp+fn) if tp+fn else 0.0
    return {"accuracy":(tp+tn)/len(y) if len(y) else float("nan"),"precision":precision,
            "recall":recall,"f1":2*precision*recall/(precision+recall) if precision+recall else 0.0,
            "confusion_matrix":[[tn,fp],[fn,tp]]}
