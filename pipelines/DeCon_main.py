"""
Modules for machine learning using scikit-learn.

Dependencies:
python module: numpy, pandas, sklearn
"""
################################################################################
"""
Classifiers
Functions with name *_classifier takes arguments:
  training_sample: Array-like of shape (n_samples, n_features).
                   The training input samples.
  training_results: Array-like of shape (1, n_samples).
                    The target labels, real numbers or strings.
  testing: Array-like of shape (m_samples, n_features).
           The input samples.
and returns:
  Array-like of shape (1, m_samples).
  The predicted labels for samples defined by rows in argument "testing".
"""
def decision_tree_classifier(training_sample,training_results,testing):
  # """Decision tree"""
  import numpy as np
  import pandas as pd
  from sklearn.tree import DecisionTreeClassifier
  weight=pd.DataFrame({"class":training_results}).groupby("class",as_index=False).value_counts()
  weight=dict(weight.values)
  classifier=DecisionTreeClassifier(class_weight=weight).fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
def random_forest_classifier(training_sample,training_results,testing):
  # """Random forest"""
  import numpy as np
  import pandas as pd
  from sklearn.ensemble import RandomForestClassifier
  weight=pd.DataFrame({"class":training_results}).groupby("class",as_index=False).value_counts()
  weight=dict(weight.values)
  classifier=RandomForestClassifier(class_weight=weight).fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
def bagging_classifier(training_sample,training_results,testing):
  # """Bagging classifier"""
  import numpy as np
  import pandas as pd
  from sklearn.ensemble import BaggingClassifier
  classifier=BaggingClassifier().fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
def linear_SV_classifier(training_sample,training_results,testing):
  # """Linear support vector"""
  import numpy as np
  import pandas as pd
  from sklearn.svm import LinearSVC
  weight=pd.DataFrame({"class":training_results}).groupby("class",as_index=False).value_counts()
  weight=dict(weight.values)
  classifier=LinearSVC(max_iter=10000,class_weight=weight).fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
def nu_SV_classifier(training_sample,training_results,testing):
  """nu support vector"""
  import numpy as np
  import pandas as pd
  from sklearn.svm import NuSVC
  weight=pd.DataFrame({"class":training_results}).groupby("class",as_index=False).value_counts()
  weight=dict(weight.values)
  classifier=NuSVC(class_weight=weight).fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
def k_neighbors_classifier(training_sample,training_results,testing):
  """k neighbors"""
  import numpy as np
  import pandas as pd
  from sklearn.neighbors import KNeighborsClassifier
  classifier=KNeighborsClassifier().fit(training_sample,training_results)
  value=classifier.predict(testing)
  return value
################################################################################
"""
Outlier detectors
The training data contains outliers, i.e. observations that are far from the others.
Outlier detectors try to fit the regions where the training data is the most
concentrated, ignoring the deviant observations.

Functions with name *_outlier_detector takes arguments:
  training: Array-like of shape (n_samples, n_features).
            The training input samples.
  contamination: float in (0, 0.5]
                 The proportion of outliers in training.
                 It defines the threshold on the scores of the samples.
  testing: Array-like of shape (m_samples, n_features).
           The input samples.
and returns:
  Array-like of shape (1, m_samples).
  Tell whether samples defined by rows in argument "testing" are inliers/outliers
  of "training". -1 for outliers and 1 for inliers.
"""
def isolation_forest_outlier_detector(training,contamination,testing):
  """Isolation forest"""
  import numpy as np
  import pandas as pd
  from sklearn.ensemble import IsolationForest
  outlier_detector=IsolationForest(contamination=contamination).fit(training)
  value=outlier_detector.predict(testing)
  return value
def lof_outlier_detector(training,contamination,testing):
  """Local outlier factor"""
  import numpy as np
  import pandas as pd
  from sklearn.neighbors import LocalOutlierFactor
  outlier_detector=LocalOutlierFactor(contamination=contamination).fit(training)
  value=outlier_detector.predict(testing)
  return value
################################################################################
"""
Clustering
Functions with name *_clustering takes arguments:
  samples: Array-like of shape (n_samples, n_features).
           The samples to be clustered.
and returns
  Array-like of shape (1, m_samples).
  Index of the cluster each sample belongs to.
"""
def kmean_clustering(samples,k):
  """kmean clustering"""
  import numpy as np
  import pandas as pd
  from sklearn.cluster import KMeans
  clustering=KMeans(n_clusters=k).fit(samples)
  values=clustering.predict(samples)
  return values
################################################################################
"""
Clustering evaluation
Functions with name *_clustering_evaluator takes argument:
  samples: Array-like of shape (n_samples, n_features).
           Samples clustered.
  clusters: array-like of shape (1, n_samples).
            Predicted cluster indexes for samples.
"""
def SC(samples,clusters):
  """mean Silhouette Coefficient
     The best value is 1 and the worst value is -1.
     Values near 0 indicate overlapping clusters.
     Negative values generally indicate that a sample
     has been assigned to the wrong cluster, as a
     different cluster is more similar."""
  import numpy as np
  import pandas as pd
  from sklearn import metrics
  sc=metrics.silhouette_score(samples,clusters)
  return sc
def CH(samples,clusters):
  """Calinski and Harabasz score
     It is ratio of the sum of between-cluster dispersion
     and of within-cluster dispersion. Bigger value means
     better clustering."""
  import numpy as np
  import pandas as pd
  from sklearn import metrics
  ch=metrics.calinski_harabasz_score(samples,clusters)
  return ch
def DB(samples,clusters):
  """Davies-Bouldin score.
     The minimum score is zero, with lower values indicating better clustering."""
  from sklearn import metrics
  db=metrics.davies_bouldin_score(samples,clusters)
  return(db)
