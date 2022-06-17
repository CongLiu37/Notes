Modules for machine learning using scikit-learn.

Dependencies:
python module: numpy, pandas, sklearn

Classifiers:
(1) decision_tree_classifier
(2) random_forest_classifier
(3) bagging_classifier
(4) linear_SV_classifier
(5) nu_SV_classifier
(6) k_neighbors_classifier
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

Outlier detectors
(1) isolation_forest_outlier_detector
(2) lof_outlier_detector
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

Clustering
(1) kmean_clustering
Functions with name *_clustering takes arguments:
  samples: Array-like of shape (n_samples, n_features).
           The samples to be clustered.
and returns
   Array-like of shape (1, m_samples).
  Index of the cluster each sample belongs to.

Clustering evaluation
(1) SC
(2) CH
(3) DB
Functions with name *_clustering_evaluator takes argument:
  samples: Array-like of shape (n_samples, n_features).
           Samples clustered.
  clusters: array-like of shape (1, n_samples).
            Predicted cluster indexes for samples.
