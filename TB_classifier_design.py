#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Martínez-Enguita (2023)

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.feature_selection import SelectFromModel, RFE
from sklearn.linear_model import LogisticRegression, ElasticNet
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve, auc
from sklearn.utils import resample

# %%

# Set main and input directories
main_dir = "./"
input_dir = "Data/"
output_dir = "Results/"

# %%

# Load DMC beta values and metadata
dmc_input = pd.read_csv((main_dir + input_dir + "DMC_Linkoping.csv"),
                        index_col=0)

dmc_meta = pd.read_csv((main_dir + input_dir + "TB_metadata.csv"),
                       index_col=0)

# Randomize DMC order,  and binarize sample group
dmc_input = dmc_input.sample(frac = 1)

dmc_meta["Group_bin"] = (dmc_meta["Group"] == "Patient")*1

# Divide into training/validation set and test set
print(all(dmc_input.columns == dmc_meta.loc[:, "Sample_Name"]))

trainval_meta = dmc_meta[dmc_meta["Country"].str.contains("Sweden", case=False)]
trainval_set = dmc_input.loc[:, trainval_meta.loc[:, "Sample_Name"]]

test_meta = dmc_meta[(dmc_meta["Country"].str.contains("Kenya", case=False)) | (dmc_meta["Country"].str.contains("Peru", case=False))]
test_set = dmc_input.loc[:, test_meta.loc[:, "Sample_Name"]]

# Split trainval data
X_train, X_test, y_train, y_test = train_test_split(np.array(np.transpose(trainval_set)),
                                                    np.transpose(np.array(trainval_meta.loc[:, "Group_bin"])),
                                                    test_size=0.2,
                                                    random_state=7777)

# %%

## Univariate Logistic Regression for Feature Selection
univariate_model = SelectFromModel(LogisticRegression(penalty="l1", solver="liblinear", C=10))
X_train_selected_ = univariate_model.fit_transform(X_train, y_train)

# Selected DMCs for top features
selected_dmc_indices = univariate_model.get_support(indices=True)
selected_dmc_names = np.array(dmc_input.index[selected_dmc_indices])


## Multivariate Elastic Net Regression for Feature Selection
elastic_net_model = ElasticNet(alpha=1.0, l1_ratio=0.5)
elastic_net_model.fit(X_train, y_train)

# Get feature importances from the Elastic Net model
feature_importances = np.abs(elastic_net_model.coef_)

# Select top features based on importance scores
num_top_features = 10
selected_dmc_indices = np.argsort(feature_importances)[-num_top_features:]
selected_dmc_names = np.array(dmc_input.index[selected_dmc_indices])

# %%

## C parameter grid search
# Grid of C values for logistic regression
c_values = np.arange(1, 100)

# Perform grid search (Random Forest CV) to find best C value
gridsearch_results = []

for c in c_values:
    logistic_model = LogisticRegression(penalty="l1",
                                        solver="liblinear",
                                        C=c)
    univariate_model = SelectFromModel(logistic_model)
    X_train_selected = univariate_model.fit_transform(X_train, y_train)

    # Train and evaluate Random Forest on selected features
    clf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=777)
    clf.fit(X_train_selected, y_train)
    y_pred_prob = clf.predict_proba(np.array(np.transpose(test_set))[:, univariate_model.get_support()])[:, 1]
    auc_out = roc_auc_score(np.transpose(np.array(test_meta.loc[:, "Group_bin"])), y_pred_prob)

    # Store iteration results
    it_out_indices = univariate_model.get_support(indices=True)
    it_out_names = np.array(dmc_input.index[selected_dmc_indices])

    it_out = {"C": c,
              "Selected_Indices": ",".join(map(str, it_out_indices)),
              "Selected_Names": ",".join(map(str, it_out_names)),
              "Number": len(it_out_indices),
              "AUC": auc_out}
    gridsearch_results.append(it_out)

gridsearch_results_df = pd.DataFrame(gridsearch_results)

# Find the row with the highest AUC and its corresponding information
best_result = gridsearch_results_df.loc[gridsearch_results_df["AUC"].idxmax()]

print("Best C value:", best_result["C"])
print("Best AUC:", best_result["AUC"])
print("Selected Names:", best_result["Selected_Names"])
print("Selected Indices:", best_result["Selected_Indices"])

# Store grid search output
gridsearch_results_df.to_csv((main_dir + output_dir + "gridsearch_logreg_univariate.csv"), index=False)

# %%

## C parameter and number of features grid search
# Grid of number of features and C values for logistic regression
c_values = np.arange(1, 100)
gridsearch_results = []
n_sel_feat = 20

# Perform grid search with cross-validation to find best C
for c in c_values:
    logistic_model = LogisticRegression(penalty='l1', solver='liblinear', C=c)
    rfe_model = RFE(logistic_model, n_features_to_select=n_sel_feat)
    X_train_selected = rfe_model.fit(X_train, y_train)

    # Transform the training data to retain only the selected features
    X_train_selected = rfe_model.transform(X_train)

    # Obtain indices and names of selected features
    selected_features_mask = rfe_model.support_
    selected_feature_indices = [index for index, mask_value in enumerate(selected_features_mask) if mask_value]
    selected_feature_names = [dmc_input.index[index] for index in selected_feature_indices]
    print("Selected feature names:", selected_feature_names)

    # Train and evaluate Random Forest on selected features
    clf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=777)
    clf.fit(X_train_selected, y_train)
    y_pred_prob = clf.predict_proba(np.array(np.transpose(test_set))[:, selected_features_mask])[:, 1]
    auc_out = roc_auc_score(np.transpose(np.array(test_meta.loc[:, "Group_bin"])), y_pred_prob)

    # Store iteration results
    it_out_indices = selected_feature_indices
    it_out_names = selected_feature_names

    it_out = {"C": c,
              "Selected_Indices": ",".join(map(str, it_out_indices)),
              "Selected_Names": ",".join(map(str, it_out_names)),
              "Number": len(it_out_indices),
              "AUC": auc_out}
    gridsearch_results.append(it_out)

# Convert results to a DataFrame for easier analysis
gridsearch_results_df_RFE = pd.DataFrame(gridsearch_results)

# Find the row with the highest AUC and its corresponding information
best_result = gridsearch_results_df_RFE.loc[gridsearch_results_df_RFE["AUC"].idxmax()]

print("Best C value:", best_result["C"])
print("Best AUC:", best_result["AUC"])
print("Selected Names:", best_result["Selected_Names"])
print("Selected Indices:", best_result["Selected_Indices"])

gridsearch_results_df_RFE.to_csv((main_dir + output_dir + "gridsearch_results_RFE" + str(n_sel_feat) + "_validation.csv"), index=False)

# %%
## Evaluate performance of selected features
selected_dmc_indices = np.array([int(i) for i in best_result["Selected_Indices"].split(",")])
selected_dmc_names = np.array(dmc_input.index[selected_dmc_indices])

# Define classifiers
classifiers = [
    ("SVM", SVC(probability=True, C=1.0, kernel="rbf")),
    ("KNN", KNeighborsClassifier(n_neighbors=5)),
    ("Random Forest", RandomForestClassifier(n_estimators=100, max_depth=10, random_state=777)),
    ("XGBoost", XGBClassifier(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=777))
]

# Train and evaluate classifiers
for clf_name, clf in classifiers:
    clf.fit(X_train[:, selected_dmc_indices], y_train)
    y_pred_prob = clf.predict_proba(X_test[:, selected_dmc_indices])[:, 1]
    y_pred = (y_pred_prob > 0.5).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    specificity = tn / (tn + fp)
    sensitivity = tp / (tp + fn)
    auc_out = roc_auc_score(y_test, y_pred_prob)

    print(f"Classifier: {clf_name}")
    print(f"Specificity: {specificity:.2f}")
    print(f"Sensitivity: {sensitivity:.2f}")
    print(f"AUC: {auc_out:.2f}")
    print("--------------------------")

# Evaluate validation set
test_sets = [
    ("Training Set", np.array(np.transpose(trainval_set)), np.transpose(np.array(trainval_meta.loc[:, "Group_bin"]))),
    ("Test Set", np.array(np.transpose(test_set)), np.transpose(np.array(test_meta.loc[:, "Group_bin"])))
]

# Define the number of bootstrap samples
num_bootstrap_samples = 1000

# Calculate the 95% CI for AUC using bootstrapping
results = []

for test_set_name, test_data, test_labels in test_sets:
    print(f"Evaluating on {test_set_name}:")
    for clf_name, clf in classifiers:
        auc_values = []
        sensitivity_values = []
        specificity_values = []
        for _ in range(num_bootstrap_samples):
            # Bootstrap resampling of the test set
            boot_indices = np.random.choice(range(len(test_data)), size=len(test_data), replace=True)
            boot_data = test_data[boot_indices]
            boot_labels = test_labels[boot_indices]

            try:
                y_pred_prob = clf.predict_proba(boot_data[:, selected_dmc_indices])[:, 1]
                auc_out = roc_auc_score(boot_labels, y_pred_prob)
                auc_values.append(auc_out)

                y_pred = (y_pred_prob > 0.5).astype(int)
                tn, fp, fn, tp = confusion_matrix(boot_labels, y_pred, labels=[0, 1]).ravel()
                sensitivity = tp / (tp + fn)
                specificity = tn / (tn + fp)
                sensitivity_values.append(sensitivity)
                specificity_values.append(specificity)
            except ValueError:
                # Skip the iteration if only one class is present in y_true
                pass

        if auc_values:
            lower_bound_auc = np.percentile(auc_values, 2.5)
            upper_bound_auc = np.percentile(auc_values, 97.5)
            avg_auc = np.mean(auc_values)

            lower_bound_sens = np.percentile(sensitivity_values, 2.5)
            upper_bound_sens = np.percentile(sensitivity_values, 97.5)
            avg_sens = np.mean(sensitivity_values)

            lower_bound_spec = np.percentile(specificity_values, 2.5)
            upper_bound_spec = np.percentile(specificity_values, 97.5)
            avg_spec = np.mean(specificity_values)

            print(f"Classifier: {clf_name}")
            print(f"AUC 95% CI: ({lower_bound_auc:.2f}, {upper_bound_auc:.2f}), Avg AUC: {avg_auc:.2f}")
            print(f"Sensitivity 95% CI: ({lower_bound_sens:.2f}, {upper_bound_sens:.2f}), Avg Sensitivity: {avg_sens:.2f}")
            print(f"Specificity 95% CI: ({lower_bound_spec:.2f}, {upper_bound_spec:.2f}), Avg Specificity: {avg_spec:.2f}")

            results.append({
                "Test Set": test_set_name,
                "Classifier": clf_name,
                "AUC 95% CI (Lower)": lower_bound_auc,
                "AUC 95% CI (Upper)": upper_bound_auc,
                "Avg AUC": avg_auc,
                "Sensitivity 95% CI (Lower)": lower_bound_sens,
                "Sensitivity 95% CI (Upper)": upper_bound_sens,
                "Avg Sensitivity": avg_sens,
                "Specificity 95% CI (Lower)": lower_bound_spec,
                "Specificity 95% CI (Upper)": upper_bound_spec,
                "Avg Specificity": avg_spec,
                "Selected DMCs": ", ".join(selected_dmc_names)
                })
        else:
            print(f"Classifier: {clf_name}")
            print("Evaluation skipped due to single class in y_true")

        print("--------------------------")

results_df_RFE = pd.DataFrame(results)

# Save the test set evaluation results
results_df_RFE.to_csv((main_dir + output_dir + "eval_test_set_RFE" + str(n_sel_feat) + "_validation.csv"), index=False)

# %%

## Line plot of results
# Sample data: AUC values for each number of selected CpGs (1 to 20)
num_selected_cpgs = list(range(1, 21))

auc_values = []
auc_values_val = []
sens_values = []
sens_values_val = []
spec_values = []
spec_values_val = []

out_pattern = "eval_test_set_RFE{0}_validation.csv"

# Iterate through numbers from 1 to 20 (classifier average)
for i in range(1, 21):
    # Construct the full path to the file
    filename = out_pattern.format(i)
    file_path = os.path.join((main_dir + output_dir), filename)

    # Check if the file exists before attempting to open it
    if os.path.isfile(file_path):
        rfe_out = pd.read_csv(file_path, index_col=False, header=0)
        auc_values.append(rfe_out.iloc[0:4, 4].mean())
        auc_values_val.append(rfe_out.iloc[4:8, 4].mean())
        sens_values.append(rfe_out.iloc[0:4, 7].mean())
        sens_values_val.append(rfe_out.iloc[4:8, 7].mean())
        spec_values.append(rfe_out.iloc[0:4, 10].mean())
        spec_values_val.append(rfe_out.iloc[4:8, 10].mean())

# Iterate through numbers from 1 to 20 (best classifier)
for i in range(1, 21):
    # Construct the full path to the file
    filename = out_pattern.format(i)
    file_path = os.path.join((main_dir + output_dir), filename)

    # Check if the file exists before attempting to open it
    if os.path.isfile(file_path):
        rfe_out = pd.read_csv(file_path, index_col=False, header=0)
        auc_values.append(rfe_out.iloc[0:4, 4].max())
        auc_values_val.append(rfe_out.iloc[4:8, 4].max())
        sens_values.append(rfe_out.iloc[0:4, 7].max())
        sens_values_val.append(rfe_out.iloc[4:8, 7].max())
        spec_values.append(rfe_out.iloc[0:4, 10].max())
        spec_values_val.append(rfe_out.iloc[4:8, 10].max())


## AUC lineplot
plt.figure(figsize=(8, 6))
plt.plot(num_selected_cpgs, auc_values, color='#236AB9', marker='o', label='Training Cohort')
plt.plot(num_selected_cpgs, auc_values_val, color='#FC7307', marker='o', label='Validation Cohort')

# Set Y-axis breaks from 0 to 1 in increments of 0.1
plt.yticks(np.arange(0.80, 1.01, 0.05))

# Set X-axis breaks from 5 to 20 in increments of 1
plt.xticks(np.arange(1, 21, 1))

# Set plot labels and title
plt.xlabel('Number of CpG Sites Selected')
plt.ylabel('AUC')
plt.title('Feature Selection Analysis - TB Classifier [Multivariate Regression - Best Classifier]\n\nAUC')
plt.legend()

# Set white background
plt.gca().set_facecolor('white')

# Show the plot
plt.grid(True)
auc_lineplot = plt.gcf()
plt.show()
auc_lineplot.savefig("./Results/AUC_evaluation_validation_lineplot.png", dpi=600)


## Sensitivity lineplot
plt.figure(figsize=(8, 6))
plt.plot(num_selected_cpgs, sens_values, color='#236AB9', marker='o', label='Training Cohort')
plt.plot(num_selected_cpgs, sens_values_val, color='#FC7307', marker='o', label='Validation Cohort')

# Set Y-axis breaks from 0 to 1 in increments of 0.1
plt.yticks(np.arange(0.6, 1.01, 0.05))

# Set X-axis breaks from 5 to 20 in increments of 1
plt.xticks(np.arange(1, 21, 1))

# Set plot labels and title
plt.xlabel('Number of CpG Sites Selected')
plt.ylabel('Sensitivity')
plt.title('Feature Selection Analysis - TB Classifier [Multivariate Regression - Best Classifier]\n\nSensitivity')
plt.legend()

# Set white background
plt.gca().set_facecolor('white')

# Show the plot
plt.grid(True)
sens_lineplot = plt.gcf()
plt.show()
sens_lineplot.savefig("./Results/Sens_evaluation_validation_lineplot.png", dpi=600)


## Specificity lineplot
plt.figure(figsize=(9, 6))
plt.plot(num_selected_cpgs, spec_values, color='#236AB9', marker='o', label='Training Cohort')
plt.plot(num_selected_cpgs, spec_values_val, color='#FC7307', marker='o', label='Validation Cohort')

# Set Y-axis breaks from 0 to 1 in increments of 0.1
plt.yticks(np.arange(0.8, 1.01, 0.05))

# Set X-axis breaks from 5 to 20 in increments of 1
plt.xticks(np.arange(1, 21, 1))

# Set plot labels and title
plt.xlabel('Number of CpG Sites Selected')
plt.ylabel('Specificity')
plt.title('Feature Selection Analysis - TB Classifier [Multivariate Regression - Best Classifier]\n\nSpecificity')
plt.legend()

# Set white background
plt.gca().set_facecolor('white')

# Show the plot
plt.grid(True)
spec_lineplot = plt.gcf()
plt.show()
spec_lineplot.savefig("./Results/Spec_evaluation_validation_lineplot.png", dpi=600)

#%%

## AUC plots
## Evaluate performance of selected features
selected_dmc_indices = np.array([""])
selected_dmc_names = np.array(dmc_input.index[selected_dmc_indices])

# Define evaluation sets
test_sets = [
    ("Train Set", np.array(np.transpose(trainval_set)), np.transpose(np.array(trainval_meta.loc[:, "Group_bin"])))
    ("Test Set", np.array(np.transpose(test_set)), np.transpose(np.array(test_meta.loc[:, "Group_bin"])))
]

# Calculate AUC for specific classifier
classifiers = [
    ("Random Forest", RandomForestClassifier(n_estimators=100, max_depth=10, random_state=7))
]

# Train and evaluate classifiers
for clf_name, clf in classifiers:
    clf.fit(X_train[:, selected_dmc_indices], y_train)
    y_pred_prob = clf.predict_proba(X_test[:, selected_dmc_indices])[:, 1]
    y_pred = (y_pred_prob > 0.5).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
    specificity = tn / (tn + fp)
    sensitivity = tp / (tp + fn)
    auc_out = roc_auc_score(y_test, y_pred_prob)

    print(f"Classifier: {clf_name}")
    print(f"Specificity: {specificity:.2f}")
    print(f"Sensitivity: {sensitivity:.2f}")
    print(f"AUC: {auc_out:.2f}")
    print("--------------------------")

# Plot AUC
test_set_name, test_data, test_labels = test_sets[0]
boot_indices = np.random.choice(range(len(test_data)), size=len(test_data), replace=True)
boot_data = test_data[boot_indices]
boot_labels = test_labels[boot_indices]

y_pred_prob = clf.predict_proba(boot_data[:, selected_dmc_indices])[:, 1]
fpr, tpr, thresholds = roc_curve(boot_labels, y_pred_prob)
roc_auc = auc(fpr, tpr)

f = plt.figure(figsize=(8, 6))
#plt.plot(fpr, tpr, color='#236AB9', lw=2, label="ROC Curve Training Cohort (area = %0.2f)" % auc(fpr, tpr))
plt.plot(fpr, tpr, color='#FC7307', lw=2, label="ROC Curve Validation Cohort (area = %0.2f)" % auc(fpr, tpr))
plt.plot([0, 1], [0, 1], color='grey', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc='lower right')
plt.grid(True)
plt.show()
f.savefig(main_dir + "Results/AUC_RF_validation.pdf")

