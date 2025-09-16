import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.metrics import confusion_matrix, classification_report
from tqdm import tqdm
import joblib
import gzip
import os

# --- Step 1: Load VCF Data and Process Genotypes ---
def parse_vcf(vcf_file):
    """
    Parse a VCF file (plain or gzipped) to extract CHROM, POS and genotype calls.
    Genotypes are converted as follows:
        - 0/0 or 0|0 -> 0
        - 0/1, 1/0, 0|1 or 1|0 -> 1
        - 1/1 or 1|1 -> 2
    """
    # Choose the correct open function based on file extension.
    open_func = gzip.open if vcf_file.endswith('.gz') else open

    # Read the file to extract the header line with sample names.
    with open_func(vcf_file, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                header = line.strip().split("\t")
                break
    # In a VCF, the first 9 columns are fixed; sample names start at column 10.
    sample_names = header[9:]
    
    # Determine the compression type for pandas.read_csv
    compression = 'gzip' if vcf_file.endswith('.gz') else None
    
    # Load the VCF data into a DataFrame (skip header lines starting with '#')
    df = pd.read_csv(vcf_file, sep="\t", comment="#", header=None, compression=compression)
    # Extract CHROM and POS columns.
    chrom_pos_df = df.iloc[:, [0, 1]]
    # Genotype data starts at column index 9.
    genotype_df = df.iloc[:, 9:].copy()
    genotype_df.columns = sample_names

    def convert_genotype(gt_str):
        """
        Convert genotype string to numeric code.
        """
        if pd.isnull(gt_str):
            return np.nan
        # Genotype field is usually the first field before a colon.
        gt = gt_str.split(":")[0]
        if gt in {"./.", ".|."}:
            return np.nan
        if gt in ["0/0", "0|0"]:
            return 0
        elif gt in ["0/1", "1/0", "0|1", "1|0"]:
            return 1
        elif gt in ["1/1", "1|1"]:
            return 2
        else:
            return np.nan

    # Apply the conversion for every cell in the genotype DataFrame.
    genotype_df = genotype_df.applymap(convert_genotype)
 
    # Combine the chromosome/position data with the converted genotype calls.
    df_vcf = pd.concat([chrom_pos_df, genotype_df], axis=1)
    return df_vcf, sample_names
print("Loading VCF file...")
vcf_file = "filtered_snps.vcf.gz"  # change this to your VCF file path
df_vcf, sample_names = parse_vcf(vcf_file)

# Extract cultivar labels from sample names (assuming cultivar is the prefix before an underscore)
cultivar_labels = [name.split('_')[0].rstrip('.') for name in sample_names]
sample_to_cultivar = dict(zip(sample_names, cultivar_labels))

# Prepare SNP dataset from the VCF data:
chrom_positions = df_vcf.iloc[:, :2]  # CHROM and POS columns
snp_data = df_vcf.iloc[:, 2:]
snp_data.columns = sample_names  # Ensure correct sample names

print("Loading hotspots...")
hotspots = pd.read_csv("filtered_regions.bed", sep="\t", header=None, names=["Chromosome", "Start", "End","Values","Flanks"])
print(hotspots)
# --- Step 2: Compute CGR Encoding ---
def get_final_cgr_position(snp_vector, grid_size=8):
    x, y = grid_size // 2, grid_size // 2
    for snp in snp_vector:
        if snp == 0:
            x, y = x // 2, y // 2
        elif snp == 1:
            x, y = (x + grid_size) // 2, (y + grid_size) // 2
        elif snp == 2:
            x, y = (x + grid_size) // 2, y // 2
    distance = np.sqrt((x - (grid_size // 2))**2 + (y - (grid_size // 2))**2)
    return x, y, distance

print("Generating CGR features...")
feature_dict, x_positions, y_positions = {}, {}, {}
for sample in sample_names:
    feature_dict[sample], x_positions[sample], y_positions[sample] = [], [], []

tqdm_hotspots = tqdm(hotspots.iterrows(), total=len(hotspots), desc="Processing Hotspots")
for idx, hotspot in tqdm_hotspots:
    
    chrom, start, end = hotspot["Chromosome"], hotspot["Start"], hotspot["End"]
    print(chrom)
    hotspot_snps = chrom_positions[
        (chrom_positions[0] == chrom) &
        (chrom_positions[1] >= (start + 1)) &  # Adjusting from 0-based to 1-based
        (chrom_positions[1] <= end)
    ].index

    hotspot_snp_data = snp_data.iloc[hotspot_snps]
    for sample in sample_names:
        # Get the SNP vector for this hotspot for the current sample
        snp_vector = hotspot_snp_data[sample].values if sample in hotspot_snp_data.columns else []
        final_x, final_y, final_distance = get_final_cgr_position(snp_vector)
        feature_dict[sample].append(final_distance)
        x_positions[sample].append(final_x)
        y_positions[sample].append(final_y)

feature_matrix = pd.DataFrame.from_dict(feature_dict, orient="index")
feature_matrix["Cultivar"] = [sample_to_cultivar[s] for s in feature_matrix.index]
feature_matrix.to_csv("hotspot_cgr_encoded_matrix.txt", sep="\t", index=False)
print("Feature matrix saved.")
cultivar_counts = feature_matrix["Cultivar"].value_counts()
print("Cultivar counts before exclusion:")
print(cultivar_counts)
common_cultivars = cultivar_counts[cultivar_counts >= 3].index
# Filter the feature matrix to only include samples from common cultivars
feature_matrix_filtered = feature_matrix[feature_matrix["Cultivar"].isin(common_cultivars)]
print("Cultivar counts after exclusion:")
print(feature_matrix_filtered["Cultivar"].value_counts())

# --- Step 3: Train Classifier ---
print("Training classifier...")
# Prepare data for training
X = feature_matrix_filtered.drop(columns=["Cultivar"]).values
y = feature_matrix_filtered["Cultivar"].values


rare_cultivars = cultivar_counts[cultivar_counts < 3]
print("Cultivars with fewer than 3 samples:")
print(rare_cultivars)
param_grid = {"n_estimators": [50, 100, 200, 300, 400, 500, 600]}
model = RandomForestClassifier(random_state=42)

outer_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=42)

nested_scores = []
for train_idx, test_idx in outer_cv.split(X, y):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    grid_search = GridSearchCV(model, param_grid, cv=inner_cv, scoring='accuracy')
    grid_search.fit(X_train, y_train)
    best_model = grid_search.best_estimator_

    test_score = best_model.score(X_test, y_test)
    print(test_score)
    nested_scores.append(test_score)

mean_nested_accuracy = np.mean(nested_scores)
std_nested_accuracy = np.std(nested_scores)
print(f"Nested Cross-Validation Accuracy: {mean_nested_accuracy:.4f} (+/- {std_nested_accuracy:.4f})")

# --- Step 4: Retrain on Whole Dataset with New Grid Search ---
print("Performing new grid search on full dataset...")
grid_search_final = GridSearchCV(model, param_grid, cv=3, scoring='accuracy', n_jobs=-1)
grid_search_final.fit(X, y)
best_model_final = grid_search_final.best_estimator_

# Extract feature importance
feature_importances = pd.DataFrame({
    "Feature": range(X.shape[1]),
    "Importance": best_model_final.feature_importances_
})
feature_importances = feature_importances.sort_values(by="Importance", ascending=False)

# Select Top 30 Features
top_30_features = feature_importances.iloc[:15]["Feature"].values
X_top_30 = X[:, top_30_features]
print(top_30_features)

# Retrain model with top 30 features
print("Retraining with top 30 features...")
best_model_final.fit(X_top_30, y)
y_pred = best_model_final.predict(X_top_30)

# Compute and display Confusion Matrix
conf_matrix = confusion_matrix(y, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=np.unique(y), yticklabels=np.unique(y))
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix - Top 30 Features")
plt.savefig("confusion.png")

# Print Classification Report
print("Classification Report:")
print(classification_report(y, y_pred))

# --- Step 5: Map Top 30 Features to Hotspots ---
print("Mapping top 30 features to their respective hotspots...")
top_30_hotspots = hotspots.iloc[top_30_features].copy()
top_30_hotspots.loc[:, "Start"] = top_30_hotspots["Start"]
top_30_hotspots = top_30_hotspots[['Chromosome', 'Start', 'End']]
top_30_hotspots.to_csv("top_15_hotspots.bed", sep="\t", header=False, index=False)
print("Top 30 hotspots saved in BED format as 'top_30_hotspots.bed'.")

# Save the trained model with top 30 features
print("Saving trained model with top 30 features...")
joblib.dump(best_model_final, "best_model_top15.pkl")
print("Model saved as 'best_model_top30.pkl'.")
