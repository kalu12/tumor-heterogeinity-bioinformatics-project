import pandas as pd
import vcf
import io
import os
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("cnv.txt", delimiter="\t")

# Izdvajanje zadatih kolona
selected_rows = df[(df["cn"] == 2) & (df["cn1"] == 1) & (df["cn2"] == 1)]
selected_rows1 = df[
    (df["cn"] == 2)
    & (
        (df["cn1"] == 1)
        | (df["cn2"] == 1)
        | (pd.isna(df["cn1"])) & (pd.isna(df["cn2"]))
    )
]


regioni = ["chromosome", "start", "end", "cn", "cn1", "cn2"]
regionidf = selected_rows1[regioni]

print(regionidf)


def read_vcf(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})


vcf_file = read_vcf("Strelka.vcf")
vcf_file1 = read_vcf("Strelka.vcf")
vcf_file["VAF"] = 0
vcf_file["Base number"] = 0


for index, row in vcf_file.iterrows():
    # Pristupanje vrednostima u svakom redu
    chromosome = row["CHROM"]
    position = row["POS"]

    # Iteriraj kroz redove
    match_found = False
    for index1, row1 in regionidf.iterrows():
        chromosome1 = row1["chromosome"]

        # Provera da li se vrednosti poklapaju
        if chromosome == chromosome1 and row1["start"] <= position <= row1["end"]:
            match_found = True

            # REF i ALT
            ref_value = row["REF"]
            alt_value = row["ALT"]

            if ref_value == "A":
                ref_value = 5
            elif ref_value == "C":
                ref_value = 6
            elif ref_value == "G":
                ref_value = 7
            elif ref_value == "T":
                ref_value = 8

            if alt_value == "A":
                alt_value = 5
            elif alt_value == "C":
                alt_value = 6
            elif alt_value == "G":
                alt_value = 7
            elif alt_value == "T":
                alt_value = 8

            VAF = row["TUMOR"]
            j = 0
            REF = ""
            check = 0
            for i, char in enumerate(VAF):
                if check == 1:
                    break
                if char == ":":
                    j += 1
                if (ref_value - 1) == j:
                    for k in range(i + 1, len(VAF)):
                        if VAF[k] == "," or VAF[k] == ":":
                            check = 1
                            break
                        REF += VAF[k]

            REF = int(REF)
            j = 0
            ALT = ""
            check = 0
            for i, char in enumerate(VAF):
                if check == 1:
                    break
                if char == ":":
                    j += 1
                if (alt_value - 1) == j:
                    for k in range(i + 1, len(VAF)):
                        if VAF[k] == "," or VAF[k] == ":":
                            check = 1
                            break
                        ALT += VAF[k]

            ALT = int(ALT)
            vcf_file.at[index, "Base number"] = REF + ALT

            Variant_Allele_Frequency = ALT / (ALT + REF)
            vcf_file.at[index, "VAF"] = Variant_Allele_Frequency

            if Variant_Allele_Frequency <= 0.6:
                match_found2 = True

            break

    # Dropovanje reda ako se nisu poklapale vrednosti
    if not match_found or not match_found2:
        vcf_file.drop(index, inplace=True)

import hdbscan


X = vcf_file[["VAF", "Base number"]]

clusterer = hdbscan.HDBSCAN(
    min_cluster_size=15, cluster_selection_epsilon=0.5, gen_min_span_tree=False
)
cluster_labels = clusterer.fit_predict(X)


X["Cluster"] = cluster_labels


sns.scatterplot(x="VAF", y="Base number", hue="Cluster", data=X)


unique_clusters = X["Cluster"].unique()
for cluster in unique_clusters:
    if cluster != -1:
        cluster_center = X[X["Cluster"] == cluster][["VAF", "Base number"]].mean()
        plt.scatter(
            cluster_center["VAF"],
            cluster_center["Base number"],
            marker="X",
            s=100,
            color="black",
        )

plt.show()
