
# CPT-BASED LIQUEFACTION ASSESSMENT
# Boulanger & Idriss (2014) procedure applied to six sites
# from the 2011 Christchurch Earthquake
"""
Created on Mon Jul 27 2026

@author: Nyarko Kennedy
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# STAGE 0: LOAD DATA

df = pd.read_pickle("CANTERBURYDATASET.pkl")


# STAGE 1: SITE SELECTION
# Filter the dataset and select one CPT per manifestation class (0–5)


# Extract 2011 manifestation class (index 1 of [2010, 2011, 2016])
df["2011_Manifestation"] = df["Manifestation"].apply(lambda x: x[1])

# Remove sites with unknown manifestation (code 10)
df_unknown = df[df["2011_Manifestation"] != 10].copy()
df_unknown.loc[:, "2011_GWT"] = df_unknown["GWT"].apply(lambda x: x[1])

# Keep only sites where pre-drill is above the groundwater table
df_clean = df_unknown[df_unknown["pd"] < df_unknown["2011_GWT"]].copy()
df_clean.loc[:, "2011_PGA"] = df_clean["PGA"].apply(lambda x: x[1])

# Remove sites with very weak shaking
df_final = df_clean[df_clean["2011_PGA"] >= 0.075].copy()

def site_selection(manifestation_code):
    """This function takes the manifestation code(0-5) and randomly selects one CPT site from the data set
    with that same maifestation."""
    df_m = df_final[df_final["2011_Manifestation"] == manifestation_code].copy()
    df_m.loc[:, "Max_Depth"] = df_m["depth"].apply(lambda x: max(x))
    df_m = df_m[df_m["Max_Depth"] >= 6]          # need sufficient depth
    site = df_m.sample(n=1, random_state=42)     # reproducible selection
    return site

# One site for each manifestation class
site_0 = site_selection(0)
site_1 = site_selection(1)
site_2 = site_selection(2)
site_3 = site_selection(3)
site_4 = site_selection(4)
site_5 = site_selection(5)

sites = [site_0, site_1, site_2, site_3, site_4, site_5]
classes = [0, 1, 2, 3, 4, 5]


# STAGE 2: STRESS PROFILE + SOIL BEHAVIOUR TYPE INDEX (Ic)


def analyze_Ic_profile(site, qc_col="qc", fs_col="fs", a=0.80, tolerance=0.001):
    """
    This function calculates both total and effective stresses and iteratively compute Robertson Ic.
    It works for both raw/measured (qc, fs) and thin-layer/corrected (qc_inv, fs_inv) data.
    """
    depth = site["depth"].iloc[0]
    qc = site[qc_col].iloc[0]
    fs = site[fs_col].iloc[0]
    u2 = site["u2"].iloc[0]
    pd_depth = site["pd"].iloc[0]
    gwt = site["2011_GWT"].iloc[0]

    # Corrected tip resistance (unequal end area)
    qt = qc + u2 * (1 - a)

    # Friction ratio for unit-weight correlation
    Rf = (fs / qt) * 100

    gamma_w = 9.81          # kN/m³
    Pa = 100                # kPa

    # Unit weight (Robertson-style correlation) 
    gamma = gamma_w * (0.27 * np.log10(Rf) + 0.36 * np.log10(qt / Pa) + 1.236)

    # Total vertical stress
    delta_z = np.diff(depth)
    layer_stress = gamma[:-1] * delta_z
    sigma_v = np.insert(np.cumsum(layer_stress), 0, 0)

    # Pore pressure (hydrostatic below GWT)
    u = np.where(depth > gwt, gamma_w * (depth - gwt), 0)
    sigma_v_eff = sigma_v - u

    # Protect against zero/negative effective stress
    sigma_v_eff = np.maximum(sigma_v_eff, 1.0)
    qt = np.maximum(qt, 1.0)

    # Determine reliable start depth (below pre-drill and first non-zero qc)
    reliable_start = np.where(depth > pd_depth)[0][0]
    qc_diff = np.diff(qc)
    true_reliable_start = np.where(qc_diff[reliable_start:] != 0)[0][0]
    true_start = reliable_start + true_reliable_start

    # Iterative calculation of Ic and n
    capped_n = np.ones(len(depth))
    max_difference = 0.5
    max_iterations = 20
    iteration_count = 0

    while max_difference > tolerance and iteration_count < max_iterations:
        n_previous = capped_n.copy()

        Qtn = ((qt - sigma_v) / Pa) * (Pa / sigma_v_eff) ** capped_n
        Fr = (fs / (qt - sigma_v)) * 100

        # Numerical safeguards
        Qtn = np.maximum(Qtn, 0.1)
        Fr = np.maximum(Fr, 0.01)

        Ic = np.sqrt((3.47 - np.log10(Qtn))**2 + (np.log10(Fr) + 1.22)**2)

        n_new = 0.381 * Ic + 0.05 * (sigma_v_eff / Pa) - 0.15
        n_damped = 0.5 * n_previous + 0.5 * n_new
        capped_n = np.where(n_damped >= 1.0, 1.0, n_damped)

        difference = np.abs(capped_n[true_start:] - n_previous[true_start:])
        max_difference = np.max(difference)
        iteration_count += 1

    return {"Ic": Ic,
            "sigma_v": sigma_v,
            "sigma_v_eff": sigma_v_eff,
            "true_start": true_start,
            "Qtn": Qtn}


# STAGE 3: FINES CONTENT ESTIMATION


def calculate_fines_content(Ic, CFC=0.2):
    """This function estimate fines content from Ic (Christchurch-calibrated CFC)."""
    FC = 80 * (Ic + CFC) - 137
    FC = np.clip(FC, 0, 100)
    return FC


# STAGE 4: CYCLIC RESISTANCE RATIO (CRR) – Boulanger & Idriss 2014


def calculate_crr(qc, FC, sigma_v_eff, M, Pa=100, tolerance=0.05):
    """
    This function uses Boulanger & Idriss (2014) CRR calculation.
    It iteratively computes qc1Ncs from raw/measured qc, then CRR with 0.6 cap.
    """
    qc1Ncs = np.full_like(qc, 80.0)
    max_difference = 1.0
    max_iterations = 20
    iteration_count = 0

    while max_difference > tolerance and iteration_count < max_iterations:
        qc1Ncs_previous = qc1Ncs.copy()

        m = 1.338 - 0.249 * (qc1Ncs ** 0.264)
        m = np.maximum(m, 0.0)

        CN = (Pa / np.maximum(sigma_v_eff, 1.0)) ** m
        CN = np.minimum(CN, 1.7)

        qc1N = CN * (qc / Pa)

        delta = (11.9 + qc1N / 14.6) * np.exp(1.63 - 9.7/(FC + 2) - (15.7/(FC + 2))**2)
        qc1Ncs = qc1N + delta
        qc1Ncs = np.clip(qc1Ncs, 1.0, 300.0)

        difference = np.abs(qc1Ncs - qc1Ncs_previous)
        max_difference = np.nanmax(difference)
        iteration_count += 1

    # CRR at M=7.5, σ'v=1 atm (with mandatory 0.6 cap)
    CRR_7 = np.exp(
        (qc1Ncs/113) + (qc1Ncs/1000)**2
        - (qc1Ncs/140)**3 + (qc1Ncs/137)**4 - 2.8)
    CRR_7 = np.minimum(CRR_7, 0.6)

    # Magnitude Scaling Factor
    MSF_max = 1.09 + (qc1Ncs/180)**3
    MSF_max = np.where(MSF_max > 2.2, 2.2, MSF_max)
    MSF = 1 + (MSF_max - 1) * (8.64 * np.exp(-M/4) - 1.325)

    # Overburden correction Kσ
    C_sigma = 1 / (37.3 - 8.27 * np.minimum(qc1Ncs, 211)**0.264)
    C_sigma = np.minimum(C_sigma, 0.3)
    
    K_sigma = 1 - C_sigma * np.log(np.maximum(sigma_v_eff, 1.0) / Pa)
    K_sigma = np.minimum(K_sigma, 1.1)

    CRR = CRR_7 * MSF * K_sigma
    return CRR


# STAGE 5: CYCLIC STRESS RATIO (CSR)


def calculate_csr(depth, sigma_v, sigma_v_eff, PGA, M):
    """This function uses Boulanger & Idriss (2014) CSR with rd 
    dependent on magnitude and depth."""
    alpha = -1.012 - 1.126 * np.sin(depth/11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(depth/11.28 + 5.142)
    rd = np.exp(alpha + beta * M)

    CSR = 0.65 * PGA * (sigma_v / sigma_v_eff) * rd
    return CSR


# STAGE 6: RUN ANALYSIS ON ALL SIX SITES


results = []

for site, m_class in zip(sites, classes):
    M = site["Magnitude"].iloc[0][1]
    PGA = site["2011_PGA"].iloc[0]
    depth = site["depth"].iloc[0]

    # Measured (raw qc)
    res_m = analyze_Ic_profile(site)
    FC_m = calculate_fines_content(res_m["Ic"])
    CRR_m = calculate_crr(site["qc"].iloc[0], FC_m, res_m["sigma_v_eff"], M)
    CSR_m = calculate_csr(depth, res_m["sigma_v"], res_m["sigma_v_eff"], PGA, M)
    FS_m = CRR_m / CSR_m
    FS_m = np.where(res_m["Ic"] > 2.6, np.nan, FS_m)

    # Thin-layer corrected (qc_inv) 
    res_c = analyze_Ic_profile(site, qc_col="qc_inv", fs_col="fs_inv")
    FC_c = calculate_fines_content(res_c["Ic"])
    CRR_c = calculate_crr(site["qc_inv"].iloc[0], FC_c, res_c["sigma_v_eff"], M)
    CSR_c = calculate_csr(depth, res_c["sigma_v"], res_c["sigma_v_eff"], PGA, M)
    FS_c = CRR_c / CSR_c
    FS_c = np.where(res_c["Ic"] > 2.6, np.nan, FS_c)

    results.append({"Class": m_class,
                    "Min_FS_Measured": np.nanmin(FS_m),
                    "Min_FS_Corrected": np.nanmin(FS_c)})

summary = pd.DataFrame(results)
print("\n=== Summary of Minimum Factors of Safety ===")
print(summary.round(3))


# STAGE 7: PLOTTING


# Plot 1 – Grouped bar chart (main result)
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(classes))
width = 0.35

ax.bar(x - width/2, summary["Min_FS_Measured"], width,
       label="Measured (raw qc)", color="steelblue")
ax.bar(x + width/2, summary["Min_FS_Corrected"], width,
       label="Corrected (qc_inv)", color="darkorange")
ax.axhline(1.0, color="red", linestyle="--", linewidth=1.2, label="FS = 1.0")

ax.set_xlabel("Observed Manifestation Class")
ax.set_ylabel("Minimum Factor of Safety")
ax.set_title("Minimum FS vs Observed Liquefaction Severity\n(2011 Christchurch Earthquake)")
ax.set_xticks(x)
ax.set_xticklabels(classes)
ax.set_ylim(0, 0.6)
ax.legend()
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig("plot1_minFS_bar.png", dpi=300, bbox_inches="tight")
plt.show()

# Plot 2 – Trend view
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(classes, summary["Min_FS_Measured"], "o-", color="steelblue",
        label="Measured (raw qc)", markersize=8)
ax.plot(classes, summary["Min_FS_Corrected"], "s--", color="darkorange",
        label="Corrected (qc_inv)", markersize=8)
ax.axhline(1.0, color="red", linestyle="--", linewidth=1.2, label="FS = 1.0")

ax.set_xlabel("Observed Manifestation Class")
ax.set_ylabel("Minimum Factor of Safety")
ax.set_title("Trend of Minimum FS with Manifestation Severity")
ax.set_xticks(classes)
ax.set_ylim(0, 0.6)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plot2_minFS_trend.png", dpi=300, bbox_inches="tight")
plt.show()

# Plot 3 – FS profiles for Class 0 and Class 5

def get_fs_profile(site):
    M = site["Magnitude"].iloc[0][1]
    PGA = site["2011_PGA"].iloc[0]
    depth = site["depth"].iloc[0]

    # Measured CPT
    res_m = analyze_Ic_profile(site)
    FC_m = calculate_fines_content(res_m["Ic"])
    CRR_m = calculate_crr(site["qc"].iloc[0], FC_m, res_m["sigma_v_eff"], M)
    CSR_m = calculate_csr(depth, res_m["sigma_v"], res_m["sigma_v_eff"], PGA, M)
    FS_m = CRR_m / CSR_m
    FS_m = np.where(res_m["Ic"] > 2.6, np.nan, FS_m)

    # Thin-layer-corrected CPT
    res_c = analyze_Ic_profile(site, qc_col="qc_inv", fs_col="fs_inv")
    FC_c = calculate_fines_content(res_c["Ic"])
    CRR_c = calculate_crr(site["qc_inv"].iloc[0], FC_c, res_c["sigma_v_eff"], M)
    CSR_c = calculate_csr(depth, res_c["sigma_v"], res_c["sigma_v_eff"], PGA, M)
    FS_c = CRR_c / CSR_c
    FS_c = np.where(res_c["Ic"] > 2.6, np.nan, FS_c)

    return depth, FS_m, FS_c


depth0, FS0_m, FS0_c = get_fs_profile(site_0)
depth5, FS5_m, FS5_c = get_fs_profile(site_5)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6), sharey=True)

# Site Class 0
ax1.plot(FS0_m, depth0, color="steelblue", linewidth=1.5, linestyle="-", label="Measured FS")
ax1.plot(FS0_c, depth0, color="darkorange", linewidth=1.5, linestyle="--", label="Thin-layer-corrected FS")
ax1.axvline(1.0, color="black", linestyle="-.", linewidth=1.0, label="FS = 1.0")

ax1.set_xlabel("Factor of Safety")
ax1.set_ylabel("Depth (m)")
ax1.set_title("Site Class 0")
ax1.set_xlim(0, 1.5)
ax1.grid(True, alpha=0.3)
ax1.legend()

# Site Class 5
ax2.plot(FS5_m, depth5, color="steelblue", linewidth=1.5, linestyle="-", label="Measured FS")
ax2.plot(FS5_c, depth5, color="darkorange", linewidth=1.5, linestyle="--", label="Thin-layer-corrected FS")
ax2.axvline( 1.0, color="black", linestyle="-.", linewidth=1.0, label="FS = 1.0")

ax2.set_xlabel("Factor of Safety")
ax2.set_title("Site Class 5")
ax2.set_xlim(0, 1.5)
ax2.grid(True, alpha=0.3)
ax2.legend()


# Invert the shared depth axis once
ax1.invert_yaxis()

fig.suptitle("Factor of Safety Profiles – Class 0 vs Class 5", fontsize=13)
plt.tight_layout()
plt.savefig("plot3_FS_profiles.png", dpi=300, bbox_inches="tight")
plt.show()

# Plot 4 – FS profiles for all six sites
titles = ["Class 0 – No manifestation",
          "Class 1 – Minor",
          "Class 2 – Moderate",
          "Class 3 – Severe",
          "Class 4 – Lateral spreading",
          "Class 5 – Severe lateral spreading"]

fig, axes = plt.subplots(2, 3, figsize=(14, 9), sharey=True)
axes = axes.flatten()

for i, (site, title) in enumerate(zip(sites, titles)):
    M = site["Magnitude"].iloc[0][1]
    PGA = site["2011_PGA"].iloc[0]
    depth = site["depth"].iloc[0]

    # Measured CPT
    res_m = analyze_Ic_profile(site)
    FC_m = calculate_fines_content(res_m["Ic"])
    CRR_m = calculate_crr(site["qc"].iloc[0], FC_m, res_m["sigma_v_eff"], M)
    CSR_m = calculate_csr(depth, res_m["sigma_v"], res_m["sigma_v_eff"], PGA, M)
    FS_m = CRR_m / CSR_m
    FS_m = np.where(res_m["Ic"] > 2.6, np.nan, FS_m)

    # Thin-layer corrected CPT
    res_c = analyze_Ic_profile(site, qc_col="qc_inv", fs_col="fs_inv")
    FC_c = calculate_fines_content(res_c["Ic"])
    CRR_c = calculate_crr(site["qc_inv"].iloc[0], FC_c, res_c["sigma_v_eff"], M)
    CSR_c = calculate_csr(depth, res_c["sigma_v"], res_c["sigma_v_eff"], PGA, M)
    FS_c = CRR_c / CSR_c
    FS_c = np.where(res_c["Ic"] > 2.6, np.nan, FS_c)

    ax = axes[i]
    # Measured FS plot
    ax.plot(FS_m, depth, color="steelblue", linewidth=1.5, label="Measured FS")
    
    # Corrected FS plot
    ax.plot(FS_c, depth, color="darkorange", linewidth=1.5, linestyle="-.", label="Corrected FS")
    
    # FS =1.0 reference line
    ax.axvline(1.0, color="red", linestyle="-.", linewidth=1.0)
    
    
    ax.set_xlabel("Factor of Safety")
    if i % 3 == 0:
        ax.set_ylabel("Depth (m)")
    ax.set_title(title, fontsize=10)
    
    ax.set_xlim(0, 1.5)
    ax.grid(True, alpha=0.3)
axes[0].invert_yaxis()
fig.suptitle("Factor of Safety Profiles – All Six Sites (2011 Christchurch)", fontsize=13, y=1.02)

# One shared legend for the entire figure
handles, labels = axes[0].get_legend_handles_labels()

fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=3, frameon=False)

plt.tight_layout()
plt.savefig("FS_profiles_all_six_sites.png", dpi=300, bbox_inches="tight")
plt.show()

print("\nAll analyses and plots completed successfully.")