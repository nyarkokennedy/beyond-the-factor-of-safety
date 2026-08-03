# Beyond the Factor of Safety

**A CPT-Based Investigation of Liquefaction Triggering and Surface Manifestation Using the 2011 Christchurch Earthquake Case Histories**

Python implementation of the Boulanger & Idriss (2014) CPT liquefaction-triggering procedure, applied to six real case histories from the 2011 Christchurch earthquake — one site for each observed surface manifestation severity (Class 0 = no damage → Class 5 = severe lateral spreading).

Independent technical exercise · July 2026  
Author: Nyarko Kennedy · Materials Laboratory, Ghana Highway Authority

---

### What this project does

- Implements the full Boulanger & Idriss (2014) CPT-based triggering workflow in pure Python (stress profile → Ic → fines content → qc1Ncs → CRR → CSR → Factor of Safety).
- Runs the analysis on both **as-measured** CPT data and **thin-layer-corrected** data (Boulanger & DeJong 2018 inverse filtering).
- Compares predicted minimum FS against the actual observed manifestation class at each site.

### Key finding

Every site (including the Class 0 “no manifestation” site) produced a minimum FS < 1.0.  
The calculated FS values do **not** rank in the same order as the real-world surface damage severity.

This is the deliberate, counter-intuitive result the project is built around — and the reason for the title.

Thin-layer correction changed the numerical FS values (sometimes substantially) but never moved any site across the FS = 1.0 threshold, and did not restore a monotonic severity ranking.

---

### Results summary (minimum FS)

| Class | Observed Manifestation       | Min FS (Measured) | Min FS (Corrected) |
|-------|------------------------------|-------------------|--------------------|
| 0     | No manifestation             | 0.324             | 0.349              |
| 1     | Minor                        | 0.310             | 0.334              |
| 2     | Moderate                     | 0.259             | 0.264              |
| 3     | Severe                       | 0.269             | 0.347              |
| 4     | Lateral spreading            | 0.424             | 0.450              |
| 5     | Severe lateral spreading     | 0.327             | 0.317              |

---

### Repository contents

- Full Python implementation (`liquefaction_assessment.py`)
- Site selection, stress calculation, Ic iteration, CRR/CSR, FS profiles
- Comparison of measured vs thin-layer-corrected CPT data
- Publication-quality plots

---

### Requirements

pandas  
numpy  
matplotlib  

The Canterbury CPT dataset (`CANTERBURYDATASET.pkl`) from Geyin et al. (2021) is required to re-run the analysis.

---

### How to run

python liquefaction_assessment.py

The script will:

1. Filter and select one CPT per manifestation class (0–5) with a fixed random seed for reproducibility  
2. Compute the full FS profiles for both measured and corrected data  
3. Print the summary table  
4. Generate and save the plots

---

### Citation / Context

- Method: Boulanger & Idriss (2014)  
- Thin-layer correction: Boulanger & DeJong (2018)  
- Case histories: Geyin et al. (2021) Canterbury CPT liquefaction dataset  

### Disclaimer

This is a focused demonstration of method implementation and engineering understanding. The six sites are treated as individual case studies, not a statistical sample.
