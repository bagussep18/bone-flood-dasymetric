# bone-flood-dasymetric

Code and data for a dasymetric refinement of flood exposure mapping in Bone Regency, South Sulawesi, Indonesia — comparing choropleth versus dasymetric approaches across 372 administrative villages.

---

## Repository Structure

```
bone-flood-dasymetric/
│
├── README.md
│
├── scripts/
│   ├── bone_gee_v11_vulnerable.js     # Google Earth Engine — dasymetric pipeline
│   └── statistical_analysis.py        # Python — statistical validation & SSEDI
│
└── data/
    └── bone_choropleth_dasymetric_vulnerable.csv  # Main dataset — 372 villages, 37 variables
```

---

## Scripts

### 1. `bone_gee_v11_vulnerable.js` — Google Earth Engine Pipeline

**Platform:** [Google Earth Engine](https://earthengine.google.com/) (free academic account required)

**What it does:**
- Loads BPS SP2025 census data for 372 desa/kelurahan in Bone Regency
- Performs dasymetric disaggregation anchored to **Google Open Buildings v3** footprints, cross-validated against the RDTR spatial planning dataset
- Computes GFI flood hazard classification from SRTM DEM: `GFI = ln(A/d²)`
- Overlays GFI hazard class masks (Low / Moderate / High) onto four population grids:
  - Total population
  - Umur Rentan (Vulnerable Age Group: children < 15 + elderly ≥ 65)
  - Miskin (Poor Households — BPS P3KE classification)
  - Disabilitas (Persons with Disability)
- Runs `reduceRegions()` at `tileScale=8` to compute pixel-level exposure per village per class
- Exports results as CSV

**How to run:**
1. Open [code.earthengine.google.com](https://code.earthengine.google.com)
2. Create a new script and paste the contents of `bone_gee_v11_vulnerable.js`
3. Update the asset paths at the top of the script to point to your uploaded assets
4. Click **Run** → export task will appear in the Tasks panel
5. Click **Run** on the export task to save to Google Drive

---

### 2. `statistical_analysis.py` — Statistical Validation

**Requirements:**
```
pip install pandas numpy scipy matplotlib seaborn
```

**What it does:**
- Wilcoxon signed-rank test (choropleth vs dasymetric per GFI class)
- Mann-Whitney U test (High class distributions)
- Spearman rank correlation (flood fraction vs MAUP error)
- OLS regression (spatial predictors of MAUP error)
- Social-Spatial Exposure Disparity Index (SSEDI) computation for three sub-groups
- Generates statistical figures

**How to run:**
```bash
python statistical_analysis.py bone_choropleth_dasymetric_vulnerable.csv
```

**Outputs:**
- `fig_statistics.png` — 4-panel statistical comparison
- `fig_vulnerable.png` — SSEDI by sub-group
- `stats_results.csv` — Full statistical test results

---

## Data

### `bone_choropleth_dasymetric_vulnerable.csv`

Main dataset exported from Google Earth Engine. Contains choropleth and dasymetric flood exposure estimates for all 372 administrative villages in Bone Regency.

**Dimensions:** 372 rows × 37 columns

**Key columns:**

| Column | Description |
|---|---|
| `nama_desa` | Village name |
| `kode_desa` | Village code (BPS) |
| `choro_R` / `dasy_R` | Choropleth / dasymetric exposure — Low hazard (persons) |
| `choro_S` / `dasy_S` | Choropleth / dasymetric exposure — Moderate hazard (persons) |
| `choro_T` / `dasy_T` | Choropleth / dasymetric exposure — High hazard (persons) |
| `ur_R/S/T` | Vulnerable Age Group exposure by class (dasymetric) |
| `mis_R/S/T` | Poor Households exposure by class (dasymetric) |
| `disab_R/S/T` | Persons with Disability exposure by class (dasymetric) |
| `flood_frac` | Proportion of village area in High hazard zone |
| `gap_jiwa` | MAUP gap: dasymetric − choropleth (High class, persons) |
| `gap_pct` | Relative MAUP gap (%) |
| `total_pop` | Total village population (BPS SP2025) |

---

## Study Area

**Bone Regency (Kabupaten Bone)**, South Sulawesi, Indonesia

- 372 administrative villages (desa/kelurahan)
- Total population: 838,161 persons (BPS SP2025)
- Coordinate system: UTM Zone 51S (EPSG:32751)
- Raster resolution: 30 metres

**GFI hazard classification:**

| Class | GFI range | Area (ha) | % of flood extent |
|---|---|---|---|
| Low | 0.000–0.333 | 41,664 | 45.3% |
| Moderate | 0.333–0.666 | 36,827 | 40.1% |
| High | 0.666–1.000 | 13,465 | 14.6% |

---

## Data Sources

| Dataset | Source | Access |
|---|---|---|
| Population census SP2025 | BPS Kabupaten Bone | Public |
| Building footprints | Google Open Buildings v3 | [sites.research.google/open-buildings](https://sites.research.google/open-buildings/) |
| SRTM DEM (30m) | NASA / USGS via GEE | `USGS/SRTMGL1_003` |
| GFI flood hazard | Derived (this study) | This repository |
| RDTR spatial plan | Dinas PUPR Bone Regency | Restricted — available on request |
| Sub-group population rasters | Derived (this study) | Available on request |

---

## License

Code: [MIT License](LICENSE)  
Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

---

## Contact

Open an issue in this repository for questions about the GEE pipeline or data.
