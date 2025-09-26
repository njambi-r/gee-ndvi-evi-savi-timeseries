# NDVI / EVI / SAVI Time Series with Google Earth Engine

This repository contains a Google Earth Engine (GEE) script for comparing **NDVI, EVI, and SAVI** across time.  
The workflow automates data preprocessing, vegetation index calculations, visualization, and GIF/chart generation for temporal vegetation analysis.

---

## 🌍 Project Context
- **Author**: Rose Njambi  
- **Date**: September 2025  
- **Blog Reference**: Developed for the [Geohub Blog](https://geohubkenya.wordpress.com/)) 

The script builds upon the work of [Source of Flavor](https://github.com/Source-of-Flavor/google-earth-engine/blob/main/NDVI_Time_Series.j), with modifications for multi-index comparisons, histograms, etc.

---

## 📌 Features
- Define **Area of Interest (AOI)** for forest, agriculture, or custom study areas.
- Load **Sentinel-2 Level-2A Harmonized Surface Reflectance** data.
- Apply **cloud masking** using:
  - Sentinel-2 Scene Classification Layer (SCL)
  - Sentinel-2 Cloud Probability dataset
- Generate **monthly composites** using median values.
- Compute vegetation indices:
  - **NDVI** – Normalized Difference Vegetation Index
  - **EVI** – Enhanced Vegetation Index
  - **SAVI** – Soil Adjusted Vegetation Index
- Normalize indices between **0–1** for comparability across months.
- Visualization outputs:
  - 📊 **Time series charts** for NDVI, EVI, SAVI, and reflectance bands (B4, B8).
  - 📉 **Image count** and **cloud contamination charts**.
  - 🖼️ **Histograms** of vegetation indices.
  - 🎞️ **Animated GIFs** of indices, RGB composites, and side-by-side comparisons.

---

## 🚀 Getting Started

### Prerequisites
- A [Google Earth Engine](https://earthengine.google.com/) account.
- Familiarity with GEE Code Editor (`code.earthengine.google.com`).

### Usage
1. Open the GEE **Code Editor**.
2. Copy-paste the script from [`ndvi_evi_savi_timeseries.js`](./ndvi_evi_savi_timeseries.js).
3. Define your AOI:
   ```javascript
   // Example AOI
   var aoi = Karura_outline;   // Forest
   // var aoi = Morendat;      // Agricultural
4. Run the script to:
   - Generate charts.
   - Print image counts.
   - Export images to Google Drive.
   - Create GIF animations.

## 📊 Example Outputs
### NDVI vs. EVI Time Series
<img width="1326" height="541" alt="ee-chart_NDVI and EVI Monthly Time Series 2024" src="https://github.com/user-attachments/assets/202ecff1-aad7-41cc-b28d-e7faab86712d" />

### Monthly Image Count
<img width="1326" height="541" alt="ee-chart_Monthly Image Count" src="https://github.com/user-attachments/assets/c36ccb5e-169a-4f8b-905e-5cfe2d4a734c" />

### Cloud/Shadow Contamination (%)
<img width="1326" height="541" alt="Monthly Cloud_Shadow Contamination" src="https://github.com/user-attachments/assets/46d7af9c-59fe-4164-a4b4-cfba6d9a8b7b" />

### NDVI/EVI/SAVI GIF animations   
![Sat_NDVI_EVI_SAVI](https://github.com/user-attachments/assets/e39f2b19-e109-4ddd-8aeb-3194ba449a1a)

## 🙌 Acknowledgments
Built with Google Earth Engine

## 📜 License
This project is open source. Please provide appropriate attribution when using or modifying this code.

## Contributing
Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests to improve the functionality and documentation.
