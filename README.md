# Estimasi Debit Sungai Menggunakan Metode CMHC Berbasis Citra Sentinel-2

**Penulis:** Anissa Zuhrita  

## 📌 Deskripsi Proyek
Repositori ini berisi kode sumber untuk penelitian Tesis mengenai estimasi debit sungai pada skala segmen menggunakan metode *Calibration Measurement Hierarchical Classification* (CMHC). Studi kasus dilakukan di tiga sungai dengan karakteristik morfologi berbeda: Sungai Progo, Sungai Oyo, dan Sungai Serang.

## 📂 Struktur Kode

Proyek ini terdiri dari 4 langkah utama:

### 1. Langkah A (Google Earth Engine)
- **File:** `01_Langkah_A_Ekstraksi_Tanggal.js`
- **Fungsi:** Melakukan seleksi citra Sentinel-2 secara ketat (*Strict Filtering*) berdasarkan tutupan awan dan menentukan label musim (Hujan/Kemarau) berdasarkan data ZOM historis.

### 2. Langkah B (Python / Google Colab)
- **File:** `02_Langkah_B_Preprocessing_Jenks.ipynb`
- **Fungsi:** Membersihkan data debit lapangan (outlier removal), melakukan klasifikasi debit menjadi 3 kelas (Rendah, Menengah, Tinggi) menggunakan metode *Jenks Natural Breaks*, dan visualisasi distribusi data.

### 3. Langkah C (Google Earth Engine)
- **File:** `03_Langkah_C_Ekstraksi_Fitur.js`
- **Fungsi:** Mengekstraksi fitur spektral (Reflektansi C, M, W, dan rasio) dari citra Sentinel-2. Menggunakan algoritma otomatis untuk memilih piksel 'M' terbaik berdasarkan korelasi Spearman tertinggi.

### 4. Langkah D (Python / Google Colab)
- **File:** `04_Langkah_D_Pemodelan_Evaluasi.ipynb`
- **Fungsi:**
  - Melatih model klasifikasi *Random Forest* (RF).
  - Membangun model regresi hierarkis (Cubic/Power/Exponential) untuk setiap kelas debit.
  - Evaluasi akurasi menggunakan NSE, RMSE, dan RRMSE.
  - Visualisasi Hidrograf dan Scatter Plot.

## 🛠️ Cara Menggunakan
1. **Langkah A & C:** Copy kode ke [Google Earth Engine Code Editor](https://code.earthengine.google.com/). Pastikan asset Shapefile ROI sudah diimport.
2. **Langkah B & D:** Buka file `.ipynb` menggunakan Google Colab atau Jupyter Notebook. Sesuaikan path Google Drive (`GEE_EXPORT_FOLDER`) dengan lokasi data Anda.

