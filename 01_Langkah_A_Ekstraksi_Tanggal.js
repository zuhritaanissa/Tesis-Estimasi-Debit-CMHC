/**
 * ============================================================================================
 * JUDUL       : Pre-processing & Ekstraksi Tanggal Citra Sentinel-2 (Metode CMHC)
 * DESKRIPSI   : Script ini melakukan filtering ketat (Strict Quality Control) pada koleksi
 * Sentinel-2 Surface Reflectance untuk mendapatkan tanggal citra bebas awan
 * dan pelabelan musim presisi berdasarkan data historis ZOM.
 * ============================================================================================
 */

/* * --------------------------------------------------------------------------------------------
 * [PENTING] INSTRUKSI PENGGUNAAN:
 * Sebelum menjalankan script ini, pastikan Anda telah mengimpor geometri/shapefile
 * batas daerah kajian (ROI) ke dalam script editor GEE dengan nama variabel:
 * 1. kalibawang
 * 2. kedungmiri
 * 3. bendungan
 * --------------------------------------------------------------------------------------------
 */

/* ========== 1. KONFIGURASI PARAMETER & DATABASE ========== */

// Periode Penelitian
var START_FULL = '2019-01-01';
var END_FULL   = '2024-12-31';

// Daftar Lokasi Stasiun (ROI)
var STATIONS = {
  'Kalibawang': kalibawang, 
  'Kedungmiri': kedungmiri,
  'Bendungan' : bendungan
};

// Parameter Filter Kualitas (Strict Mode)
// Tujuannya meminimalkan noise atmosfer pada analisis debit
var CLOUD_PROBABILITY_THRESHOLD = 20;  // Batas toleransi probabilitas awan (%)
var VALID_PIXEL_THRESHOLD       = 95;  // Minimum piksel valid dalam ROI (%)
var FINAL_CLOUD_COVER_PERCENT   = 10;  // Batas akhir tutupan awan yang diizinkan (%)

// Database Musim Hujan Presisi (Diselaraskan dengan Grafik Curah Hujan & ZOM)
var rainySeasonDatabase = [
    // --- KALIBAWANG (ZOM_01) ---
    { 'station':'Kalibawang', 'start':'2019-11-01', 'end':'2020-05-20' }, 
    { 'station':'Kalibawang', 'start':'2020-10-11', 'end':'2021-06-10' }, 
    { 'station':'Kalibawang', 'start':'2021-10-11', 'end':'2022-06-10' }, 
    { 'station':'Kalibawang', 'start':'2022-10-01', 'end':'2023-05-20' },
    { 'station':'Kalibawang', 'start':'2023-11-11', 'end':'2024-05-20' }, 
    { 'station':'Kalibawang', 'start':'2024-10-11', 'end':'2024-12-31' },

    // --- BENDUNGAN (ZOM_06) ---
    { 'station':'Bendungan', 'start':'2019-10-21', 'end':'2020-05-01' }, 
    { 'station':'Bendungan', 'start':'2020-10-11', 'end':'2021-04-30' }, 
    { 'station':'Bendungan', 'start':'2021-10-11', 'end':'2022-05-20' }, 
    { 'station':'Bendungan', 'start':'2022-10-01', 'end':'2023-06-30' },
    { 'station':'Bendungan', 'start':'2023-11-11', 'end':'2024-05-01' }, 
    { 'station':'Bendungan', 'start':'2024-10-21', 'end':'2024-12-31' },

    // --- KEDUNGMIRI (ZOM_07) ---
    { 'station':'Kedungmiri', 'start':'2019-11-01', 'end':'2020-04-30' }, 
    { 'station':'Kedungmiri', 'start':'2020-10-11', 'end':'2021-04-20' }, 
    { 'station':'Kedungmiri', 'start':'2021-10-11', 'end':'2022-05-01' }, 
    { 'station':'Kedungmiri', 'start':'2022-10-01', 'end':'2023-06-30' },
    { 'station':'Kedungmiri', 'start':'2023-11-01', 'end':'2024-05-01' }, 
    { 'station':'Kedungmiri', 'start':'2024-11-01', 'end':'2024-12-31' } 
];

/* ========== 2. KOLEKSI CITRA ========== */

var S2_SR    = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var S2_CLOUD = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');

/* ========== 3. FUNGSI UTAMA ========== */

/**
 * Fungsi untuk menentukan label musim berdasarkan tanggal citra dan database.
 * @param {String} dateStr - Tanggal citra (YYYY-MM-DD)
 * @param {String} stationName - Nama stasiun
 * @param {List} seasonDB - Database musim
 * @returns {String} 'Hujan' atau 'Kemarau'
 */
var determineSeason = function(dateStr, stationName, seasonDB) {
  var date = ee.Date.parse('YYYY-MM-dd', dateStr);
  
  var isRainy = ee.List(seasonDB).map(function(record){
    var rec = ee.Dictionary(record);
    var isSameStation = ee.String(rec.get('station')).compareTo(stationName).eq(0);
    var start = ee.Date(rec.get('start'));
    var end   = ee.Date(rec.get('end'));
    
    // Cek apakah tanggal berada dalam rentang musim hujan
    var inRange = date.millis().gte(start.millis()).and(date.millis().lte(end.millis()));
    return isSameStation.and(inRange);
  }).reduce(ee.Reducer.max());
  
  return ee.Algorithms.If(isRainy, 'Hujan', 'Kemarau');
};

/**
 * Fungsi inti untuk memfilter citra berdasarkan kualitas piksel dan menghitung statistik.
 */
function extractCleanDates(roiFC, start, end, roiName) {
  var geom = ee.FeatureCollection(roiFC).geometry();
  
  // 1. Ambil semua tanggal ketersediaan awal
  var dates = S2_SR.filterBounds(geom).filterDate(start, end)
    .aggregate_array('system:time_start')
    .map(function(t){ return ee.Date(t).format('YYYY-MM-dd'); })
    .distinct();

  // 2. Iterasi per tanggal untuk evaluasi kualitas
  var feats = ee.FeatureCollection(ee.List(dates).map(function(d) {
    var day  = ee.Date.parse('YYYY-MM-dd', d);
    var next = day.advance(1, 'day');
    
    // Mosaic harian (penting untuk menangani overlap orbit)
    var sr = S2_SR.filterBounds(geom).filterDate(day, next).mosaic();
    var cp = S2_CLOUD.filterBounds(geom).filterDate(day, next).mosaic();

    // --- Masking Awan (Cloud Probability) ---
    var cloudProbScore = ee.Algorithms.If(
      cp.bandNames().length().gt(0),
      cp.select('probability').gt(CLOUD_PROBABILITY_THRESHOLD), ee.Image(0));

    // --- Masking SCL (Scene Classification Map) ---
    // Kelas yang dihindari: 3 (Shadow), 8 (Medium), 9 (High), 10 (Cirrus)
    var sclBadParams = [3, 8, 9, 10];
    var sclBadPixel = ee.Algorithms.If(
      sr.bandNames().contains('SCL'),
      sr.select('SCL').remap(sclBadParams, [1, 1, 1, 1], 0), ee.Image(0));

    var combinedBadPixel = ee.Image(cloudProbScore).or(ee.Image(sclBadPixel));

    // --- Hitung Persentase Area ---
    // Area Buruk (Awan/Bayangan)
    var badArea = combinedBadPixel.multiply(ee.Image.pixelArea())
        .reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: 10, maxPixels: 1e10})
        .values().get(0);
    badArea = ee.Algorithms.If(badArea, badArea, 0); 
    var badPct = ee.Number(badArea).divide(geom.area({maxError:1})).multiply(100);

    // Area Valid (Data tersedia/tidak bolong)
    var validArea = ee.Algorithms.If(
      sr.bandNames().length().gt(0),
      sr.select('B4').mask().multiply(ee.Image.pixelArea())
        .reduceRegion({reducer: ee.Reducer.sum(), geometry: geom, scale: 10, maxPixels: 1e10})
        .get('B4'), 0);
    var validPct = ee.Number(validArea).divide(geom.area({maxError:1})).multiply(100);

    // --- Labeling Musim ---
    var seasonDB_EE = ee.List(rainySeasonDatabase);
    var seasonLabel = determineSeason(d, roiName, seasonDB_EE);

    return ee.Feature(null, {
      date: d,
      BAD_PIXEL_PERCENTAGE: badPct,
      VALID_PIXEL_PERCENTAGE_ROI:  validPct,
      STATION: roiName,
      season: seasonLabel
    });
  }));

  // 3. Terapkan Filter Threshold
  return feats
    .filter(ee.Filter.lt('BAD_PIXEL_PERCENTAGE', FINAL_CLOUD_COVER_PERCENT)) 
    .filter(ee.Filter.gt('VALID_PIXEL_PERCENTAGE_ROI',  VALID_PIXEL_THRESHOLD));
}

/* ========== 4. EKSEKUSI UTAMA ========== */

print('🚀 MEMULAI PROSES EKSTRAKSI DATA...');
print('----------------------------------------------------');

var stationNames = Object.keys(STATIONS);

stationNames.forEach(function(name) {
  var currentROI = STATIONS[name];
  
  // Proses Ekstraksi
  var resultFC = extractCleanDates(currentROI, START_FULL, END_FULL, name);
  
  // Tampilkan Ringkasan di Console
  var stats = resultFC.aggregate_histogram('season');
  var total = resultFC.size();
  
  total.evaluate(function(n) {
    stats.evaluate(function(s) {
      print('📍 STASIUN: ' + name.toUpperCase());
      print('   • Total Citra Bersih : ' + n);
      if (s) {
        print('   • 🌧️ Musim Hujan      : ' + (s['Hujan'] || 0));
        print('   • ☀️ Musim Kemarau    : ' + (s['Kemarau'] || 0));
      } else {
        print('   • (Tidak ada data lolos filter)');
      }
      print('----------------------------------------------------');
    });
  });
  
  // Ekspor Data ke Google Drive (CSV)
  var fileName = 'MASTER_TANGGAL_STRICT_' + name.toUpperCase() + 
                 '_' + START_FULL.slice(0,4) + '_' + END_FULL.slice(0,4);
                 
  Export.table.toDrive({
    collection: resultFC.select(['date', 'season', 'STATION']), 
    description: 'Ekspor_' + name,
    folder: 'GEE_Exports_Tesis_Master',
    fileNamePrefix: fileName,
    fileFormat: 'CSV',
    selectors: ['date', 'season', 'STATION']
  });
});

print('⏳ Silakan cek tab "Tasks" di sebelah kanan untuk menjalankan proses ekspor.');