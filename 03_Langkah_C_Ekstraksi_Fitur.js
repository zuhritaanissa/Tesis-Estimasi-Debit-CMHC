/**
 * ============================================================================================
 * JUDUL       : Ekstraksi Fitur Spektral & Pembentukan Aset CMHC (Langkah C)
 * DESKRIPSI   : Script ini mengekstraksi nilai reflektansi piksel Kalibrasi (C), Air (W),
 * dan Pengukuran (M) dari citra Sentinel-2. Pemilihan piksel M dilakukan
 * secara otomatis berdasarkan korelasi Spearman tertinggi terhadap debit.
 * ============================================================================================
 */

/* ============================================================================================
 * [PRASYARAT UPLOAD ASSET KE GEE]
 * Sebelum menjalankan script, Anda WAJIB meng-upload file CSV berikut ke tab 'Assets' di GEE:
 *
 * 1. FILE TANGGAL LATIH (Dari Output Langkah B):
 * - Nama File Asli: tanggal_latih_kalibawang_2019_2021_FIX.csv
 * - Rename Asset ID di GEE menjadi: tanggal_latih_kalibawang_2019_2021_FIX
 * (Ulangi untuk Bendungan dan Kedungmiri)
 *
 * 2. FILE TANGGAL UJI (Dari Output Langkah B):
 * - Nama File Asli: tanggal_uji_kalibawang_2022_2024_FIX.csv
 * - Rename Asset ID di GEE menjadi: tanggal_uji_kalibawang_2022_2024_FIX
 * (Ulangi untuk Bendungan dan Kedungmiri)
 *
 * 3. FILE DEBIT PENUH (File Mentah Input Langkah B):
 * - Nama File Asli: debit_kalibawang_2019_2024.csv
 * - Rename Asset ID di GEE menjadi: debit_kalibawang_2019_2024
 * (Pastikan kolomnya: 'date' dan 'debit')
 * ============================================================================================
 */

/* ========== 1. KONFIGURASI UTAMA ========== */

var USERNAME     = 'zuhritaanissa';      // <--- GANTI DENGAN USERNAME GEE ANDA
var USE_SEASONAL = false;                // [PENTING] false = Tanpa Musim, true = Musiman (Hujan/Kemarau)
var DRIVE_FOLDER = 'GEE_Exports_Tesis_Master'; 

// Daftar Stasiun (Pastikan geometri ROI sudah diimport: kalibawang, bendungan, kedungmiri)
var STATION_LIST = ['Kalibawang', 'Bendungan', 'Kedungmiri'];

// Database Threshold Jenks (Hasil Output Langkah B - SINKRON)
var THRESHOLDS_DB = {
  'Kalibawang': {
    'Hujan':    {q1: 55.300, q2: 94.300},
    'Kemarau':  {q1: 13.700, q2: 33.400},
    'Gabungan': {q1: 38.100, q2: 81.600}
  },
  'Bendungan': {
    'Hujan':    {q1: 0.110, q2: 3.570},
    'Kemarau':  {q1: 0.990, q2: 2.290},
    'Gabungan': {q1: 0.990, q2: 2.730}
  },
  'Kedungmiri': {
    'Hujan':    {q1: 10.400, q2: 26.300},
    'Kemarau':  {q1: 5.020,  q2: 34.600},
    'Gabungan': {q1: 10.400, q2: 26.300}
  }
};

// Parameter Ekstraksi
var THRESH_KORIDOR = 0.30; // Frekuensi air minimum
var M_BUFFER_PIX   = 3;    // Buffer piksel (30m)

// Periode Data
var START_LATIH = '2019-01-01'; var END_LATIH = '2021-12-31';
var START_UJI   = '2022-01-01'; var END_UJI   = '2024-12-31';
var PERIODE_TAG_Latih = '2019_2021';
var PERIODE_TAG_Uji   = '2022_2024';

/* ========== 2. FUNGSI HELPER (BAND MATH & STATISTIK) ========== */

var safeReduceFirst = function(img, band, geom, scale){
  return ee.Algorithms.If(
    ee.Algorithms.IsEqual(geom, null),
    null,
    img.select(band).reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: geom,
      scale: scale,
      maxPixels: 1e9
    }).get(band)
  );
};

var hitungMW = function(M, W){
  var mIsNull = ee.Algorithms.IsEqual(M, null); 
  var wIsNull = ee.Algorithms.IsEqual(W, null);
  return ee.Algorithms.If(
    mIsNull, null, 
    ee.Algorithms.If(
      wIsNull, null, 
      ee.Algorithms.If(
        ee.Number(M).add(ee.Number(W)).neq(0), 
        ee.Number(M).subtract(ee.Number(W)).divide(ee.Number(M).add(ee.Number(W))), 
        null
      )
    )
  );
};

var hitungRasio = function(C, M){
  var cIsNull = ee.Algorithms.IsEqual(C, null); 
  var mIsNull = ee.Algorithms.IsEqual(M, null);
  return ee.Algorithms.If(
    cIsNull, null, 
    ee.Algorithms.If(
      mIsNull, null, 
      ee.Algorithms.If(
        ee.Number(M).neq(0), 
        ee.Number(C).divide(ee.Number(M)), 
        null
      )
    )
  );
};

var hitungAWEIsh = function(img) {
  return img.addBands(img.expression(
    'B2 + 2.5*B3 - 1.5*(B8+B11) - 0.25*B12', 
    {'B2':img.select('B2'),'B3':img.select('B3'),'B8':img.select('B8'),'B11':img.select('B11'),'B12':img.select('B12')}
  ).rename('AWEIsh'));
};

var buatMaskerHarian = function(img){ return img.select('AWEIsh').gt(0).rename('waterMask'); };

/* ========== 3. DATABASE MUSIM PRESISI ========== */

var rois = {'Bendungan': bendungan, 'Kalibawang': kalibawang, 'Kedungmiri': kedungmiri};
var roi_zom_map = {'Bendungan':'ZOM_06','Kalibawang':'ZOM_01','Kedungmiri':'ZOM_07'};

var rainySeasonDatabase = [
  { zom:'ZOM_01', station:'Kalibawang', start:'2019-11-01', end:'2020-05-20' },
  { zom:'ZOM_01', station:'Kalibawang', start:'2020-10-11', end:'2021-06-10' },
  { zom:'ZOM_01', station:'Kalibawang', start:'2021-10-11', end:'2022-06-10' },
  { zom:'ZOM_01', station:'Kalibawang', start:'2022-10-01', end:'2023-05-20' },
  { zom:'ZOM_01', station:'Kalibawang', start:'2023-11-11', end:'2024-05-20' },
  { zom:'ZOM_01', station:'Kalibawang', start:'2024-10-11', end:'2024-12-31' },
  { zom:'ZOM_06', station:'Bendungan', start:'2019-10-21', end:'2020-05-01' },
  { zom:'ZOM_06', station:'Bendungan', start:'2020-10-11', end:'2021-04-30' },
  { zom:'ZOM_06', station:'Bendungan', start:'2021-10-11', end:'2022-05-20' },
  { zom:'ZOM_06', station:'Bendungan', start:'2022-10-01', end:'2023-06-30' },
  { zom:'ZOM_06', station:'Bendungan', start:'2023-11-11', end:'2024-05-01' },
  { zom:'ZOM_06', station:'Bendungan', start:'2024-10-21', end:'2024-12-31' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2019-11-01', end:'2020-04-30' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2020-10-11', end:'2021-04-20' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2021-10-11', end:'2022-05-10' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2022-10-01', end:'2023-06-30' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2023-11-01', end:'2024-05-10' },
  { zom:'ZOM_07', station:'Kedungmiri', start:'2024-11-01', end:'2024-12-31' }
];

var rainySeasonFC = ee.FeatureCollection(rainySeasonDatabase.map(function(r){return ee.Feature(null,r);}));

var getRainyIntervals = function(zomCode, startStr, endStr) {
  var baseStart = ee.Date(startStr), baseEnd = ee.Date(endStr);
  return rainySeasonFC.filter(ee.Filter.eq('zom', zomCode)).map(function(f){
    var rs = ee.Date(f.get('start')), re = ee.Date(f.get('end'));
    var keep = re.millis().gte(baseStart.millis()).and(baseEnd.millis().gte(rs.millis()));
    return f.set('keep', keep);
  }).filter(ee.Filter.eq('keep', 1));
};

var seasonLabel = function(dateStr, zomCode, startStr, endStr) {
  var d = ee.Date(dateStr);
  var intervals = getRainyIntervals(zomCode, startStr, endStr);
  var isRainy = intervals.map(function(f){
    return f.set('hit', d.millis().gte(ee.Date(f.get('start')).millis()).and(ee.Date(f.get('end')).millis().gte(d.millis())));
  }).aggregate_max('hit');
  return ee.String(ee.Algorithms.If(isRainy, 'Hujan', 'Kemarau'));
};

/* ========== 4. FUNGSI PROSES PER STASIUN ========== */

function processStation(currentStationName) {
  print('==================================================');
  print('🔄 MEMPROSES STASIUN: ' + currentStationName.toUpperCase());

  var roi_lower = currentStationName.toLowerCase();
  var roi_upper = currentStationName.toUpperCase();
  var current_roi = rois[currentStationName];
  var current_zom = roi_zom_map[currentStationName];
  var current_thresh = THRESHOLDS_DB[currentStationName];

  // Path Aset (Sesuai instruksi upload di atas)
  var pathTglLatih = 'projects/ee-' + USERNAME + '/assets/tanggal_latih_' + roi_lower + '_' + PERIODE_TAG_Latih + '_FIX';
  var pathTglUji   = 'projects/ee-' + USERNAME + '/assets/tanggal_uji_' + roi_lower + '_' + PERIODE_TAG_Uji + '_FIX';
  var pathDebit    = 'projects/ee-' + USERNAME + '/assets/debit_' + roi_lower + '_2019_2024';
  var ASSET_PROJECT_PATH = 'projects/ee-' + USERNAME + '/assets';

  // Load Feature Collection
  var fcLatih = ee.FeatureCollection(pathTglLatih);
  var fcUji = ee.FeatureCollection(pathTglUji);
  var listLatih = fcLatih.aggregate_array('date');
  var listUji = fcUji.aggregate_array('date');

  // Load Debit
  var debitRaw = ee.FeatureCollection(pathDebit).filter(ee.Filter.notNull(['date','debit']));
  var debitList = debitRaw.toList(5000).map(function(f){
    var fe = ee.Feature(f); var q = fe.get('debit');
    q = ee.Algorithms.If(ee.Algorithms.ObjectType(q).equals('String'), ee.Number.parse(q), q);
    return fe.set('debit_num', ee.Number(q).max(0));
  });
  var debitDict = ee.Dictionary.fromLists(
    debitList.map(function(f){return ee.Date(ee.Feature(f).get('date')).format('YYYY-MM-dd')}), 
    debitList.map(function(f){return ee.Feature(f).get('debit_num')})
  );

  // Koleksi Citra S2
  var S2_SR = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
  var buatKoleksi = function(l) { 
    return ee.ImageCollection.fromImages(l.map(function(d){
      var s=ee.Date(d); 
      return S2_SR.filterBounds(current_roi).filterDate(s, s.advance(1,'day'))
                  .mosaic().set('date',d).set('system:time_start',s.millis());
    })).sort('system:time_start'); 
  };
  var kolLatih = buatKoleksi(listLatih);
  var kolUji = buatKoleksi(listUji);

  // Badan Air (Water Mask Frequency)
  var maps = kolLatih.map(hitungAWEIsh).map(buatMaskerHarian);
  var freq = maps.reduce(ee.Reducer.mean());
  var BADAN_AIR = freq.gte(THRESH_KORIDOR);
  Map.addLayer(BADAN_AIR.selfMask().clip(current_roi), {palette:['blue']}, 'Air Perm '+currentStationName, false);

  // --- Fungsi Inti: Pilih M Terbaik & Ekstraksi ---
  var defineAssetsForSeason = function(slug, collUnik, collLengkap, q1, q2) {
    print('   -> Memproses ' + slug + ' (' + currentStationName + ')');
    
    // Identifikasi C (Kalibrasi) dan W (Air)
    var air = BADAN_AIR.rename('air_'+slug);
    var dil = air.focalMax({radius: M_BUFFER_PIX, units:'pixels'});
    var kandM = dil.rename('kandM_'+slug).selfMask();

    var meanB8 = collUnik.select('B8').mean();
    var stdB8  = collUnik.select('B8').reduce(ee.Reducer.stdDev());
    
    // Kandidat W (Dark Pixel)
    var valW = meanB8.multiply(stdB8).rename('valW_'+slug).updateMask(air);
    var thW  = valW.reduceRegion({reducer: ee.Reducer.percentile([5]), geometry: current_roi.geometry(), scale:10, maxPixels:1e9}).get('valW_'+slug);
    var kandW = valW.lte(ee.Number(thW).max(1e-5)).rename('W').selfMask();

    // Kandidat C (Stable Pixel)
    var cv = stdB8.divide(meanB8).rename('cv_'+slug).updateMask(air.eq(0));
    var thC = cv.reduceRegion({reducer: ee.Reducer.percentile([5]), geometry: current_roi.geometry(), scale:10, maxPixels:1e9}).get('cv_'+slug);
    var kandC = cv.lte(ee.Number(thC).max(1e-5)).rename('C').selfMask();

    // Sample M
    var vecM = kandM.sample({region: current_roi.geometry(), scale:10, geometries:true});
    var ids = ee.List.sequence(0, vecM.size().subtract(1));
    var vecM_ID = ee.FeatureCollection(ids.zip(vecM.toList(5000)).map(function(l){ 
      return ee.Feature(ee.Feature(ee.List(l).get(1)).geometry()).set('ID_PIKSEL_M', ee.List(l).get(0)); 
    }));

    // Assign Kelas Debit Sementara
    var collClass = collLengkap.map(function(img){
      var d = ee.Number(img.get('debit'));
      var k = ee.Algorithms.If(d.gte(q2), 'Tinggi', ee.Algorithms.If(d.gte(q1), 'Menengah', 'Rendah'));
      return img.set('kelas_aliran_sementara', k);
    });

    var ekstrakDanGabung = function(image, kandidatC, vektorM_ID){
      var t = image.date().format('YYYY-MM-dd'); var d = image.get('debit'); var k = image.get('kelas_aliran_sementara');
      var nilaiC = ee.Number(image.select('B8').updateMask(kandidatC).reduceRegion({reducer: ee.Reducer.mean(), geometry: current_roi.geometry(), scale: 10, maxPixels: 1e9}).get('B8'));
      var Mplus = image.select('B8').reduceRegions({collection: vektorM_ID, reducer: ee.Reducer.first(), scale: 10});
      return Mplus.map(function(f){
        var nilaiM = ee.Number(f.get('first')); var rasio  = hitungRasio(nilaiC, nilaiM);
        return f.set({'tanggal': t, 'debit': d, 'kelas_aliran': k, 'nilai_C': nilaiC, 'nilai_M': nilaiM, 'rasio_C_per_M': rasio})
                .select(['ID_PIKSEL_M','tanggal','debit','kelas_aliran','nilai_C','nilai_M','rasio_C_per_M'], null, false);
      });
    };

    var dataM = collClass.map(function(img){ return ekstrakDanGabung(img, kandC, vecM_ID); }).flatten();

    // Pilih M Terbaik per Kelas (Spearman)
    var pilihM = function(fcRaw, kls){
      var fc = fcRaw.filter(ee.Filter.eq('kelas_aliran', kls)).filter(ee.Filter.notNull(['ID_PIKSEL_M','debit','rasio_C_per_M']));
      var out = ee.Dictionary(ee.Algorithms.If(fc.size().gt(0), 
        fc.reduceColumns({selectors:['rasio_C_per_M','debit','ID_PIKSEL_M'], reducer: ee.Reducer.spearmansCorrelation().group({groupField:2, groupName:'id'})}), 
        {'groups':[]}
      ));
      var listRho = ee.List(out.get('groups')).map(function(g){
        g=ee.Dictionary(g); var isValid = g.contains('correlation');
        return ee.Algorithms.If(isValid, ee.Feature(null, {'ID': g.get('id'), 'rho': g.get('correlation')}), null);
      }).removeAll([null]);
      var fcRho = ee.FeatureCollection(listRho);
      return ee.Feature(ee.Algorithms.If(fcRho.size().gt(0), fcRho.sort('rho', false).first(), null));
    };

    var mLow = pilihM(dataM, 'Rendah'); var mMed = pilihM(dataM, 'Menengah'); var mHi = pilihM(dataM, 'Tinggi');
    var setProp = function(feat, kls) { return ee.Feature(feat).set({'kelas': kls, 'musim': slug}); };

    var listBest = ee.List([
       ee.Algorithms.If(mLow, setProp(vecM_ID.filter(ee.Filter.eq('ID_PIKSEL_M', mLow.get('ID'))).first(), 'Rendah'), null),
       ee.Algorithms.If(mMed, setProp(vecM_ID.filter(ee.Filter.eq('ID_PIKSEL_M', mMed.get('ID'))).first(), 'Menengah'), null),
       ee.Algorithms.If(mHi,  setProp(vecM_ID.filter(ee.Filter.eq('ID_PIKSEL_M', mHi.get('ID'))).first(), 'Tinggi'), null)
    ]).removeAll([null]);

    return {C: kandC, W: kandW, M: ee.FeatureCollection(listBest)};
  };

  var extractFeaturesStation = function(img, assetSet, slug) {
    var d = img.get('date');
    var mLow = assetSet.M.filter(ee.Filter.eq('kelas','Rendah'));
    var mMed = assetSet.M.filter(ee.Filter.eq('kelas','Menengah'));
    var mHi  = assetSet.M.filter(ee.Filter.eq('kelas','Tinggi'));
    var valW = img.select('B8').updateMask(assetSet.W).reduceRegion(ee.Reducer.mean(), current_roi, 10).get('B8');
    var valC = img.select('B8').updateMask(assetSet.C).reduceRegion(ee.Reducer.mean(), current_roi, 10).get('B8');
    var vMLow = safeReduceFirst(img, 'B8', mLow.geometry(), 10);
    var vMMed = safeReduceFirst(img, 'B8', mMed.geometry(), 10);
    var vMHi  = safeReduceFirst(img, 'B8', mHi.geometry(), 10);
    return ee.Feature(null, {
      'Tanggal': d, 'Musim': slug, 'debit': img.get('debit'),
      'Nilai C': valC, 'Nilai W': valW, 'Nilai Mlow': vMLow, 'Nilai Mmedium': vMMed, 'Nilai Mhigh': vMHi,
      'Rasio C/Mlow': hitungRasio(valC, vMLow), 'Rasio C/Mmedium': hitungRasio(valC, vMMed), 'Rasio C/Mhigh': hitungRasio(valC, vMHi),
      'MWlow': hitungMW(vMLow, valW), 'MWmedium': hitungMW(vMMed, valW), 'MWhigh': hitungMW(vMHi, valW)
    });
  };

  // --- EKSEKUSI PIPELINE ---
  var collLatihReady = kolLatih.map(function(img){ var t = img.get('date'); return ee.Algorithms.If(debitDict.contains(t), img.set('debit', debitDict.get(t)), null); }, true);
  var collUjiReady = kolUji.map(function(img){ var t = img.get('date'); return ee.Algorithms.If(debitDict.contains(t), img.set('debit', debitDict.get(t)), null); }, true);
  var colSel = ['Tanggal','Musim','debit','Nilai C','Nilai W','Nilai Mlow','Nilai Mmedium','Nilai Mhigh','Rasio C/Mlow','Rasio C/Mmedium','Rasio C/Mhigh','MWlow','MWmedium','MWhigh'];

  if (USE_SEASONAL) {
    // Mode Musiman (Jalankan saat USE_SEASONAL = true)
    var collLatihS = collLatihReady.map(function(i){ return i.set('musim', seasonLabel(i.get('date'), current_zom, START_LATIH, END_LATIH)); });
    var collUjiS   = collUjiReady.map(function(i){ return i.set('musim', seasonLabel(i.get('date'), current_zom, START_UJI, END_UJI)); });

    var colH = collLatihS.filter(ee.Filter.eq('musim','Hujan'));
    var colK = collLatihS.filter(ee.Filter.eq('musim','Kemarau'));

    var asetH = defineAssetsForSeason('Hujan', colH, colH, current_thresh.Hujan.q1, current_thresh.Hujan.q2);
    var asetK = defineAssetsForSeason('Kemarau', colK, colK, current_thresh.Kemarau.q1, current_thresh.Kemarau.q2);

    // Export Aset (Opsional)
    var M_H = ASSET_PROJECT_PATH + '/Aset_M_Hujan_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';
    var M_K = ASSET_PROJECT_PATH + '/Aset_M_Kemarau_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';
    var C_H = ASSET_PROJECT_PATH + '/Aset_C_Hujan_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';
    var C_K = ASSET_PROJECT_PATH + '/Aset_C_Kemarau_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';

    Export.table.toAsset({collection: asetH.M, description: 'M_Hujan_'+roi_upper, assetId: M_H});
    Export.table.toAsset({collection: asetK.M, description: 'M_Kemarau_'+roi_upper, assetId: M_K});
    Export.image.toAsset({image: asetH.C, description: 'C_Hujan_'+roi_upper, assetId: C_H, scale: 10, region: current_roi.geometry()});
    Export.image.toAsset({image: asetK.C, description: 'C_Kemarau_'+roi_upper, assetId: C_K, scale: 10, region: current_roi.geometry()});

    var ftLatihH = colH.map(function(i){ return extractFeaturesStation(i, asetH, 'Hujan'); });
    var ftLatihK = colK.map(function(i){ return extractFeaturesStation(i, asetK, 'Kemarau'); });
    var ftUjiH   = collUjiS.filter(ee.Filter.eq('musim','Hujan')).map(function(i){ return extractFeaturesStation(i, asetH, 'Hujan'); });
    var ftUjiK   = collUjiS.filter(ee.Filter.eq('musim','Kemarau')).map(function(i){ return extractFeaturesStation(i, asetK, 'Kemarau'); });

    // Ekspor ke Drive
    Export.table.toDrive({ collection: ftLatihH, description: 'Latih_H_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_latih_'+roi_lower+'_Hujan_'+PERIODE_TAG_Latih+'_individu', selectors: colSel });
    Export.table.toDrive({ collection: ftLatihK, description: 'Latih_K_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_latih_'+roi_lower+'_Kemarau_'+PERIODE_TAG_Latih+'_individu', selectors: colSel });
    Export.table.toDrive({ collection: ftUjiH, description: 'Uji_H_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_uji_'+roi_lower+'_Hujan_'+PERIODE_TAG_Uji+'_individu', selectors: colSel });
    Export.table.toDrive({ collection: ftUjiK, description: 'Uji_K_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_uji_'+roi_lower+'_Kemarau_'+PERIODE_TAG_Uji+'_individu', selectors: colSel });

  } else {
    // Mode Tanpa Musim (Jalankan saat USE_SEASONAL = false)
    var asetG = defineAssetsForSeason('TANPA_MUSIM', collLatihReady, collLatihReady, current_thresh.Gabungan.q1, current_thresh.Gabungan.q2);

    var M_G = ASSET_PROJECT_PATH + '/Aset_M_GABUNGAN_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';
    var C_G = ASSET_PROJECT_PATH + '/Aset_C_GABUNGAN_' + roi_upper + '_' + PERIODE_TAG_Latih + '_Bersih';

    Export.table.toAsset({collection: asetG.M, description: 'M_Gabungan_'+roi_upper, assetId: M_G});
    Export.image.toAsset({image: asetG.C, description: 'C_Gabungan_'+roi_upper, assetId: C_G, scale: 10, region: current_roi.geometry()});

    var ftLatihG = collLatihReady.map(function(i){ return extractFeaturesStation(i, asetG, 'TANPA_MUSIM'); });
    var ftUjiG   = collUjiReady.map(function(i){ return extractFeaturesStation(i, asetG, 'TANPA_MUSIM'); });

    Export.table.toDrive({ collection: ftLatihG, description: 'Latih_G_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_latih_'+roi_lower+'_TANPA_MUSIM_'+PERIODE_TAG_Latih+'_individu', selectors: colSel });
    Export.table.toDrive({ collection: ftUjiG, description: 'Uji_G_'+roi_upper, folder: DRIVE_FOLDER, fileNamePrefix: 'tabel_fitur_uji_'+roi_lower+'_TANPA_MUSIM_'+PERIODE_TAG_Uji+'_individu', selectors: colSel });
  }
}

/* ========== 5. EKSEKUSI UTAMA ========== */

STATION_LIST.forEach(function(station) {
  processStation(station);
});

print('\n🎉 Selesai. Cek Tab TASKS untuk menjalankan ekspor.');