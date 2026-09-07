// ============================================================================
//  CHOROPLETH vs DASYMETRIC — 3 KELAS GFI + ANALISIS SUB-KELOMPOK RENTAN
//  Kabupaten Bone, Sulawesi Selatan, Indonesia  |  v11
// ============================================================================
//
//  IMPORT (klik + → tambah semua aset, pastikan nama variabel persis sama):
//    Banjir_Bone                    → ee.Image  (GFI raster)
//    distribusi_kepadatan_penduduk  → ee.Image  (populasi umum)
//    Jumlah_Penduduk                → ee.FeatureCollection (batas desa)
//    distribusi_umur_rentan              → ee.Image  (Distribusi_Penduduk_Umur_Rentan30)
//    distribusi_miskin              → ee.Image  (Distribusi_Penduduk_Miskin30)
//    distribusi_disabilitas         → ee.Image  (Distribusi_Penduduk_Disabilitas30)
// ============================================================================

var SCALE      = 30;
var MAX_PIX    = 1e11;
var TILE_SCALE = 8;
var NODATA     = -3.4e+38;
var T_LOW      = 0.333;
var T_MED      = 0.666;

// ── STEP 1: AREA STUDI ───────────────────────────────────────────────────────
var boneGeom = Jumlah_Penduduk.geometry().dissolve(1);
Map.centerObject(boneGeom, 9);
print('Jumlah desa:', Jumlah_Penduduk.size());

var desaFC = Jumlah_Penduduk.map(function(f) {
  var area = f.geometry().area(1);
  var pop  = ee.Number(f.get('Populasi'));
  return f.set({ 'area_m2': area, 'bps_pop': pop,
                 'density': pop.divide(area) });
});

// ── STEP 2: RASTER ───────────────────────────────────────────────────────────
var pixArea = ee.Image.pixelArea();
var maskNodata = function(img) {
  return img.select([0]).updateMask(img.select([0]).gt(NODATA + 1e28));
};

var gfi       = maskNodata(Banjir_Bone).clip(boneGeom).rename('gfi');
var popDasy   = maskNodata(distribusi_kepadatan_penduduk).clip(boneGeom).rename('pop');
var popUmurRentan = maskNodata(distribusi_umur_rentan).clip(boneGeom).rename('umur_rentan');
var popMiskin = maskNodata(distribusi_miskin).clip(boneGeom).rename('miskin');
var popDisab  = maskNodata(distribusi_disabilitas).clip(boneGeom).rename('disab');

// Verifikasi total sub-kelompok
print('Total Umur Rentan (jiwa):', popUmurRentan.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: boneGeom,
  scale: SCALE, maxPixels: MAX_PIX }).get('umur_rentan'));
print('Total Miskin (jiwa):', popMiskin.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: boneGeom,
  scale: SCALE, maxPixels: MAX_PIX }).get('miskin'));
print('Total Disabilitas (jiwa):', popDisab.reduceRegion({
  reducer: ee.Reducer.sum(), geometry: boneGeom,
  scale: SCALE, maxPixels: MAX_PIX }).get('disab'));

// ── STEP 3: MASKER 3 KELAS ───────────────────────────────────────────────────
var mkR = gfi.gte(0.000).and(gfi.lt(T_LOW));  // Rendah
var mkS = gfi.gte(T_LOW).and(gfi.lt(T_MED));  // Sedang
var mkT = gfi.gte(T_MED);                      // Tinggi

var gfiClass = ee.Image(0)
  .where(mkR, 1).where(mkS, 2).where(mkT, 3)
  .updateMask(gfi.mask()).rename('class');

// ── STEP 4: CHOROPLETH RASTER ─────────────────────────────────────────────
var choroPop = ee.Image().float()
  .paint(desaFC, 'density').clip(boneGeom)
  .multiply(pixArea).rename('choro');

// ── STEP 5: MULTI-BAND IMAGE — POPULASI UMUM + SUB-KELOMPOK ──────────────
var allBands = popDasy                                   // 'pop'
  .addBands(choroPop)                                    // 'choro'
  // Populasi umum per kelas
  .addBands(popDasy.updateMask(mkR).rename('dR'))
  .addBands(popDasy.updateMask(mkS).rename('dS'))
  .addBands(popDasy.updateMask(mkT).rename('dT'))
  .addBands(choroPop.updateMask(mkR).rename('cR'))
  .addBands(choroPop.updateMask(mkS).rename('cS'))
  .addBands(choroPop.updateMask(mkT).rename('cT'))
  // Sub-kelompok: Umur Rentan (Children + Elderly) — pixel-level spatial overlay
  .addBands(popUmurRentan.rename('umur_rentan_tot'))
  .addBands(popUmurRentan.updateMask(mkR).rename('ur_R'))
  .addBands(popUmurRentan.updateMask(mkS).rename('ur_S'))
  .addBands(popUmurRentan.updateMask(mkT).rename('ur_T'))
  // Sub-kelompok: Miskin (Poor)
  .addBands(popMiskin.rename('miskin_tot'))
  .addBands(popMiskin.updateMask(mkR).rename('mis_R'))
  .addBands(popMiskin.updateMask(mkS).rename('mis_S'))
  .addBands(popMiskin.updateMask(mkT).rename('mis_T'))
  // Sub-kelompok: Disabilitas
  .addBands(popDisab.rename('disab_tot'))
  .addBands(popDisab.updateMask(mkR).rename('dis_R'))
  .addBands(popDisab.updateMask(mkS).rename('dis_S'))
  .addBands(popDisab.updateMask(mkT).rename('dis_T'))
  // Luas per kelas
  .addBands(mkR.unmask(0).multiply(pixArea).rename('aR'))
  .addBands(mkS.unmask(0).multiply(pixArea).rename('aS'))
  .addBands(mkT.unmask(0).multiply(pixArea).rename('aT'))
  .addBands(pixArea.rename('area'));

// ── STEP 6: SATU reduceRegions() ─────────────────────────────────────────
var raw = allBands.reduceRegions({
  collection : desaFC,
  reducer    : ee.Reducer.sum(),
  scale      : SCALE,
  tileScale  : TILE_SCALE
});

// ── STEP 7: HITUNG METRICS PER DESA ──────────────────────────────────────
var exposureFC = raw.map(function(f) {
  var get = function(k) { return ee.Number(f.get(k)); };

  // Populasi umum
  var dTot=get('pop');   var cTot=get('choro');
  var dR=get('dR'); var dS=get('dS'); var dT=get('dT');
  var cR=get('cR'); var cS=get('cS'); var cT=get('cT');
  var aR=get('aR').divide(10000); var aS=get('aS').divide(10000);
  var aT=get('aT').divide(10000); var aA=get('area').divide(10000);

  // Sub-kelompok
  var lTot=get('umur_rentan_tot'); var lR=get('ur_R'); var lS=get('ur_S'); var lT=get('ur_T');
  var mTot=get('miskin_tot'); var mR=get('mis_R'); var mS=get('mis_S'); var mT=get('mis_T');
  var xTot=get('disab_tot');  var xR=get('dis_R'); var xS=get('dis_S'); var xT=get('dis_T');

  // MAUP metrics populasi umum
  var gap    = dT.subtract(cT);
  var gapPct = gap.divide(dT.add(1)).multiply(100);
  var dPct   = dT.divide(dTot.add(1)).multiply(100);
  var cPct   = cT.divide(cTot.add(1)).multiply(100);
  var ssedi  = dPct.divide(cPct.add(0.001));
  var fFrac  = aT.divide(aA.add(0.001)).multiply(100);

  // SSEDI sub-kelompok (dasymetric — pixel-level)
  var lPct    = lT.add(lS).add(lR).divide(lTot.add(1)).multiply(100);
  var mPct    = mT.add(mS).add(mR).divide(mTot.add(1)).multiply(100);
  var xPct    = xT.add(xS).add(xR).divide(xTot.add(1)).multiply(100);
  var refRate = dT.add(dS).add(dR).divide(dTot.add(1)).multiply(100);
  var ssedi_ur = lPct.divide(refRate.add(0.001));
  var ssedi_m = mPct.divide(refRate.add(0.001));
  var ssedi_x = xPct.divide(refRate.add(0.001));

  return f.set({
    'nama_desa'    : f.get('Desa'),
    'bps_pop'      : f.get('Populasi'),
    // Populasi umum
    'dasy_tot'     : dTot.round(), 'choro_tot': cTot.round(),
    'dasy_R'       : dR.round(),   'choro_R'  : cR.round(),
    'dasy_S'       : dS.round(),   'choro_S'  : cS.round(),
    'dasy_T'       : dT.round(),   'choro_T'  : cT.round(),
    'dasy_pct'     : dPct,         'choro_pct': cPct,
    'gap_jiwa'     : gap.round(),  'gap_pct'  : gapPct,
    'ssedi'        : ssedi,
    // Sub-kelompok (pixel-level dasymetric)
    'ur_tot': lTot.round(), 'ur_R': lR.round(), 'ur_S': lS.round(), 'ur_T': lT.round(),
    'mis_tot': mTot.round(), 'mis_R': mR.round(), 'mis_S': mS.round(), 'mis_T': mT.round(),
    'dis_tot': xTot.round(), 'dis_R': xR.round(), 'dis_S': xS.round(), 'dis_T': xT.round(),
    'ssedi_umur_rentan'  : ssedi_ur,
    'ssedi_miskin'  : ssedi_m,
    'ssedi_disab'   : ssedi_x,
    // Luas (ha)
    'ha_R': aR, 'ha_S': aS, 'ha_T': aT, 'desa_ha': aA,
    'flood_frac'   : fFrac,
  });
});

// ── STEP 8: STATISTIK RINGKASAN ──────────────────────────────────────────
var agg = function(p) { return exposureFC.aggregate_sum(p); };

var totPop  = agg('bps_pop');
var totExp  = ee.Number(agg('dasy_R')).add(agg('dasy_S')).add(agg('dasy_T'));
var refRate = totExp.divide(totPop).multiply(100);

print('══════════════════════════════════════════════════');
print('KABUPATEN BONE — Populasi Umum');
print('══════════════════════════════════════════════════');
print('Choro total:', agg('choro_R'), '+', agg('choro_S'), '+', agg('choro_T'));
print('Dasy  total:', agg('dasy_R'),  '+', agg('dasy_S'),  '+', agg('dasy_T'));
print('Reference exposure rate (%):', refRate);

print('\n══════════════════════════════════════════════════');
print('SUB-KELOMPOK RENTAN — Pixel-Level Dasymetric');
print('══════════════════════════════════════════════════');

var printSubgroup = function(label, tot, r, s, t) {
  var totalG  = agg(tot);
  var expAll  = ee.Number(agg(r)).add(agg(s)).add(agg(t));
  var rateG   = expAll.divide(totalG).multiply(100);
  var ssediG  = rateG.divide(refRate);
  print(label);
  print('  Total sub-kelompok:', totalG);
  print('  Terekspos (semua kelas):', expAll, '(', rateG, '%)');
  print('  SSEDI:', ssediG);
};

printSubgroup('▸ Umur Rentan (Children + Elderly):',
  'ur_tot','ur_R','ur_S','ur_T');
printSubgroup('▸ Miskin (Poor):',
  'mis_tot','mis_R','mis_S','mis_T');
printSubgroup('▸ Disabilitas:',
  'dis_tot','dis_R','dis_S','dis_T');

// ── STEP 9: CHARTS ───────────────────────────────────────────────────────
// Chart 1: Choropleth vs Dasymetric — populasi umum
var classFC = ee.FeatureCollection([
  ee.Feature(null,{kelas:'1 Rendah',  choro:agg('choro_R'), dasy:agg('dasy_R')}),
  ee.Feature(null,{kelas:'2 Sedang',  choro:agg('choro_S'), dasy:agg('dasy_S')}),
  ee.Feature(null,{kelas:'3 Tinggi',  choro:agg('choro_T'), dasy:agg('dasy_T')}),
]);
print(ui.Chart.feature.byFeature(classFC,'kelas',['choro','dasy'])
  .setSeriesNames(['Choropleth','Dasymetric'])
  .setChartType('ColumnChart')
  .setOptions({
    title:'Flood Exposure: Choropleth vs Dasymetric — 3 GFI Classes',
    colors:['#E24B4A','#2E75B6'],
    bar:{groupWidth:'65%'}, legend:{position:'top'},
    vAxis:{title:'Population exposed (persons)',format:'#,###'},
    hAxis:{title:'GFI Hazard Class'},
    backgroundColor:'#FAFAFA',
  }));

// Chart 2: SSEDI per sub-kelompok
ee.List([
  agg('ur_tot'), ee.Number(agg('ur_R')).add(agg('ur_S')).add(agg('ur_T')),
  agg('mis_tot'), ee.Number(agg('mis_R')).add(agg('mis_S')).add(agg('mis_T')),
  agg('dis_tot'), ee.Number(agg('dis_R')).add(agg('dis_S')).add(agg('dis_T')),
  totPop, totExp, refRate
]).evaluate(function(v) {
  if (!v) return;
  var ssediL = (v[1]/v[0]*100) / v[8];
  var ssediM = (v[3]/v[2]*100) / v[8];
  var ssediX = (v[5]/v[4]*100) / v[8];

  var ssediFC = ee.FeatureCollection([
    ee.Feature(null,{group:'Umur Rentan',      ssedi:ssediL, exposed:v[1], total:v[0]}),
    ee.Feature(null,{group:'Poor',         ssedi:ssediM, exposed:v[3], total:v[2]}),
    ee.Feature(null,{group:'Disability',   ssedi:ssediX, exposed:v[5], total:v[4]}),
  ]);
  print(ui.Chart.feature.byFeature(ssediFC,'group',['ssedi'])
    .setSeriesNames(['SSEDI (dasymetric)'])
    .setChartType('ColumnChart')
    .setOptions({
      title:'Social-Spatial Exposure Disparity Index (SSEDI) by Sub-Group',
      colors:['#2E75B6'], legend:{position:'none'},
      vAxis:{title:'SSEDI (1.0 = parity with general population)'},
      hAxis:{title:'Vulnerable Sub-Group'},
      backgroundColor:'#FAFAFA',
    }));

  statLbl.setValue(
    'General pop exposure rate : ' + v[8].toFixed(2) + '%\n' +
    'SSEDI Elderly   : ' + ssediL.toFixed(3) + '\n' +
    'SSEDI Poor      : ' + ssediM.toFixed(3) + '\n' +
    'SSEDI Disability: ' + ssediX.toFixed(3)
  );
});

// ── STEP 10: SPLIT PANEL MAP ──────────────────────────────────────────────
var palPop = ['#FFF7EC','#FDD49E','#FC8D59','#D7301F','#7F0000'];
var palGFI = ['#64B5F6','#FFB74D','#E53935'];
var desaBdr = ee.Image().paint(Jumlah_Penduduk,0,0.7);

var leftMap = ui.Map();
leftMap.setOptions('TERRAIN'); leftMap.setCenter(120.35,-4.35,9);
leftMap.addLayer(choroPop.visualize({min:0,max:30,palette:palPop}),
  {opacity:0.85},'Choropleth (persons/pixel)');
leftMap.addLayer(gfiClass.visualize({min:1,max:3,palette:palGFI}),
  {opacity:0.50},'GFI 3 Classes');
leftMap.addLayer(desaBdr,{palette:['#FFFFFF'],opacity:0.55},'Desa Boundary');

var rightMap = ui.Map();
rightMap.setOptions('TERRAIN');
rightMap.addLayer(popDasy.visualize({min:0,max:30,palette:palPop}),
  {opacity:0.85},'Dasymetric (persons/pixel)');
rightMap.addLayer(gfiClass.visualize({min:1,max:3,palette:palGFI}),
  {opacity:0.50},'GFI 3 Classes');
rightMap.addLayer(desaBdr,{palette:['#FFFFFF'],opacity:0.55},'Desa Boundary');

var linker   = ui.Map.Linker([leftMap, rightMap]);
var splitMap = ui.SplitPanel({
  firstPanel:leftMap, secondPanel:rightMap,
  orientation:'horizontal', wipe:true
});

// Info panel
var info = ui.Panel({style:{position:'bottom-left',padding:'10px',
  width:'270px',backgroundColor:'rgba(255,255,255,0.93)',border:'1px solid #CCC'}});
info.add(ui.Label('◀ Choropleth  |  Dasymetric ▶',
  {fontWeight:'bold',fontSize:'13px',color:'#1F3864'}));
info.add(ui.Label('Kabupaten Bone · 3 GFI Classes + Vulnerable Sub-Groups',
  {fontSize:'10px',color:'#888',margin:'2px 0 8px'}));
var statLbl = ui.Label('Computing...',{fontSize:'10px',color:'#444'});
info.add(statLbl);
info.add(ui.Label('GFI Legend:',
  {fontWeight:'bold',fontSize:'10px',color:'#1F3864',margin:'8px 0 3px'}));
var leg = function(col,lbl) {
  return ui.Panel([
    ui.Label('',{backgroundColor:col,padding:'5px 9px',margin:'1px 4px 1px 0'}),
    ui.Label(lbl,{fontSize:'10px',color:'#333',margin:'1px 0'})
  ],ui.Panel.Layout.flow('horizontal'));
};
[['#64B5F6','Low   0.000–0.333'],
 ['#FFB74D','Moderate 0.333–0.666'],
 ['#E53935','High  0.666–1.000']
].forEach(function(r){info.add(leg(r[0],r[1]));});
leftMap.add(info);
ui.root.clear();
ui.root.add(splitMap);

// ── STEP 11: EXPORT ──────────────────────────────────────────────────────
Export.table.toDrive({
  collection: exposureFC.select([
    'nama_desa','bps_pop',
    'dasy_tot','dasy_R','dasy_S','dasy_T','dasy_pct',
    'choro_tot','choro_R','choro_S','choro_T','choro_pct',
    'gap_jiwa','gap_pct','ssedi',
    'ur_tot','ur_R','ur_S','ur_T',
    'mis_tot','mis_R','mis_S','mis_T',
    'dis_tot','dis_R','dis_S','dis_T',
    'ssedi_umur_rentan','ssedi_miskin','ssedi_disab',
    'ha_R','ha_S','ha_T','desa_ha','flood_frac'
  ]),
  description   :'bone_exposure_vulnerable_v11',
  folder        :'GEE_Exports',
  fileNamePrefix:'bone_choropleth_dasymetric_vulnerable',
  fileFormat    :'CSV',
});

Export.image.toDrive({
  image:gfiClass, description:'bone_gfi_3class',
  folder:'GEE_Exports', fileNamePrefix:'bone_gfi_classified',
  region:boneGeom, scale:SCALE, maxPixels:MAX_PIX,
  fileFormat:'GeoTIFF', crs:'EPSG:32751'
});

print('✓ Done. Tasks → Run to export.');
