const pptxgen = require('pptxgenjs');
const fs = require('fs');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'Balázs Kakuk';
pptx.subject = 'EHV-1 dynamic transcriptome conference presentation';
pptx.title = 'Long-read transcriptomics reveals the dynamic EHV-1 transcriptional program';
pptx.company = 'University of Szeged';
pptx.lang = 'en-US';
pptx.theme = {
  headFontFace: 'Aptos Display',
  bodyFontFace: 'Aptos',
  lang: 'en-US'
};
pptx.defineSlideMaster({
  title: 'MASTER',
  background: { color: 'F7F9FC' },
  objects: [
    { rect: { x: 0, y: 0, w: 13.333, h: 0.12, fill: { color: '1F8A8A' }, line: { color: '1F8A8A' } } },
    { text: { text: 'EHV-1 dynamic transcriptome', options: { x: 0.55, y: 7.12, w: 4.8, h: 0.18, fontFace: 'Aptos', fontSize: 8.5, color: '8A97A6', margin: 0 } } },
    { text: { text: 'Tombácz et al. · PLOS ONE 2025 · doi:10.1371/journal.pone.0320439', options: { x: 7.15, y: 7.12, w: 5.6, h: 0.18, fontFace: 'Aptos', fontSize: 8.5, color: '8A97A6', align: 'right', margin: 0 } } },
  ],
  slideNumber: { x: 12.82, y: 7.11, w: 0.25, h: 0.18, fontFace: 'Aptos', fontSize: 8, color: '8A97A6', align: 'right' }
});

const C = {
  navy: '12263A',
  teal: '1F8A8A',
  cyan: '56BFC2',
  blue: '4E79A7',
  orange: 'F28E2B',
  red: 'E15759',
  green: '59A14F',
  purple: '9C6FB6',
  gray: '667085',
  mid: 'A7B1BE',
  light: 'E8EEF3',
  paleTeal: 'E8F6F5',
  paleBlue: 'EAF2F8',
  paleOrange: 'FFF2E3',
  paleRed: 'FCEBEC',
  white: 'FFFFFF',
  bg: 'F7F9FC',
};
const shadow = { type: 'outer', color: '000000', opacity: 0.13, blur: 1.5, angle: 45, distance: 0.7 };

function addTitle(slide, title, kicker) {
  if (kicker) slide.addText(kicker.toUpperCase(), { x: 0.65, y: 0.38, w: 4.5, h: 0.25, fontSize: 10.5, bold: true, color: C.teal, charSpacing: 1.6, margin: 0 });
  slide.addText(title, { x: 0.65, y: kicker ? 0.66 : 0.48, w: 12.0, h: 0.54, fontSize: 28, bold: true, color: C.navy, margin: 0, breakLine: false });
}
function addSource(slide, text) {
  slide.addText(text, { x: 0.7, y: 6.82, w: 11.9, h: 0.20, fontSize: 8.5, italic: true, color: '7B8794', margin: 0, align: 'left' });
}
function addNote(slide, text) { slide.addNotes(text); }
function addCallout(slide, x, y, w, h, big, label, color=C.teal, fill=C.white) {
  slide.addText([
    { text: big + '\n', options: { bold: true, fontSize: 25, color } },
    { text: label, options: { fontSize: 12.5, color: C.gray } }
  ], { x, y, w, h, margin: 0.16, valign: 'mid', align: 'center', fill: { color: fill }, line: { color: C.light, width: 1 }, radius: 0.08, shadow });
}
function addPill(slide, text, x, y, w, color, fill) {
  slide.addText(text, { x, y, w, h: 0.34, fontSize: 11.2, bold: true, color, align: 'center', valign: 'mid', margin: 0.04, fill: { color: fill }, line: { color: fill }, radius: 0.17 });
}
function addSectionLabel(slide, text, x, y, w=2.1) {
  slide.addText(text, { x, y, w, h: 0.31, fontSize: 10.5, bold: true, color: C.teal, margin: 0, charSpacing: 1.0 });
}
function lineChart(slide, x, y, w, h, series, categories, opts={}) {
  slide.addChart(pptx.ChartType.line, series.map(s => ({ name: s.name, labels: categories, values: s.values })), {
    x, y, w, h,
    showTitle: false,
    showLegend: opts.showLegend !== false,
    legendPos: opts.legendPos || 'b',
    legendFontSize: 10,
    catAxisLabelFontSize: 10,
    valAxisLabelFontSize: 9,
    catAxisTitle: opts.catTitle || 'hours post-infection',
    valAxisTitle: opts.valTitle || 'relative abundance',
    showCatName: false,
    showValAxisTitle: true,
    showCatAxisTitle: true,
    valAxisMinVal: opts.minVal === undefined ? 0 : opts.minVal,
    valAxisMaxVal: opts.maxVal,
    valAxisMajorUnit: opts.majorUnit,
    chartColors: opts.colors || [C.blue, C.orange, C.green, C.red, C.purple],
    showValue: false,
    lineSize: 2.3,
    markerSize: 5,
    showBorder: false,
    showGridLines: true,
    gridLine: { color: 'D9E0E7', width: 1 },
    catAxisLineColor: 'AAB5C1',
    valAxisLineColor: 'AAB5C1',
    showCatAxisLabel: true,
    showValAxisLabel: true,
  });
}
function simpleBar(slide, x, y, w, h, cats, vals, colors) {
  slide.addChart(pptx.ChartType.bar, [{ name: 'relative expression', labels: cats, values: vals }], {
    x, y, w, h,
    showLegend: false, showTitle: false,
    showValue: true, dataLabelPosition: 'outEnd', dataLabelColor: C.navy, dataLabelFormatCode: '0.0x',
    catAxisLabelFontSize: 13, valAxisLabelFontSize: 9,
    chartColors: colors,
    valAxisMinVal: 0, valAxisMaxVal: 13, valAxisMajorUnit: 2,
    showGridLines: true, gridLine: { color: 'E1E6EC' },
    showBorder: false,
  });
}
function transcript(slide, x, y, w, label, color, startFrac=0, endFrac=1, dashed=false) {
  const sx = x + startFrac*w, ex = x + endFrac*w;
  slide.addShape(pptx.ShapeType.line, { x: sx, y: y+0.12, w: ex-sx, h: 0, line: { color, width: 2.4, dash: dashed ? 'dash' : 'solid', beginArrowType: 'none', endArrowType: 'triangle' } });
  slide.addText(label, { x: sx, y: y-0.04, w: Math.max(0.7, ex-sx), h: 0.18, fontSize: 9.5, color, margin: 0 });
}
function geneArrow(slide, x, y, w, label, color, reverse=false) {
  slide.addShape(pptx.ShapeType.chevron, { x, y, w, h: 0.36, rotate: reverse ? 180 : 0, fill: { color }, line: { color }, radius: 0.04 });
  slide.addText(label, { x, y: y+0.05, w, h: 0.2, fontSize: 9.5, bold: true, align: 'center', color: C.white, margin: 0 });
}

// 1. Title slide
{
  const s = pptx.addSlide();
  s.background = { color: C.navy };
  s.addShape(pptx.ShapeType.rect, { x: 0, y: 0, w: 13.333, h: 0.16, fill: { color: C.teal }, line: { color: C.teal } });
  s.addText('LONG-READ VIRAL TRANSCRIPTOMICS', { x: 0.78, y: 0.70, w: 5.1, h: 0.28, fontSize: 11.5, bold: true, color: '75D1CE', charSpacing: 1.8, margin: 0 });
  s.addText('Long-read transcriptomics reveals\nthe dynamic EHV-1 transcriptional program', { x: 0.78, y: 1.12, w: 8.8, h: 1.42, fontSize: 33, bold: true, color: C.white, margin: 0, breakLine: false });
  s.addText('Integrating nanopore direct cDNA sequencing, CAGE-Seq and native RNA validation across infection', { x: 0.82, y: 2.77, w: 8.7, h: 0.64, fontSize: 18, color: 'D6E2EA', margin: 0 });
  // abstract transcript motif
  for (let i=0;i<6;i++) {
    const yy=4.05+i*0.38;
    const start=7.8+(i%2)*0.3;
    s.addShape(pptx.ShapeType.line, { x: start, y: yy, w: 4.25-(i*0.24), h: 0, line: { color: i%2? '56BFC2':'4E79A7', width: 2.1, transparency: 10, endArrowType: 'triangle' } });
    s.addShape(pptx.ShapeType.ellipse, { x: start-0.06, y: yy-0.06, w: 0.12, h: 0.12, fill: { color: 'F28E2B' }, line: { color: 'F28E2B' } });
  }
  s.addText('Balázs Kakuk · Dóra Tombácz · Gábor Torma · Ádám Fülöp · Ákos Dörmő · Gábor Gulyás · Zsolt Csabai · Zsolt Boldogkői', { x: 0.82, y: 5.40, w: 10.8, h: 0.55, fontSize: 13.2, color: 'D6E2EA', margin: 0 });
  s.addText('Department of Medical Biology · University of Szeged', { x: 0.82, y: 6.03, w: 7.5, h: 0.28, fontSize: 12.5, color: 'AFC0CC', margin: 0 });
  s.addText('Based on Tombácz et al., PLOS ONE 20(4): e0320439 (2025)', { x: 0.82, y: 6.47, w: 6.9, h: 0.24, fontSize: 10.5, color: '8FA8B7', margin: 0 });
  addNote(s, 'Opening: EHV-1 is a compact DNA virus, but its transcriptome is anything but simple. In this study we combined full-length nanopore sequencing with cap-specific CAGE data to follow transcript architecture over the entire infection cycle. The key point of the talk is that gene expression is not only changing in amount — the identity and structure of the RNA molecules changes as well.');
}

// 2. Why hard
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Why is a 150 kb viral genome transcriptomically difficult?', 'Problem');
  s.addText('EHV-1 packages dense transcriptional logic into a small genome.', { x: 0.7, y: 1.30, w: 6.2, h: 0.42, fontSize: 19, color: C.gray, margin: 0 });
  // genome strip
  const gx=0.85, gy=2.02, gw=7.05;
  s.addText('~150 kb dsDNA genome · 76 unique protein-coding genes', { x: gx, y: gy-0.38, w: gw, h: 0.25, fontSize: 12.5, bold: true, color: C.navy, margin: 0 });
  const widths=[0.75,0.55,0.9,0.6,1.0,0.7,0.9,0.7,0.65,0.9];
  let xx=gx;
  widths.forEach((ww,i)=>{ geneArrow(s,xx,gy,ww,'ORF'+(i+1), i%2?C.blue:C.teal, i===6); xx+=ww-0.06; });
  s.addText('Genes are tandem, convergent and divergent — and many transcripts share sequence.', { x: gx, y: gy+0.58, w: gw, h: 0.42, fontSize: 14.5, color: C.gray, margin: 0 });
  // transcript architecture examples
  s.addText('One genomic region can generate…', { x: 0.85, y: 3.30, w: 3.0, h: 0.3, fontSize: 14, bold: true, color: C.navy, margin: 0 });
  transcript(s,0.95,3.78,6.5,'canonical',C.teal,0.03,0.55);
  transcript(s,0.95,4.18,6.5,'long 5′-UTR isoform',C.blue,0.03,0.78);
  transcript(s,0.95,4.58,6.5,'multicistronic RNA',C.orange,0.03,0.96);
  transcript(s,0.95,4.98,6.5,'5′-truncated / embedded mRNA',C.red,0.30,0.55);
  // right column classical model
  s.addText('Classical kinetic model', { x: 8.45, y: 1.38, w: 3.7, h: 0.3, fontSize: 15.5, bold: true, color: C.navy, margin: 0 });
  const bx=8.50, bw=3.65;
  s.addText('IE', { x: bx, y: 1.98, w: 0.78, h: 0.78, fontSize: 25, bold: true, color: C.white, align: 'center', valign: 'mid', margin: 0, fill: { color: C.red }, radius: 0.10 });
  s.addText('E', { x: bx+1.05, y: 2.26, w: 0.78, h: 0.78, fontSize: 25, bold: true, color: C.white, align: 'center', valign: 'mid', margin: 0, fill: { color: C.orange }, radius: 0.10 });
  s.addText('L', { x: bx+2.10, y: 2.54, w: 0.78, h: 0.78, fontSize: 25, bold: true, color: C.white, align: 'center', valign: 'mid', margin: 0, fill: { color: C.green }, radius: 0.10 });
  s.addShape(pptx.ShapeType.line, { x: bx+0.82, y: 2.36, w: 0.22, h: 0.18, line: { color: C.mid, width: 1.5, endArrowType: 'triangle' } });
  s.addShape(pptx.ShapeType.line, { x: bx+1.87, y: 2.64, w: 0.22, h: 0.18, line: { color: C.mid, width: 1.5, endArrowType: 'triangle' } });
  s.addText('Useful framework — but transcript boundaries, splicing and overlapping RNAs can blur it.', { x: 8.45, y: 3.65, w: 3.8, h: 1.0, fontSize: 16.5, color: C.gray, margin: 0.03, breakLine: false });
  s.addText('Question', { x: 8.45, y: 4.92, w: 1.1, h: 0.27, fontSize: 11, bold: true, color: C.teal, margin: 0, charSpacing: 1.2 });
  s.addText('Can full-length, time-resolved RNA sequencing reveal the real regulatory program?', { x: 8.45, y: 5.27, w: 3.85, h: 0.92, fontSize: 20, bold: true, color: C.navy, margin: 0 });
  addSource(s, 'Source: published PLOS ONE article; EHV-1 genome architecture and classical IE/E/L framework.');
  addNote(s, 'Set up the problem. EHV-1 has a roughly 150 kb double-stranded DNA genome and 76 unique protein-coding genes, yet adjacent genes frequently share transcript ends or are embedded in longer RNAs. The classical immediate-early, early and late scheme describes gene-level timing, but it does not necessarily describe which full-length RNA isoform is present.');
}

// 3. Complementary sequencing
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'No single sequencing method resolves the whole transcript', 'Strategy');
  s.addText('The study combines complementary strengths rather than treating one platform as ground truth.', { x: 0.7, y: 1.26, w: 9.5, h: 0.38, fontSize: 18.5, color: C.gray, margin: 0 });
  const cols = [
    {x:0.75, title:'CAGE-Seq', sub:'Illumina MiSeq', color:C.orange, fill:C.paleOrange, big:'5′ cap', bullets:['High-resolution capped TSS signal','Strong evidence for promoter usage','Does not link 5′ end to a full transcript']},
    {x:4.52, title:'Direct cDNA-Seq', sub:'ONT MinION', color:C.teal, fill:C.paleTeal, big:'full length', bullets:['Long reads link TSS → splice → TES','Adapter evidence gives orientation','Enables isoform-specific kinetics']},
    {x:8.29, title:'Direct RNA-Seq', sub:'previous dataset', color:C.blue, fill:C.paleBlue, big:'native RNA', bullets:['Independent validation of introns/TSSs','No RT/PCR artifacts','5′ truncation limits TSS precision']},
  ];
  cols.forEach(c=>{
    s.addText(c.title, { x:c.x, y:1.90, w:3.28, h:0.36, fontSize:21, bold:true, color:C.navy, margin:0 });
    s.addText(c.sub, { x:c.x, y:2.28, w:3.28, h:0.25, fontSize:11.5, color:C.gray, margin:0 });
    s.addText(c.big, { x:c.x, y:2.72, w:3.28, h:0.72, fontSize:27, bold:true, color:c.color, align:'center', valign:'mid', margin:0, fill:{color:c.fill}, line:{color:c.fill}, radius:0.12 });
    c.bullets.forEach((b,i)=>s.addText('• '+b,{x:c.x+0.05,y:3.73+i*0.62,w:3.15,h:0.48,fontSize:14.2,color:C.gray,margin:0.02,breakLine:false}));
  });
  s.addShape(pptx.ShapeType.line, { x: 3.92, y: 3.05, w: 0.42, h: 0, line: { color: C.mid, width: 2, endArrowType:'triangle' } });
  s.addShape(pptx.ShapeType.line, { x: 7.69, y: 3.05, w: 0.42, h: 0, line: { color: C.mid, width: 2, endArrowType:'triangle' } });
  s.addText('Integrated evidence', { x: 4.72, y: 6.06, w: 3.0, h: 0.38, fontSize: 17.5, bold:true, color:C.white, align:'center', valign:'mid', margin:0, fill:{color:C.navy}, radius:0.16 });
  addSource(s, 'CAGE-Seq + ONT dcDNA-Seq in this study; earlier ONT dRNA-Seq reanalyzed for validation.');
  addNote(s, 'The central methodological idea is complementarity. CAGE is excellent for capped 5-prime ends. Direct cDNA provides long reads that connect those starts to splice junctions and transcript ends. The earlier native direct-RNA dataset is especially useful as an independent validation source, although its 5-prime ends are less reliable because of truncation.');
}

// 4. Experimental design
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'A dense time course captures the complete infection cycle', 'Experimental design');
  s.addText('RK-13 cells infected with the EHV-1 MdBio field isolate', { x:0.7,y:1.28,w:7.4,h:0.38,fontSize:18.5,color:C.gray,margin:0 });
  // timeline
  const times=[1,2,4,6,8,12,18,24,48];
  const x0=0.95, x1=12.0, yy=3.05;
  s.addShape(pptx.ShapeType.line,{x:x0,y:yy,w:x1-x0,h:0,line:{color:C.navy,width:2.2,endArrowType:'triangle'}});
  times.forEach((t,i)=>{
    const frac=i/(times.length-1); const x=x0+frac*(x1-x0-0.15);
    s.addShape(pptx.ShapeType.ellipse,{x:x-0.10,y:yy-0.10,w:0.20,h:0.20,fill:{color:i<3?C.orange:(i<6?C.teal:C.blue)},line:{color:C.white,width:1}});
    s.addText(String(t),{x:x-0.24,y:yy+0.18,w:0.48,h:0.25,fontSize:12,bold:true,color:C.navy,align:'center',margin:0});
  });
  s.addText('hours post-infection', {x:5.15,y:3.54,w:2.7,h:0.28,fontSize:11.5,color:C.gray,align:'center',margin:0});
  addCallout(s,0.85,4.35,2.5,1.18,'9','time points',C.teal,C.white);
  addCallout(s,3.67,4.35,2.5,1.18,'×3','replicates / time point',C.blue,C.white);
  addCallout(s,6.49,4.35,2.5,1.18,'27','dcDNA libraries',C.orange,C.white);
  addCallout(s,9.31,4.35,2.5,1.18,'3.33 M','viral long reads',C.red,C.white);
  // CHX mini box
  s.addText('Immediate-early test', {x:0.85,y:5.92,w:2.1,h:0.30,fontSize:12,bold:true,color:C.teal,margin:0});
  s.addText('Cycloheximide: 20 or 100 µg/mL · samples at 6 and 8 hpi', {x:2.60,y:5.88,w:6.2,h:0.40,fontSize:14.5,color:C.gray,margin:0});
  s.addText('Protein synthesis blocked → transcripts that remain strongly expressed can be classified as immediate-early.', {x:2.60,y:6.28,w:8.4,h:0.40,fontSize:13.3,color:C.gray,margin:0});
  addSource(s, 'Published design: 1–48 hpi time course; 3 replicates/time point; CHX treatment at 20 and 100 µg/mL.');
  addNote(s, 'The time course is one of the strengths of the study: nine points from one hour to 48 hours post infection, with three replicates at each point, yielding 27 direct-cDNA libraries and about 3.33 million viral long reads. A separate cycloheximide experiment blocks new protein synthesis and tests which viral genes can still be expressed as true immediate-early genes.');
}

// 5. Workflow
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'From nanopore reads to validated viral transcripts', 'Integrated workflow');
  const y=1.75;
  const steps=[
    ['ONT dcDNA','Dorado → minimap2',C.teal,C.paleTeal],
    ['LoRTIA','adapter-aware transfrags',C.blue,C.paleBlue],
    ['CAGEfightR','capped TSS clusters',C.orange,C.paleOrange],
    ['TSS refinement','dcDNA 5′ peaks inside CAGE',C.red,C.paleRed],
    ['TSS ↔ TES pairing','full-length transcript models',C.purple,'F4ECF8'],
    ['dRNA / NAGATA','independent validation',C.green,'EEF7EA'],
  ];
  const sx=0.72, gap=0.18, sw=1.94;
  steps.forEach((st,i)=>{
    const x=sx+i*(sw+gap);
    s.addText(st[0],{x,y,w:sw,h:0.52,fontSize:15.5,bold:true,color:C.navy,align:'center',valign:'mid',margin:0.06,fill:{color:st[3]},line:{color:st[2],width:1.4},radius:0.08});
    s.addText(st[1],{x,y:y+0.62,w:sw,h:0.62,fontSize:11.2,color:C.gray,align:'center',valign:'top',margin:0.04});
    if(i<steps.length-1) s.addShape(pptx.ShapeType.line,{x:x+sw+0.02,y:y+0.26,w:gap-0.04,h:0,line:{color:C.mid,width:1.6,endArrowType:'triangle'}});
  });
  s.addText('Key integration rule', {x:0.80,y:3.42,w:2.1,h:0.3,fontSize:12,bold:true,color:C.teal,margin:0});
  s.addText('A novel transcript had to connect a refined CAGE-supported TSS to a validated TES on dcDNA reads with correct 5′ adapter evidence.', {x:0.80,y:3.78,w:11.4,h:0.72,fontSize:19.5,bold:true,color:C.navy,margin:0.02,breakLine:false});
  const criteria=[
    ['≥3 dcDNA reads','same TSS + TES'],
    ['correct 5′ adapter','orientation / intact end'],
    ['CAGE-supported TSS','high-resolution 5′ evidence'],
    ['validated TES','±10 nt matching'],
  ];
  criteria.forEach((c,i)=>{
    const x=0.88+i*3.02;
    s.addText([{text:c[0]+'\n',options:{bold:true,fontSize:17,color:C.teal}},{text:c[1],options:{fontSize:11.5,color:C.gray}}],{x,y:5.00,w:2.64,h:0.88,fill:{color:C.white},line:{color:C.light},margin:0.12,align:'center',valign:'mid',radius:0.08,shadow});
  });
  s.addText('For 5′-truncated ORF-carrying isoforms, additional dRNA support and ≥5% CAGE signal relative to the canonical TSS were required.', {x:0.90,y:6.18,w:11.5,h:0.42,fontSize:12.8,color:C.gray,align:'center',margin:0});
  addSource(s, 'Final revised workflow: LoRTIA + CAGEfightR + custom TSS refinement + NAGATA validation.');
  addNote(s, 'This is the methodological core. Long-read alignments were converted into nucleotide-resolved transfrags retaining adapter evidence. CAGEfightR identifies TSS clusters, and broad clusters were refined using the actual dcDNA 5-prime-end peaks. Refined TSSs are then paired to validated TESs. Novel transcripts require multiple reads, correct adapters and agreement with the independent CAGE evidence. Putative embedded mRNAs received an extra stringent filter.');
}

// 6. Reannotation outcome
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Integrated evidence substantially expands the EHV-1 transcriptome', 'Reannotation');
  addCallout(s,0.80,1.55,3.25,1.28,'169','novel transcripts',C.teal,C.paleTeal);
  addCallout(s,4.28,1.55,3.25,1.28,'3.33 M','viral dcDNA reads',C.blue,C.paleBlue);
  addCallout(s,7.76,1.55,4.10,1.28,'220 / 34 / 36','high / medium / low CAGE support\nfor previously annotated TSSs',C.orange,C.paleOrange);
  s.addText('Selected novel transcript categories reported in the paper', {x:0.82,y:3.22,w:6.2,h:0.34,fontSize:15.5,bold:true,color:C.navy,margin:0});
  const cats=[
    ['80','monocistronic\nUTR isoforms',C.teal,C.paleTeal],
    ['18','multicistronic\nUTR isoforms',C.blue,C.paleBlue],
    ['26','non-coding\nRNAs',C.purple,'F4ECF8'],
    ['29','putative embedded\nmRNAs',C.red,C.paleRed],
  ];
  cats.forEach((c,i)=>{
    const x=0.85+i*2.92;
    s.addText([{text:c[0]+'\n',options:{bold:true,fontSize:26,color:c[2]}},{text:c[1],options:{fontSize:12.5,color:C.gray}}],{x,y:3.75,w:2.55,h:1.25,fill:{color:c[3]},line:{color:c[3]},margin:0.12,align:'center',valign:'mid',radius:0.10});
  });
  s.addText('These categories are highlighted examples, not an exhaustive partition of all 169 new transcript models.', {x:0.90,y:5.22,w:11.4,h:0.28,fontSize:10.7,italic:true,color:'7B8794',align:'center',margin:0});
  // schematic before-after
  s.addText('Previous annotation', {x:1.00,y:5.75,w:2.0,h:0.25,fontSize:11.5,bold:true,color:C.gray,margin:0});
  transcript(s,0.95,6.08,4.3,'canonical',C.mid,0.02,0.62);
  s.addShape(pptx.ShapeType.line,{x:5.15,y:6.20,w:1.2,h:0,line:{color:C.teal,width:3,endArrowType:'triangle'}});
  s.addText('Refined annotation', {x:6.70,y:5.75,w:2.1,h:0.25,fontSize:11.5,bold:true,color:C.gray,margin:0});
  transcript(s,6.65,5.98,5.3,'canonical',C.teal,0.02,0.55);
  transcript(s,6.65,6.25,5.3,'long / polycistronic',C.blue,0.02,0.90);
  transcript(s,6.65,6.52,5.3,'short / embedded',C.red,0.28,0.55);
  addSource(s, 'Published article: 169 novel transcripts; final published CAGE-support counts 220 high, 34 medium, 36 low.');
  addNote(s, 'The integration has a concrete annotation payoff: 169 novel transcript models. The paper explicitly reports 80 monocistronic UTR isoforms, 18 multicistronic UTR isoforms, 26 non-coding RNAs and 29 putative embedded mRNAs. I would not add these numbers as if they exhaust the 169, because the published category list is not a complete partition. The important message is that the EHV-1 transcriptome is substantially richer than a one-gene-one-transcript view.');
}

// 7. ORF64 IE
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Protein-synthesis inhibition isolates a single immediate-early gene', 'Immediate-early program');
  s.addText('Cycloheximide blocks production of new viral regulatory proteins.', {x:0.75,y:1.30,w:7.7,h:0.38,fontSize:18.2,color:C.gray,margin:0});
  simpleBar(s,0.78,1.98,6.0,3.65,['ORF67\n(next highest)','ORF64'],[1,12.1],[C.mid,C.red]);
  s.addText('ORF64', {x:7.35,y:2.05,w:2.8,h:0.55,fontSize:32,bold:true,color:C.red,margin:0});
  s.addText('is the sole EHV-1 immediate-early gene', {x:7.35,y:2.67,w:4.4,h:0.58,fontSize:21,bold:true,color:C.navy,margin:0});
  const stats=[['12.1×','higher than ORF67'],['Z = 2.84','SD above mean'],['p < 10⁻⁴⁷','one-sample t-test']];
  stats.forEach((c,i)=>s.addText([{text:c[0]+'\n',options:{bold:true,fontSize:21,color:i===0?C.red:C.teal}},{text:c[1],options:{fontSize:11.2,color:C.gray}}],{x:7.35+i*1.72,y:3.65,w:1.52,h:1.05,fill:{color:C.white},line:{color:C.light},margin:0.10,align:'center',valign:'mid',radius:0.08,shadow}));
  s.addText('This anchors the start of the cascade — but the later E/L program is much less discrete.', {x:7.35,y:5.24,w:4.7,h:0.72,fontSize:18.5,color:C.gray,margin:0});
  addSource(s, 'Published CHX analysis: ORF64 is a strong outlier and the only gene with significant immediate-early expression.');
  addNote(s, 'Cycloheximide gives a clean biological anchor for the kinetic program. ORF64 remains dramatically more abundant than every other viral gene when new protein synthesis is blocked. Its expression is 12.1-fold higher than the second most abundant gene, ORF67, with an extreme statistical signal. This confirms ORF64 as the sole immediate-early gene.');
}

// 8. Classical waves
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Canonical full-length transcripts recover the expected early and late waves', 'Temporal expression');
  s.addText('Representative canonical transcripts, normalized to each gene’s own maximum to compare timing', {x:0.73,y:1.28,w:9.8,h:0.32,fontSize:14.5,color:C.gray,margin:0});
  const times=['1','2','4','6','8','12','18','24','48'];
  function norm(a){const m=Math.max(...a); return a.map(v=>m===0?0:v/m);}
  const ORF5=[0,0.0920865791,0.0684476086,0.0103709175,0.0077843498,0.0041208329,0.0005589175,0.0005045575,0.0004376512];
  const ORF31=[0,0.0668797261,0.0729149674,0.0068035658,0.0205094964,0.0137840846,0.0011804488,0.0037964139,0.00221006];
  const ORF22=[0,0.0026442349,0.0052025087,0.0034846206,0.0090849997,0.0091557134,0.0049407406,0.0088426706,0.0062579706];
  const ORF43=[0,0.0006565988,0.0015952328,0.0030065221,0.0126531621,0.0159805384,0.0073114854,0.0134076831,0.0106841704];
  lineChart(s,0.72,1.92,8.15,4.34,[
    {name:'ORF5 · E',values:norm(ORF5)},
    {name:'ORF31 · E',values:norm(ORF31)},
    {name:'ORF22 · L',values:norm(ORF22)},
    {name:'ORF43 · L',values:norm(ORF43)},
  ],times,{maxVal:1.05,majorUnit:0.25,valTitle:'within-gene normalized abundance',colors:[C.orange,'E3A54E',C.green,C.blue]});
  addPill(s,'EARLY-DOMINANT',9.35,2.05,2.25,C.orange,C.paleOrange);
  s.addText('ORF5 / ORF31\nstrong activation by 2–4 hpi', {x:9.35,y:2.52,w:2.8,h:0.80,fontSize:16.5,bold:true,color:C.navy,margin:0});
  addPill(s,'LATE-DOMINANT',9.35,3.58,2.25,C.green,'EEF7EA');
  s.addText('ORF22 / ORF43\nrise through 8–12 hpi', {x:9.35,y:4.05,w:2.8,h:0.80,fontSize:16.5,bold:true,color:C.navy,margin:0});
  s.addText('The broad IE → E → L logic is real.', {x:9.35,y:5.35,w:2.9,h:0.55,fontSize:20,bold:true,color:C.teal,margin:0});
  s.addText('But it is only the first layer.', {x:9.35,y:5.93,w:2.9,h:0.35,fontSize:15.5,color:C.gray,margin:0});
  addSource(s, 'Representative values redrawn from repository canonical TSS+TES abundance table; timing agrees with published kinetic analysis.');
  addNote(s, 'When we require reads to span both the canonical TSS and TES, the expected temporal waves are visible. Early genes such as ORF5 and ORF31 are strongest around two to four hours, whereas representative late genes rise around eight to twelve hours. So the classical model is not wrong — the question is how much complexity it hides.');
}

// 9. Ends alone mislead
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Why TSS-only or TES-only kinetics can be misleading', 'Transcript architecture');
  s.addText('In tandem herpesvirus loci, several RNAs can share one end while carrying different upstream structure.', {x:0.75,y:1.28,w:10.9,h:0.42,fontSize:18,color:C.gray,margin:0});
  // genes
  geneArrow(s,1.00,2.10,1.35,'A',C.blue,false);
  geneArrow(s,2.22,2.10,1.35,'B',C.teal,false);
  geneArrow(s,3.44,2.10,1.35,'C',C.orange,false);
  geneArrow(s,4.66,2.10,1.35,'D',C.red,false);
  s.addText('co-terminal gene cluster', {x:1.0,y:2.58,w:5.0,h:0.25,fontSize:11,color:C.gray,align:'center',margin:0});
  transcript(s,0.95,3.20,5.4,'A-B-C-D RNA',C.blue,0.02,0.98);
  transcript(s,0.95,3.65,5.4,'B-C-D RNA',C.teal,0.25,0.98);
  transcript(s,0.95,4.10,5.4,'C-D RNA',C.orange,0.49,0.98);
  transcript(s,0.95,4.55,5.4,'D RNA',C.red,0.72,0.98);
  // ambiguity panels
  s.addText('TSS signal', {x:7.00,y:2.02,w:2.0,h:0.34,fontSize:16,bold:true,color:C.navy,margin:0});
  s.addText('Promoter-specific\nbut does not identify which TES is used', {x:7.00,y:2.48,w:2.35,h:1.05,fontSize:15,color:C.gray,margin:0.03,fill:{color:C.paleOrange},line:{color:'F6D6AF'},radius:0.08});
  s.addText('TES signal', {x:9.62,y:2.02,w:2.0,h:0.34,fontSize:16,bold:true,color:C.navy,margin:0});
  s.addText('Can sum several\ntranscripts from different promoters', {x:9.62,y:2.48,w:2.35,h:1.05,fontSize:15,color:C.gray,margin:0.03,fill:{color:C.paleBlue},line:{color:'C8DCEB'},radius:0.08});
  s.addShape(pptx.ShapeType.line,{x:8.94,y:3.95,w:0.9,h:0,line:{color:C.teal,width:2.6,endArrowType:'triangle'}});
  s.addText('Full-length read', {x:7.58,y:4.15,w:3.0,h:0.38,fontSize:17,bold:true,color:C.teal,align:'center',margin:0});
  s.addText('links the actual 5′ start, splice pattern and 3′ end on the same molecule', {x:7.10,y:4.65,w:4.15,h:0.95,fontSize:20,bold:true,color:C.navy,align:'center',margin:0.03});
  s.addText('This is why long reads are essential for interpreting viral expression dynamics.', {x:7.10,y:5.83,w:4.15,h:0.60,fontSize:15.5,color:C.gray,align:'center',margin:0});
  addSource(s, 'Published interpretation: co-terminal and overlapping RNAs explain mismatches between TSS and TES kinetics.');
  addNote(s, 'This is an important conceptual point for the audience. A transcription end may be shared by several RNAs initiated from different promoters. Conversely, one promoter can generate alternative transcript ends. Measuring ends independently therefore mixes molecules. Long reads let us link the two ends and any introns on the same RNA, which is the only way to assign kinetics to a specific transcript structure.');
}

// 10. continuous landscape
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Data-driven clustering reveals a continuum, not three perfectly separated boxes', 'Beyond IE / E / L');
  s.addText('Most genes follow their expected phase — but several bridge or violate the traditional assignments.', {x:0.73,y:1.28,w:10.7,h:0.40,fontSize:18,color:C.gray,margin:0});
  // gradient waves
  const wy=2.03;
  s.addShape(pptx.ShapeType.arc,{x:0.92,y:wy,w:3.6,h:2.45,adjustPoint:0.35,rotate:0,fill:{color:C.orange,transparency:55},line:{color:C.orange,transparency:100}});
  s.addShape(pptx.ShapeType.arc,{x:3.24,y:wy+0.18,w:4.6,h:2.55,adjustPoint:0.35,rotate:0,fill:{color:C.teal,transparency:58},line:{color:C.teal,transparency:100}});
  s.addShape(pptx.ShapeType.arc,{x:6.42,y:wy+0.16,w:5.1,h:2.75,adjustPoint:0.35,rotate:0,fill:{color:C.blue,transparency:62},line:{color:C.blue,transparency:100}});
  s.addText('early', {x:1.10,y:2.30,w:1.0,h:0.3,fontSize:13,bold:true,color:C.orange,margin:0});
  s.addText('intermediate / mixed', {x:4.40,y:2.52,w:2.1,h:0.3,fontSize:13,bold:true,color:C.teal,margin:0});
  s.addText('late', {x:9.55,y:2.60,w:1.0,h:0.3,fontSize:13,bold:true,color:C.blue,margin:0});
  // time axis
  s.addShape(pptx.ShapeType.line,{x:1.10,y:4.46,w:10.20,h:0,line:{color:C.navy,width:1.6,endArrowType:'triangle'}});
  ['1','2','4','6','8','12','18','24','48 h'].forEach((t,i)=>s.addText(t,{x:0.98+i*1.28,y:4.58,w:0.55,h:0.24,fontSize:10.5,color:C.gray,align:'center',margin:0}));
  const examples=[
    ['ORF38','traditionally L','TSS peaks at 6 hpi','L → early-like',C.red,2.2],
    ['ORF54','traditionally E','TSS peak at 24 hpi','E → late-like',C.orange,7.7],
    ['ORF45','traditionally L','peaks at 12 and 48 hpi','bimodal late',C.blue,9.5],
  ];
  examples.forEach((e,i)=>{
    const x=0.95+i*4.05;
    s.addText(e[0],{x,y:5.20,w:1.15,h:0.40,fontSize:20,bold:true,color:e[4],margin:0});
    s.addText(e[1]+'\n'+e[2],{x:x+1.20,y:5.20,w:2.45,h:0.74,fontSize:12.8,color:C.gray,margin:0});
    s.addText(e[3],{x,y:6.12,w:3.4,h:0.34,fontSize:12,bold:true,color:C.navy,align:'center',margin:0,fill:{color:C.white},line:{color:C.light},radius:0.10});
  });
  addSource(s, 'Published examples from TSS kinetics and de novo clustering: ORF38, ORF54 and ORF45 illustrate mixed temporal behavior.');
  addNote(s, 'The de novo clustering makes the central biological message clearer. Most genes still align with early or late behavior, but the borders are fuzzy. ORF38 is traditionally late yet shows an early-like TSS peak; ORF54 is classified as early but has a late TSS peak; and ORF45 shows a bimodal pattern. The transcriptome behaves as overlapping temporal waves rather than three sealed compartments.');
}

// 11. Splicing dynamics
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Splicing itself changes through the infection cycle', 'Dynamic isoform usage');
  s.addText('ORF9 provides a clear example: splice isoforms are absent early and accumulate late.', {x:0.74,y:1.28,w:9.7,h:0.40,fontSize:18,color:C.gray,margin:0});
  // ORF9 structure
  geneArrow(s,0.95,2.02,1.8,'ORF9',C.teal,false);
  transcript(s,0.92,2.62,4.65,'canonical',C.teal,0.02,0.94);
  s.addShape(pptx.ShapeType.line,{x:1.30,y:3.12,w:0.72,h:0,line:{color:C.orange,width:2.4}});
  s.addShape(pptx.ShapeType.line,{x:2.05,y:3.12,w:1.85,h:0,line:{color:C.orange,width:2.4,endArrowType:'triangle'}});
  s.addText('TR134 / TR172 (shared intron)',{x:1.15,y:3.30,w:3.2,h:0.25,fontSize:10.5,color:C.orange,margin:0});
  // authoritative milestones
  s.addText('Combined spliced isoforms', {x:5.45,y:1.95,w:2.5,h:0.30,fontSize:14.5,bold:true,color:C.navy,margin:0});
  const xs=[5.55,6.75,7.95,9.15,10.35]; const labs=['1–8 h','12 h','18 h','24 h','48 h']; const vals=['0%','1.03%','11.49%','12.85%','17.10%'];
  s.addShape(pptx.ShapeType.line,{x:5.55,y:3.00,w:5.15,h:0,line:{color:C.mid,width:1.6}});
  xs.forEach((x,i)=>{
    s.addShape(pptx.ShapeType.ellipse,{x:x-0.10,y:2.90,w:0.20,h:0.20,fill:{color:i<2?C.orange:C.red},line:{color:C.white}});
    s.addText(labs[i],{x:x-0.40,y:3.18,w:0.80,h:0.24,fontSize:10.5,color:C.gray,align:'center',margin:0});
    s.addText(vals[i],{x:x-0.45,y:2.42,w:0.90,h:0.28,fontSize:13.5,bold:true,color:i<2?C.orange:C.red,align:'center',margin:0});
  });
  // canonical decline callout
  s.addText('Canonical ORF9', {x:1.00,y:4.34,w:2.1,h:0.30,fontSize:14,bold:true,color:C.navy,margin:0});
  s.addText('~100%',{x:1.00,y:4.78,w:1.25,h:0.55,fontSize:28,bold:true,color:C.teal,margin:0});
  s.addShape(pptx.ShapeType.line,{x:2.30,y:5.05,w:2.15,h:0,line:{color:C.mid,width:3,endArrowType:'triangle'}});
  s.addText('41.69%',{x:4.62,y:4.78,w:1.45,h:0.55,fontSize:28,bold:true,color:C.red,margin:0});
  s.addText('early (2–4 hpi)',{x:0.98,y:5.36,w:1.6,h:0.26,fontSize:10.8,color:C.gray,margin:0});
  s.addText('48 hpi',{x:4.62,y:5.36,w:1.0,h:0.26,fontSize:10.8,color:C.gray,margin:0});
  s.addText('The late transcriptome is not simply “more RNA” — it contains a different mixture of RNA structures.', {x:6.15,y:4.52,w:5.5,h:1.12,fontSize:21,bold:true,color:C.navy,align:'center',valign:'mid',margin:0.04,fill:{color:C.paleTeal},line:{color:'CBE7E6'},radius:0.12});
  addSource(s, 'Published ORF9 ratios: TR134+TR172 rise from 0% at 1–8 hpi to 17.10% at 48 hpi; canonical isoform declines late.');
  addNote(s, 'Splicing shows its own kinetics. For ORF9, the two highlighted spliced transcripts are essentially absent from one through eight hours. They appear at 12 hours, exceed eleven percent by 18 hours and reach about 17 percent at 48 hours. In parallel, the canonical ORF9 RNA loses dominance. So the identity of the transcript pool changes over time.');
}

// 12. Isoform switching
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Isoform switching can invert which RNA is dominant', 'Dynamic transcript architecture');
  s.addText('ORF54 is an especially clean switch between a short isoform and the canonical transcript.', {x:0.74,y:1.28,w:10.2,h:0.40,fontSize:18,color:C.gray,margin:0});
  const times=['1','2','4','6','8','12','18','24','48'];
  const canon=[0,0,0,0.2777778,0.4884780,0.8117061,0.9636016,0.9257059,0.9316761];
  const short=[0,1,0.9761905,0.6,0.2247037,0.0148631,0.0107202,0.0012667,0.0084559];
  lineChart(s,0.72,1.90,7.30,4.43,[{name:'ORF54 canonical',values:canon},{name:'ORF54-S',values:short}],times,{maxVal:1.05,majorUnit:0.25,valTitle:'fraction of ORF54 isoforms',colors:[C.teal,C.orange]});
  s.addText('2 hpi', {x:8.48,y:2.00,w:1.0,h:0.26,fontSize:11,bold:true,color:C.orange,margin:0});
  s.addText('ORF54-S ≈ 100%', {x:8.48,y:2.36,w:3.1,h:0.42,fontSize:21,bold:true,color:C.orange,margin:0});
  s.addShape(pptx.ShapeType.line,{x:8.50,y:3.04,w:2.5,h:0,line:{color:C.mid,width:3,endArrowType:'triangle'}});
  s.addText('18 hpi', {x:8.48,y:3.35,w:1.0,h:0.26,fontSize:11,bold:true,color:C.teal,margin:0});
  s.addText('canonical ≈ 96%', {x:8.48,y:3.72,w:3.2,h:0.42,fontSize:21,bold:true,color:C.teal,margin:0});
  s.addText('A second example: ORF19', {x:8.48,y:4.62,w:3.0,h:0.30,fontSize:14.5,bold:true,color:C.navy,margin:0});
  s.addText('canonical: 100% early → 14.59% at 48 hpi\ncomplex isoforms: → 72.21% at 48 hpi', {x:8.48,y:5.02,w:3.65,h:0.85,fontSize:15.5,color:C.gray,margin:0});
  s.addText('Same gene. Different RNA architecture. Different stage.', {x:8.48,y:6.00,w:3.6,h:0.54,fontSize:17.5,bold:true,color:C.navy,margin:0});
  addSource(s, 'ORF54 curves redrawn from repository mean isoform ratios; ORF19 endpoint values from the published article.');
  addNote(s, 'ORF54 is one of the strongest visual examples for a talk. The short isoform dominates extremely early, but by eight hours the canonical transcript overtakes it and by 18 hours the canonical RNA is over 96 percent of the detected ORF54 isoform pool. ORF19 shows another dramatic shift: the canonical RNA falls from complete early dominance to about 15 percent by 48 hours, while complex isoforms become dominant.');
}

// 13. overlaps
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Transcriptional overlap becomes progressively denser', 'Genome-wide architecture');
  s.addText('Raw long reads reveal increasing overlap as more promoters and transcript isoforms become active.', {x:0.74,y:1.28,w:10.6,h:0.40,fontSize:18,color:C.gray,margin:0});
  const groups=[{y:2.05,label:'2 hpi',n:3},{y:3.48,label:'8 hpi',n:6},{y:4.91,label:'24–48 hpi',n:9}];
  groups.forEach((g,gi)=>{
    s.addText(g.label,{x:0.85,y:g.y+0.2,w:1.05,h:0.3,fontSize:12,bold:true,color:gi===0?C.orange:(gi===1?C.teal:C.blue),margin:0});
    // gene backbone
    for(let j=0;j<5;j++) geneArrow(s,2.05+j*1.25,g.y,1.32,String.fromCharCode(65+j),j%2?C.teal:C.blue,j===3);
    // transcripts
    for(let k=0;k<g.n;k++){
      const start=(k%5)*0.12; const end=Math.min(0.98,start+0.38+(k%4)*0.12);
      const col=[C.orange,C.teal,C.blue,C.red,C.purple][k%5];
      transcript(s,2.05,g.y+0.48+(k%3)*0.22,6.25,'',col,start,end, false);
    }
  });
  s.addText('Potential consequence', {x:9.20,y:2.05,w:2.5,h:0.32,fontSize:13,bold:true,color:C.teal,margin:0});
  s.addText('Convergent and divergent transcription can create collisions or promoter interference.', {x:9.20,y:2.48,w:3.15,h:1.18,fontSize:19.5,bold:true,color:C.navy,margin:0.03});
  s.addText('The study proposes transcriptional interference as one possible regulatory layer — an interpretation that still requires direct mechanistic testing.', {x:9.20,y:4.12,w:3.15,h:1.35,fontSize:15.5,color:C.gray,margin:0.02});
  addPill(s,'HYPOTHESIS',9.20,5.78,1.42,C.red,C.paleRed);
  addSource(s, 'Published observation: transcriptional overlaps increase as infection progresses; regulatory interference is a proposed interpretation.');
  addNote(s, 'The raw-read view emphasizes another property of herpesvirus transcription: overlaps become more extensive as infection progresses. This means neighboring polymerases can potentially interfere with one another. The paper discusses transcriptional interference as a plausible regulatory mechanism, but I would present that explicitly as a hypothesis rather than a demonstrated mechanism.');
}

// 14. model
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'EHV-1 regulation operates at multiple transcript-level layers', 'Integrated model');
  s.addText('Time changes not only how much a gene is transcribed, but which RNA molecule is produced.', {x:0.74,y:1.28,w:10.8,h:0.42,fontSize:19,bold:true,color:C.navy,margin:0});
  const cards=[
    {x:0.83,title:'Promoter choice',tag:'TSS',color:C.orange,fill:C.paleOrange,text:'Alternative initiation changes 5′ UTRs and can expose embedded ORFs.'},
    {x:3.83,title:'Termination choice',tag:'TES',color:C.blue,fill:C.paleBlue,text:'Alternative TES usage and co-termination change transcript length and overlap.'},
    {x:6.83,title:'Splicing',tag:'INTRONS',color:C.purple,fill:'F4ECF8',text:'Splice isoform fractions shift strongly with infection stage.'},
    {x:9.83,title:'Overlap',tag:'TOPOLOGY',color:C.teal,fill:C.paleTeal,text:'Long, complex and antisense RNAs reshape the transcriptional neighborhood.'},
  ];
  cards.forEach(c=>{
    s.addText(c.tag,{x:c.x,y:2.10,w:1.30,h:0.34,fontSize:10.5,bold:true,color:c.color,align:'center',valign:'mid',margin:0,fill:{color:c.fill},line:{color:c.fill},radius:0.15});
    s.addText(c.title,{x:c.x,y:2.62,w:2.55,h:0.40,fontSize:20,bold:true,color:C.navy,margin:0});
    s.addText(c.text,{x:c.x,y:3.18,w:2.55,h:1.25,fontSize:15.2,color:C.gray,margin:0.02});
    s.addShape(pptx.ShapeType.line,{x:c.x+1.27,y:4.70,w:0,h:0.70,line:{color:c.color,width:2.5,endArrowType:'triangle'}});
  });
  s.addText('TIME-RESOLVED TRANSCRIPT ARCHITECTURE', {x:2.40,y:5.60,w:8.55,h:0.72,fontSize:24,bold:true,color:C.white,align:'center',valign:'mid',margin:0,fill:{color:C.navy},radius:0.18,shadow});
  s.addText('A more informative unit of viral regulation than gene counts alone', {x:3.08,y:6.38,w:7.2,h:0.34,fontSize:14.2,color:C.gray,align:'center',margin:0});
  addSource(s, 'Synthesis of the published findings: alternative TSS/TES, splicing, isoform switching and overlap all vary through infection.');
  addNote(s, 'This is the conceptual synthesis slide. Gene abundance is only one layer. Promoter choice changes 5-prime structure; termination changes transcript length and overlap; splicing changes exon composition; and overlapping or antisense transcription changes the genomic environment itself. Long reads make it possible to study these layers together on individual RNA molecules.');
}

// 15. Take home
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Take-home messages', 'Conclusion');
  const msgs=[
    ['1','Integrated CAGE + long reads refine viral transcript boundaries','169 novel transcripts and stringent multi-platform validation.',C.teal,C.paleTeal],
    ['2','ORF64 is uniquely immediate-early','The beginning of the cascade is discrete; later expression is not.',C.red,C.paleRed],
    ['3','The EHV-1 transcriptome is dynamically re-wired','Splicing, TSS/TES choice, isoform switching and overlap change across infection.',C.blue,C.paleBlue],
  ];
  msgs.forEach((m,i)=>{
    const y=1.52+i*1.58;
    s.addText(m[0],{x:0.85,y,w:0.62,h:0.62,fontSize:25,bold:true,color:C.white,align:'center',valign:'mid',margin:0,fill:{color:m[3]},radius:0.31});
    s.addText(m[1],{x:1.75,y:y-0.01,w:7.7,h:0.38,fontSize:20.5,bold:true,color:C.navy,margin:0});
    s.addText(m[2],{x:1.75,y:y+0.46,w:8.65,h:0.48,fontSize:15.3,color:C.gray,margin:0});
  });
  s.addText('The central message', {x:9.75,y:1.68,w:2.35,h:0.28,fontSize:11,bold:true,color:C.teal,align:'center',margin:0,charSpacing:1.2});
  s.addText('EHV-1 gene regulation is a continuous, overlapping temporal landscape — not simply three kinetic boxes.', {x:9.38,y:2.12,w:3.15,h:2.05,fontSize:23,bold:true,color:C.white,align:'center',valign:'mid',margin:0.10,fill:{color:C.navy},radius:0.14,shadow});
  s.addText('Questions?', {x:9.75,y:4.72,w:2.4,h:0.55,fontSize:27,bold:true,color:C.teal,align:'center',margin:0});
  s.addText('PLOS ONE 20(4): e0320439\ngithub.com/Balays/EHV-1-dynamic', {x:9.55,y:5.48,w:2.8,h:0.58,fontSize:11.2,color:C.gray,align:'center',margin:0});
  addSource(s, 'Sequencing data: ENA PRJEB52190 and PRJEB6233 · analysis code: Balays/EHV-1-dynamic.');
  addNote(s, 'Close with three points. First, integrating CAGE with long reads materially improves transcript annotation. Second, ORF64 is a clean immediate-early outlier, whereas the downstream program has overlapping waves. Third, the strongest new biological insight is transcript-level regulation: the proportions of different RNAs from the same locus change through time.');
}

// 16. Backup validation
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Backup: acceptance criteria were deliberately stringent', 'Methods detail');
  const rows=[
    ['dcDNA evidence','≥3 reads sharing transcript boundaries','correct 5′ adapter / strand evidence'],
    ['CAGE evidence','5′ end within ±10 nt of refined TSS peak','confidence from CAGE support + score'],
    ['TES evidence','3′ end within ±10 nt of validated TES','poly(A) and false-priming filters'],
    ['Embedded mRNAs','dRNA/NAGATA support within 25 nt','TSS signal ≥5% of canonical transcript'],
  ];
  rows.forEach((r,i)=>{
    const y=1.50+i*1.24;
    s.addText(r[0],{x:0.82,y,w:2.05,h:0.86,fontSize:17,bold:true,color:C.white,align:'center',valign:'mid',margin:0.05,fill:{color:[C.teal,C.orange,C.blue,C.red][i]},radius:0.08});
    s.addText(r[1],{x:3.10,y,w:4.15,h:0.86,fontSize:15.2,color:C.navy,align:'center',valign:'mid',margin:0.07,fill:{color:C.white},line:{color:C.light},radius:0.08});
    s.addText(r[2],{x:7.48,y,w:4.55,h:0.86,fontSize:15.2,color:C.gray,align:'center',valign:'mid',margin:0.07,fill:{color:C.white},line:{color:C.light},radius:0.08});
  });
  s.addText('Why this matters', {x:0.85,y:6.40,w:1.8,h:0.28,fontSize:11.5,bold:true,color:C.teal,margin:0});
  s.addText('5′ truncation can arise from sequencing/library artifacts or biological recapping, so independent evidence is essential before calling short isoforms genuine viral TSSs.', {x:2.55,y:6.31,w:9.65,h:0.48,fontSize:14.2,color:C.gray,margin:0});
  addSource(s, 'Backup slide based on the final revised Methods and published filtering rationale.');
  addNote(s, 'Use this only if methods questions come up. The main point is that the pipeline does not accept every long-read boundary as a transcript. Adapter evidence, independent CAGE support, TES validation and extra filtering for potentially truncated ORF-carrying RNAs were combined to reduce false positives.');
}

// 17. Backup resources
{
  const s = pptx.addSlide('MASTER');
  addTitle(s, 'Backup: data, code and reproducibility', 'Resources');
  const items=[
    ['Published article','Tombácz et al. (2025)\nPLOS ONE 20(4): e0320439','doi:10.1371/journal.pone.0320439',C.orange,C.paleOrange],
    ['Raw sequencing data','European Nucleotide Archive','PRJEB52190 · PRJEB6233',C.blue,C.paleBlue],
    ['Analysis repository','Balays/EHV-1-dynamic','R/Rmd workflows, processed tables and publication figures',C.teal,C.paleTeal],
  ];
  items.forEach((it,i)=>{
    const y=1.62+i*1.52;
    s.addText(it[0],{x:0.90,y,w:2.25,h:0.96,fontSize:17.5,bold:true,color:C.navy,align:'center',valign:'mid',margin:0.08,fill:{color:it[4]},line:{color:it[3]},radius:0.08});
    s.addText(it[1],{x:3.45,y:y+0.05,w:4.0,h:0.72,fontSize:17,bold:true,color:C.navy,margin:0.02});
    s.addText(it[2],{x:7.70,y:y+0.05,w:4.25,h:0.72,fontSize:14.2,color:C.gray,margin:0.02});
  });
  s.addText('Repository status used for this deck', {x:0.90,y:6.26,w:2.7,h:0.28,fontSize:11.5,bold:true,color:C.teal,margin:0});
  s.addText('main branch · latest checked repository commit before deck creation: df4589e (17 Mar 2025)', {x:3.58,y:6.23,w:8.2,h:0.34,fontSize:13.5,color:C.gray,margin:0});
  addSource(s, 'The deck is grounded in the published paper and the analysis repository used to generate the publication results.');
  addNote(s, 'This backup slide is useful if someone asks where the data or code can be found. The paper provides the ENA accessions and the GitHub repository contains the workflows, processed result tables and generated publication figures.');
}


fs.mkdirSync('presentation', { recursive: true });
pptx.writeFile({ fileName: 'presentation/EHV1_dynamic_transcriptome_conference_master.pptx' });
// rebuild trigger 2026-09-15
