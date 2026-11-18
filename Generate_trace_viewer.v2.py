#!/usr/bin/env python3
"""
mTracePatho Phylogenetic Tracing Viewer Generator
==================================================
Generates an interactive HTML phylogenetic tree viewer for mTracePatho tracing results.

Usage:
    python generate_trace_viewer.py \
        --tree Jaccard_fastme_input.nwk \
        --meta Klebsiella_pneumoniae-gene.tsv \
        --trace trace_result.tsv \
        --taxname Klebsiella_pneumoniae \
        --dist-unit "Jaccard distance" \
        --output mTracePatho_viewer.html
"""

import argparse
import csv
import json
import re
import sys
import os


def read_tree(path):
    with open(path, 'r') as f:
        return f.read().strip()


def read_meta(path):
    meta = {}
    with open(path, 'r', newline='') as f:
        content = f.read().replace('\r\n', '\n').replace('\r', '\n')
    lines = content.strip().split('\n')
    reader = csv.DictReader(lines, delimiter='\t')
    cols = reader.fieldnames
    for row in reader:
        strain = row.get('strain', '')
        if strain:
            meta[strain] = dict(row)
    return meta, cols


def read_trace(path):
    with open(path, 'r') as f:
        return f.read().strip()


def format_taxname(name):
    """Convert Klebsiella_pneumoniae -> Klebsiella pneumoniae"""
    return name.replace('_', ' ')


def classify_columns(cols):
    """Auto-detect typing columns vs gene columns."""
    typing_keywords = ['type', 'mlst', 'phylo', 'group', 'lineage', 'clade', 'species', 'ST']
    gene_cols = []
    typing_cols = []
    for c in cols:
        if c == 'strain':
            continue
        is_typing = any(kw.lower() in c.lower() for kw in typing_keywords)
        if is_typing:
            typing_cols.append(c)
        else:
            gene_cols.append(c)
    return typing_cols, gene_cols


def classify_genes(gene_cols):
    """Classify gene columns into virulence vs resistance."""
    resistance_kw = [
        'bla', 'oqx', 'sul', 'fos', 'aad', 'arr', 'cat', 'rmt',
        'qnr', 'tet', 'mcr', 'cfr',
        'arm', 'aac', 'ant', 'aph', 'str', 'dfrA',
        'ere', 'erm', 'mef', 'mph', 'msr', 'lnu', 'vat',
        'van', 'optrA', 'cml', 'flor',
        'ble', 'nim', 'gar',
    ]
    vir_genes = []
    res_genes = []
    for g in gene_cols:
        if any(g.lower().startswith(kw) for kw in resistance_kw):
            res_genes.append(g)
        else:
            vir_genes.append(g)
    return vir_genes, res_genes


def generate_html(tree_nwk, meta_json, trace_result, taxname, typing_cols, vir_genes, res_genes, dist_unit='distance'):
    """Generate the complete HTML file content."""
    
    display_name = format_taxname(taxname)
    
    # Build color-by options from typing columns
    color_options = ''
    for col in typing_cols:
        label = col.replace('mlsttype', 'MLST').replace('_', ' ')
        color_options += f'        <option value="{col}">{label}</option>\n'
    color_options += '        <option value="trace">Trace Result</option>\n'

    vir_genes_js = json.dumps(vir_genes)
    res_genes_js = json.dumps(res_genes)
    typing_cols_js = json.dumps(typing_cols)

    html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>mTracePatho - ''' + display_name + r''' Tracing Viewer</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@300;400;500;600&family=Inter:wght@300;400;500;600;700&display=swap');

:root {
  --bg-primary: #f7f7fb;
  --bg-secondary: #ffffff;
  --bg-tertiary: #f0f0f6;
  --bg-card: #ffffff;
  --text-primary: #1e1e2e;
  --text-secondary: #555570;
  --text-muted: #8e8ea8;
  --accent: #7c3aed;
  --accent-light: rgba(124,58,237,0.08);
  --accent-mid: rgba(124,58,237,0.15);
  --trace-color: #e03e5c;
  --trace-light: rgba(224,62,92,0.08);
  --trace-mid: rgba(224,62,92,0.18);
  --border: #e5e5f0;
  --border-light: #ededf5;
  --success: #16a34a;
  --shadow-sm: 0 1px 3px rgba(0,0,0,0.06);
  --shadow-md: 0 4px 16px rgba(0,0,0,0.08);
  --shadow-lg: 0 8px 32px rgba(0,0,0,0.10);
  --radius: 10px;
}

* { margin:0; padding:0; box-sizing:border-box; }
html, body { height:100%; overflow:hidden; background:var(--bg-primary); color:var(--text-primary); font-family:'Inter',system-ui,sans-serif; }

#app { display:flex; height:100vh; }

/* ── Sidebar ── */
#sidebar { width:340px; min-width:340px; background:var(--bg-secondary); border-right:1px solid var(--border); display:flex; flex-direction:column; z-index:10; overflow-y:auto; box-shadow:var(--shadow-sm); }
#sidebar-header { padding:22px 24px 16px; border-bottom:1px solid var(--border); }
#sidebar-header .logo-row { display:flex; align-items:center; gap:10px; }
#sidebar-header .logo-icon { width:36px; height:36px; background:var(--accent); border-radius:8px; display:flex; align-items:center; justify-content:center; color:#fff; font-size:18px; font-weight:700; }
#sidebar-header h1 { font-size:22px; font-weight:700; letter-spacing:-0.4px; color:var(--accent); }
#sidebar-header h1 span { color:var(--accent); }
#sidebar-header .taxname { font-size:15px; color:var(--accent); margin-top:6px; font-style:italic; font-weight:500; }
#sidebar-header .sub { font-size:12px; color:var(--text-muted); margin-top:2px; font-family:'JetBrains Mono',monospace; }

.sidebar-section { padding:16px 24px; border-bottom:1px solid var(--border-light); }
.sidebar-section h3 { font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:1.5px; color:var(--text-muted); margin-bottom:10px; }

/* ── Trace Card ── */
.trace-card { background:linear-gradient(135deg, var(--trace-light), rgba(224,62,92,0.03)); border:1px solid var(--trace-mid); border-radius:var(--radius); padding:14px; }
.trace-card .label { font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:1px; color:var(--trace-color); margin-bottom:6px; }
.trace-card .value { font-family:'JetBrains Mono',monospace; font-size:13px; color:var(--text-primary); word-break:break-all; line-height:1.7; }
.trace-card .type-badge { display:inline-block; font-size:10px; font-weight:700; padding:3px 10px; border-radius:20px; margin-top:8px; }
.trace-card .type-badge.node { background:var(--trace-mid); color:var(--trace-color); }
.trace-card .type-badge.single { background:var(--accent-mid); color:var(--accent); }

/* ── Search ── */
#search-box { width:100%; padding:9px 12px; background:var(--bg-tertiary); border:1px solid var(--border); border-radius:8px; color:var(--text-primary); font-family:'JetBrains Mono',monospace; font-size:13px; outline:none; transition:border-color 0.2s, box-shadow 0.2s; }
#search-box:focus { border-color:var(--accent); box-shadow:0 0 0 3px var(--accent-light); }
#search-box::placeholder { color:var(--text-muted); }
#search-results { max-height:200px; overflow-y:auto; margin-top:6px; }
.search-result-item { padding:5px 8px; font-size:12px; font-family:'JetBrains Mono',monospace; cursor:pointer; border-radius:6px; color:var(--text-secondary); transition:background 0.15s; }
.search-result-item:hover { background:var(--accent-light); color:var(--accent); }

/* ── Controls ── */
select.ctrl-select { width:100%; padding:9px 12px; background:var(--bg-tertiary); border:1px solid var(--border); border-radius:8px; color:var(--text-primary); font-size:13px; outline:none; cursor:pointer; }
select.ctrl-select:focus { border-color:var(--accent); }

#legend { padding:8px 0; max-height:200px; overflow-y:auto; }
.legend-item { display:flex; align-items:center; gap:8px; padding:4px 0; font-size:12px; cursor:pointer; color:var(--text-secondary); transition:color 0.15s; }
.legend-item:hover { color:var(--accent); }
.legend-dot { width:11px; height:11px; border-radius:50%; flex-shrink:0; border:1px solid rgba(0,0,0,0.08); }

.stat-row { display:flex; justify-content:space-between; padding:4px 0; font-size:13px; }
.stat-row .stat-label { color:var(--text-muted); }
.stat-row .stat-value { font-family:'JetBrains Mono',monospace; color:var(--text-primary); font-weight:500; }

.btn { padding:7px 14px; border-radius:8px; border:1px solid var(--border); background:var(--bg-secondary); color:var(--text-secondary); font-size:12px; font-weight:600; cursor:pointer; transition:all 0.2s; }
.btn:hover { border-color:var(--accent); color:var(--accent); background:var(--accent-light); }
.btn.active { background:var(--accent); border-color:var(--accent); color:#fff; }
.btn-group { display:flex; gap:6px; flex-wrap:wrap; }

/* ── Main Area ── */
#main { flex:1; position:relative; overflow:hidden; background:var(--bg-primary); }
#tree-container { width:100%; height:100%; }
#tree-container svg { width:100%; height:100%; }

/* ── Tooltip ── */
#tooltip { position:fixed; pointer-events:none; background:var(--bg-card); border:1px solid var(--border); border-radius:var(--radius); padding:14px 16px; max-width:400px; font-size:13px; z-index:100; opacity:0; transition:opacity 0.15s; box-shadow:var(--shadow-lg); }
#tooltip.visible { opacity:1; }
#tooltip .tt-name { font-family:'JetBrains Mono',monospace; font-size:14px; font-weight:600; color:var(--accent); margin-bottom:8px; }
#tooltip .tt-row { display:flex; gap:8px; padding:2px 0; }
#tooltip .tt-key { color:var(--text-muted); min-width:70px; font-weight:500; }
#tooltip .tt-val { color:var(--text-primary); font-family:'JetBrains Mono',monospace; }
#tooltip .tt-genes { display:flex; flex-wrap:wrap; gap:3px; margin-top:6px; }
#tooltip .tt-gene { font-size:10px; padding:2px 6px; border-radius:4px; font-family:'JetBrains Mono',monospace; }
#tooltip .tt-gene.present { background:var(--accent-light); color:var(--accent); }
#tooltip .tt-gene.absent { background:var(--bg-tertiary); color:var(--text-muted); }
#tooltip .tt-section { font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:0.8px; color:var(--text-muted); margin:8px 0 4px; }

/* ── Detail Panel ── */
#detail-panel { position:fixed; width:370px; background:var(--bg-card); border:1px solid var(--border); border-radius:12px; padding:16px 22px 22px; z-index:20; display:none; box-shadow:var(--shadow-lg); max-height:calc(100vh - 24px); overflow-y:auto; transition:opacity 0.15s; }
#detail-panel.visible { display:block; }
.dp-toolbar { display:flex; align-items:center; justify-content:flex-end; gap:8px; margin:-16px -22px 8px; padding:8px 12px 6px 16px; cursor:grab; border-radius:12px 12px 0 0; border-bottom:1px solid var(--border-light); user-select:none; }
.dp-toolbar:active { cursor:grabbing; }
.dp-toolbar .drag-hint { flex:1; font-size:10px; color:var(--text-muted); letter-spacing:3px; text-align:center; pointer-events:none; }
.dp-opacity-label { display:flex; align-items:center; gap:4px; font-size:16px; color:var(--text-muted); cursor:pointer; user-select:none; }
.dp-opacity-label input[type=range] { width:70px; height:4px; accent-color:var(--accent); cursor:pointer; }
#detail-panel .close-btn { background:none; border:none; color:var(--text-muted); cursor:pointer; font-size:22px; line-height:1; transition:color 0.15s; position:static; }
#detail-panel .close-btn:hover { color:var(--text-primary); }
#detail-panel .dp-name { font-family:'JetBrains Mono',monospace; font-size:15px; font-weight:600; color:var(--accent); margin-bottom:14px; padding-right:28px; }
#detail-panel .dp-section { font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:1px; color:var(--text-muted); margin:14px 0 8px; padding-top:12px; border-top:1px solid var(--border-light); }
#detail-panel .dp-row { display:flex; justify-content:space-between; padding:3px 0; font-size:13px; }
#detail-panel .dp-key { color:var(--text-secondary); }
#detail-panel .dp-val { font-family:'JetBrains Mono',monospace; color:var(--text-primary); font-weight:500; }
#detail-panel .gene-grid { display:grid; grid-template-columns:repeat(3,1fr); gap:4px; margin-top:4px; }
#detail-panel .gene-chip { font-size:11px; padding:4px 8px; border-radius:6px; font-family:'JetBrains Mono',monospace; text-align:center; }
#detail-panel .gene-chip.present { background:var(--accent-light); color:var(--accent); border:1px solid var(--accent-mid); }
#detail-panel .gene-chip.absent { background:var(--bg-tertiary); color:var(--text-muted); border:1px solid var(--border-light); }
.dp-trace-badge { display:inline-block; padding:3px 12px; border-radius:20px; font-size:11px; font-weight:700; background:var(--trace-light); color:var(--trace-color); border:1px solid var(--trace-mid); margin-bottom:10px; }

/* ── Zoom Controls ── */
#zoom-controls { position:absolute; bottom:16px; right:16px; display:flex; flex-direction:column; gap:4px; z-index:20; }
#zoom-controls button { width:40px; height:40px; border-radius:var(--radius); background:var(--bg-card); border:1px solid var(--border); color:var(--text-secondary); font-size:20px; cursor:pointer; display:flex; align-items:center; justify-content:center; box-shadow:var(--shadow-sm); transition:all 0.15s; }
#zoom-controls button:hover { border-color:var(--accent); color:var(--accent); }

/* ── Info Bar ── */
#info-bar { position:absolute; bottom:16px; left:16px; display:flex; gap:10px; z-index:20; }
.info-chip { padding:5px 14px; background:var(--bg-card); border:1px solid var(--border); border-radius:20px; font-size:12px; color:var(--text-secondary); font-family:'JetBrains Mono',monospace; box-shadow:var(--shadow-sm); }

/* ── Fixed Scale Bar Overlay ── */
#scale-bar-overlay { position:absolute; top:14px; right:14px; z-index:15; background:rgba(255,255,255,0.92); border:1px solid var(--border); border-radius:10px; padding:10px 16px 8px; box-shadow:var(--shadow-sm); display:none; }
#scale-bar-overlay .sb-bar { display:flex; align-items:center; justify-content:center; }
#scale-bar-overlay .sb-line { height:2px; background:#444; position:relative; }
#scale-bar-overlay .sb-line::before, #scale-bar-overlay .sb-line::after { content:''; position:absolute; top:-4px; width:1.5px; height:10px; background:#444; }
#scale-bar-overlay .sb-line::before { left:0; }
#scale-bar-overlay .sb-line::after { right:0; }
#scale-bar-overlay .sb-label { text-align:center; font-family:'JetBrains Mono',monospace; font-size:10.5px; color:#444; font-weight:500; margin-top:5px; white-space:nowrap; }

/* ── Animations ── */
@keyframes pulse-trace { 0%,100% { r:6; opacity:1; } 50% { r:13; opacity:0.3; } }
.trace-pulse { animation: pulse-trace 2s ease-in-out infinite; }

/* ── Scrollbar ── */
::-webkit-scrollbar { width:5px; }
::-webkit-scrollbar-track { background:transparent; }
::-webkit-scrollbar-thumb { background:var(--border); border-radius:3px; }
::-webkit-scrollbar-thumb:hover { background:var(--text-muted); }

/* ── Loading ── */
#loading { position:fixed; inset:0; background:var(--bg-primary); display:flex; align-items:center; justify-content:center; z-index:1000; flex-direction:column; gap:16px; }
#loading .spinner { width:40px; height:40px; border:3px solid var(--border); border-top-color:var(--accent); border-radius:50%; animation:spin 0.8s linear infinite; }
@keyframes spin { to { transform:rotate(360deg); } }
#loading p { font-size:14px; color:var(--text-muted); }
</style>
</head>
<body>
<div id="loading"><div class="spinner"></div><p>Parsing phylogenetic tree...</p></div>
<div id="app">
  <div id="sidebar">
    <div id="sidebar-header">
      <div class="logo-row">
        <div class="logo-icon">&#119795;</div>
        <h1>mTrace<span>Patho</span></h1>
      </div>
      <div class="taxname">''' + display_name + r'''</div>
      <div class="sub">Phylogenetic Tracing Viewer</div>
    </div>
    <div class="sidebar-section">
      <h3>Tracing Result</h3>
      <div class="trace-card" id="trace-info"></div>
    </div>
    <div class="sidebar-section">
      <h3>Search Genome</h3>
      <input type="text" id="search-box" placeholder="GCF_xxx ; GCF_yyy (multi-search)" autocomplete="off">
      <div id="search-results"></div>
    </div>
    <div class="sidebar-section">
      <h3>Color By</h3>
      <select class="ctrl-select" id="color-select">
''' + color_options + r'''      </select>
      <div id="legend"></div>
    </div>
    <div class="sidebar-section">
      <h3>View Controls</h3>
      <div class="btn-group">
        <button class="btn active" id="btn-rect" onclick="setLayout('rect')">Rectangular</button>
        <button class="btn" id="btn-radial" onclick="setLayout('radial')">Radial</button>
      </div>
      <div class="btn-group" style="margin-top:8px;">
        <button class="btn" onclick="focusOnTrace()">Focus Trace</button>
        <button class="btn" onclick="resetView()">Reset View</button>
        <button class="btn" onclick="toggleLabels()">Toggle Labels</button>
      </div>
    </div>
    <div class="sidebar-section">
      <h3>Tree Statistics</h3>
      <div id="stats"></div>
    </div>
  </div>
  <div id="main">
    <div id="tree-container"></div>
    <div id="zoom-controls">
      <button onclick="zoomIn()">+</button>
      <button onclick="zoomOut()">&minus;</button>
      <button onclick="resetView()">&#x27F2;</button>
    </div>
    <div id="scale-bar-overlay"></div>
    <div id="info-bar">
      <div class="info-chip" id="info-tip-count"></div>
      <div class="info-chip" id="info-zoom"></div>
    </div>
    <div id="detail-panel">
      <div class="dp-toolbar" id="dp-drag-handle">
        <span class="drag-hint">&#x2261;&#x2261;&#x2261;</span>
        <label class="dp-opacity-label" title="Panel transparency">&#9681;
          <input type="range" id="dp-opacity-slider" min="20" max="100" value="95" oninput="setDetailOpacity(this.value)">
        </label>
        <button class="close-btn" onclick="closeDetail()">&times;</button>
      </div>
      <div id="detail-content"></div>
    </div>
  </div>
</div>
<div id="tooltip"></div>

<script>
// ============ EMBEDDED DATA ============
const NEWICK_STR = ''' + json.dumps(tree_nwk) + r''';
const META_DATA = ''' + meta_json + r''';
const TRACE_RAW = ''' + json.dumps(trace_result) + r''';
const TAX_NAME = ''' + json.dumps(display_name) + r''';
const TYPING_COLS = ''' + typing_cols_js + r''';
const VIR_GENES = ''' + vir_genes_js + r''';
const RES_GENES = ''' + res_genes_js + r''';
const DIST_UNIT = ''' + json.dumps(dist_unit) + r''';

// ============ GLOBALS ============
let treeRoot, allLeaves=[], allNodes=[];
let svg, g, zoom;
let currentLayout = 'rect';
let showLabels = true;
let tracedGenomes = TRACE_RAW.split(';').map(s=>s.trim()).filter(Boolean);
let isNodeTrace = tracedGenomes.length > 1;
let tipHeight = 16;
let lastDims = null; // store layout dims for scale bar

// ============ NEWICK PARSER (robust) ============
function parseNewick(rawStr) {
  let str = rawStr.replace(/;\s*$/, '').trim(); // strip trailing semicolons
  let i = 0;
  function skipNHX() { // skip [&NHX:...] or [...] annotations
    if (i < str.length && str[i] === '[') { while (i < str.length && str[i] !== ']') i++; i++; }
  }
  function readNode() {
    let node = { children: [], branchLength: 0, name: '' };
    if (str[i] === '(') {
      i++;
      while (true) {
        node.children.push(readNode());
        if (str[i] === ',') { i++; continue; }
        if (str[i] === ')') { i++; break; }
        break;
      }
    }
    skipNHX();
    // Read name — support 'quoted labels'
    let name = '';
    if (str[i] === "'") {
      i++; while (i < str.length && str[i] !== "'") { name += str[i]; i++; } i++;
    } else {
      while (i < str.length && str[i] !== ':' && str[i] !== ',' && str[i] !== ')' && str[i] !== ';' && str[i] !== '(' && str[i] !== '[') {
        name += str[i]; i++;
      }
    }
    node.name = name.trim();
    skipNHX();
    if (str[i] === ':') {
      i++;
      let bl = '';
      while (i < str.length && str[i] !== ',' && str[i] !== ')' && str[i] !== ';' && str[i] !== '(' && str[i] !== '[') {
        bl += str[i]; i++;
      }
      node.branchLength = Math.max(0, parseFloat(bl) || 0); // P2: clamp negative branch lengths
    }
    skipNHX();
    return node;
  }
  return readNode();
}

// ============ TREE UTILITIES ============
function assignIds(node, counter={val:0}) {
  node.id = counter.val++;
  node.isLeaf = node.children.length === 0;
  if (node.isLeaf) allLeaves.push(node);
  allNodes.push(node);
  node.children.forEach(c => { c.parent = node; assignIds(c, counter); });
}
function countLeaves(node) {
  if (node.isLeaf) { node.leafCount = 1; return 1; }
  node.leafCount = node.children.reduce((s,c) => s + countLeaves(c), 0);
  return node.leafCount;
}
function computeDepth(node, depth=0) {
  node.rootDist = depth;
  node.children.forEach(c => computeDepth(c, depth + c.branchLength));
}
function getMaxDepth(node) {
  if (node.isLeaf) return node.rootDist;
  return Math.max(...node.children.map(getMaxDepth));
}
function findLCA(root, names) {
  if (names.length === 1) {
    let found = null;
    function search(n) { if (n.name === names[0]) found = n; n.children.forEach(search); }
    search(root);
    return found;
  }
  let nameSet = new Set(names);
  function dfs(node) {
    if (node.isLeaf) return nameSet.has(node.name) ? node : null;
    let found = node.children.map(dfs).filter(Boolean);
    if (found.length > 1) return node;
    return found[0] || null;
  }
  return dfs(root);
}
function getDescendants(node) {
  let result = [];
  function walk(n) { result.push(n); n.children.forEach(walk); }
  walk(node);
  return result;
}

// ============ LAYOUT ============
function layoutRect(root) {
  let leafIdx = 0;
  let maxDist = getMaxDepth(root);
  let xScale = 800 / Math.max(maxDist, 0.001);
  function setLeafY(node) {
    if (node.isLeaf) { node.y = leafIdx * tipHeight; leafIdx++; return; }
    node.children.forEach(setLeafY);
    let ys = node.children.map(c => c.y);
    node.y = (Math.min(...ys) + Math.max(...ys)) / 2;
  }
  function setX(node) { node.x = node.rootDist * xScale; node.children.forEach(setX); }
  setLeafY(root); setX(root);
  return { width: maxDist * xScale + 300, height: leafIdx * tipHeight + 40, xScale: xScale, maxDist: maxDist };
}
function layoutRadial(root) {
  let leafIdx = 0, totalLeaves = allLeaves.length;
  let maxDist = getMaxDepth(root);
  let radius = computeRadialRadius();
  function setLeafAngle(node) {
    if (node.isLeaf) { node.angle = (leafIdx / totalLeaves) * 2 * Math.PI; leafIdx++; return; }
    node.children.forEach(setLeafAngle);
    let angles = node.children.map(c => c.angle);
    node.angle = (Math.min(...angles) + Math.max(...angles)) / 2;
  }
  function setPos(node) {
    let r = (node.rootDist / Math.max(maxDist, 0.001)) * radius;
    node.x = r * Math.cos(node.angle - Math.PI/2);
    node.y = r * Math.sin(node.angle - Math.PI/2);
    node.children.forEach(setPos);
  }
  setLeafAngle(root); setPos(root);
  let xScaleR = radius / Math.max(maxDist, 0.001);
  return { width: radius * 2 + 400, height: radius * 2 + 400, xScale: xScaleR, maxDist: maxDist };
}

// ============ COLOR SCHEMES ============
const PALETTE = ["#7c3aed","#e03e5c","#0891b2","#16a34a","#d97706","#be185d","#7c3aed","#4f46e5","#059669","#dc2626","#9333ea","#0d9488","#ca8a04"];

// P0 FIX: pre-computed trace data (set once in init, reused everywhere)
let traceLCA = null;
let traceLCADescIds = new Set();
let validTracedGenomes = [];

function computeRadialRadius() { return Math.max(400, allLeaves.length * 2.5); } // P2: single source of truth

function getColorScale(field) {
  if (field === 'trace') {
    // P0 FIX: use pre-computed sets, no per-leaf LCA calls
    return (d) => {
      if (validTracedGenomes.includes(d.name)) return '#e03e5c';
      if (traceLCADescIds.has(d.id)) return 'rgba(224,62,92,0.55)';
      return '#b0b0c0';
    };
  }
  let vals = [...new Set(allLeaves.map(l => (META_DATA[l.name]||{})[field]||'NA'))].sort();
  if (vals.length <= 13) {
    let map = {}; vals.forEach((v,i) => map[v] = PALETTE[i % PALETTE.length]);
    return (d) => map[(META_DATA[d.name]||{})[field]||'NA'] || '#b0b0c0';
  }
  let colFn = d3.interpolateTurbo;
  let map = {}; vals.forEach((v,i) => map[v] = colFn(i / vals.length));
  return (d) => map[(META_DATA[d.name]||{})[field]||'NA'] || '#b0b0c0';
}
let currentColorFn;
let currentColorField = TYPING_COLS[0] || 'trace';

// ============ RENDER ============
function render() {
  d3.select('#tree-container').selectAll('*').remove();
  let container = document.getElementById('tree-container');
  let W = container.clientWidth, H = container.clientHeight;
  let dims = currentLayout === 'rect' ? layoutRect(treeRoot) : layoutRadial(treeRoot);
  lastDims = dims;

  svg = d3.select('#tree-container').append('svg').attr('width', W).attr('height', H);
  g = svg.append('g');
  zoom = d3.zoom().scaleExtent([0.05, 12]).on('zoom', (e) => {
    g.attr('transform', e.transform);
    document.getElementById('info-zoom').textContent = 'Zoom: ' + e.transform.k.toFixed(2) + 'x';
    updateScaleBar(e.transform.k);
  });
  svg.call(zoom);
  currentColorFn = getColorScale(currentColorField);

  // Use pre-computed trace data (P0 fix: no repeated LCA calculation)
  let lca = traceLCA;
  let lcaDescIds = traceLCADescIds;

  // Clade highlight
  if (isNodeTrace && lca) {
    let descLeaves = getDescendants(lca).filter(n => n.isLeaf);
    if (descLeaves.length > 0 && currentLayout === 'rect') {
      let minY = Math.min(...descLeaves.map(l => l.y));
      let maxY = Math.max(...descLeaves.map(l => l.y));
      let minX = lca.x, maxX = Math.max(...descLeaves.map(l => l.x)) + 200;
      g.append('rect').attr('x', minX-4).attr('y', minY-tipHeight/2-2)
        .attr('width', maxX-minX+8).attr('height', maxY-minY+tipHeight+4)
        .attr('rx', 6).attr('fill','rgba(224,62,92,0.06)')
        .attr('stroke','rgba(224,62,92,0.22)').attr('stroke-width',1.5).attr('stroke-dasharray','6,3');
      g.append('text').attr('x', minX-8).attr('y', (minY+maxY)/2)
        .attr('text-anchor','end').attr('font-size','12px')
        .attr('font-family',"'JetBrains Mono',monospace")
        .attr('fill','rgba(224,62,92,0.7)').attr('font-weight','600')
        .text('Traced Clade (' + descLeaves.length + ' tips)');
    }
    if (descLeaves.length > 0 && currentLayout === 'radial') {
      let angles = descLeaves.map(l=>l.angle);
      let minA=Math.min(...angles)-0.01, maxA=Math.max(...angles)+0.01;
      let maxDist=getMaxDepth(treeRoot), radius=computeRadialRadius();
      let innerR=(lca.rootDist/Math.max(maxDist,0.001))*radius, outerR=radius*1.08;
      let x1=outerR*Math.cos(minA-Math.PI/2),y1=outerR*Math.sin(minA-Math.PI/2);
      let x2=outerR*Math.cos(maxA-Math.PI/2),y2=outerR*Math.sin(maxA-Math.PI/2);
      let ix1=innerR*Math.cos(minA-Math.PI/2),iy1=innerR*Math.sin(minA-Math.PI/2);
      let ix2=innerR*Math.cos(maxA-Math.PI/2),iy2=innerR*Math.sin(maxA-Math.PI/2);
      let la=(maxA-minA)>Math.PI?1:0;
      g.append('path').attr('d','M '+ix1+' '+iy1+' A '+innerR+' '+innerR+' 0 '+la+' 1 '+ix2+' '+iy2+' L '+x2+' '+y2+' A '+outerR+' '+outerR+' 0 '+la+' 0 '+x1+' '+y1+' Z')
        .attr('fill','rgba(224,62,92,0.06)').attr('stroke','rgba(224,62,92,0.18)').attr('stroke-width',1).attr('stroke-dasharray','6,3');
    }
  }

  // Branches
  let branchGroup = g.append('g').attr('class','branches');
  let branchNormal = 'rgba(0,0,0,0.35)';
  let branchTrace = 'rgba(224,62,92,0.75)';

  function drawBranchesRect(node) {
    if (node.children.length > 0) {
      let childYs = node.children.map(c=>c.y);
      let t = lcaDescIds.has(node.id);
      branchGroup.append('line').attr('x1',node.x).attr('y1',Math.min(...childYs)).attr('x2',node.x).attr('y2',Math.max(...childYs))
        .attr('stroke', t?branchTrace:branchNormal).attr('stroke-width', t?2.5:1.2);
    }
    node.children.forEach(child => {
      let t = lcaDescIds.has(child.id) && lcaDescIds.has(node.id);
      branchGroup.append('line').attr('x1',node.x).attr('y1',child.y).attr('x2',child.x).attr('y2',child.y)
        .attr('stroke', t?branchTrace:branchNormal).attr('stroke-width', t?2.5:1.2);
      drawBranchesRect(child);
    });
  }
  function drawBranchesRadial(node) {
    node.children.forEach(child => {
      let t = lcaDescIds.has(child.id)&&lcaDescIds.has(node.id);
      let col=t?branchTrace:branchNormal, sw=t?2.2:1.0;
      let maxD=getMaxDepth(treeRoot),radius=computeRadialRadius();
      let r1=(node.rootDist/Math.max(maxD,0.001))*radius;
      let angle=child.angle-Math.PI/2, pAngle=node.angle-Math.PI/2;
      let px1=r1*Math.cos(angle),py1=r1*Math.sin(angle);
      let arcPx=r1*Math.cos(pAngle),arcPy=r1*Math.sin(pAngle);
      if(Math.abs(child.angle-node.angle)>0.001){
        let la=Math.abs(child.angle-node.angle)>Math.PI?1:0,sweep=child.angle>node.angle?1:0;
        branchGroup.append('path').attr('d',`M ${arcPx} ${arcPy} A ${r1} ${r1} 0 ${la} ${sweep} ${px1} ${py1}`)
          .attr('fill','none').attr('stroke',col).attr('stroke-width',sw);
      }
      branchGroup.append('line').attr('x1',px1).attr('y1',py1).attr('x2',child.x).attr('y2',child.y)
        .attr('stroke',col).attr('stroke-width',sw);
      drawBranchesRadial(child);
    });
  }
  if (currentLayout==='rect') drawBranchesRect(treeRoot); else drawBranchesRadial(treeRoot);

  // Tips
  let tipGroup = g.append('g').attr('class','tips');
  allLeaves.forEach(leaf => {
    let isTraced = validTracedGenomes.includes(leaf.name);
    let color = currentColorFn(leaf);
    let lg = tipGroup.append('g').attr('transform',`translate(${leaf.x},${leaf.y})`)
      .style('cursor','pointer').datum(leaf)
      .on('mouseover',(e,d)=>showTooltip(e,d)).on('mouseout',hideTooltip).on('click',(e,d)=>showDetail(d,e));
    if (isTraced) {
      lg.append('circle').attr('r',12).attr('fill','rgba(224,62,92,0.15)').attr('class','trace-pulse');
      lg.append('circle').attr('r',6).attr('fill','#e03e5c').attr('stroke','#fff').attr('stroke-width',2);
    } else {
      lg.append('circle').attr('r',4.5).attr('fill',color).attr('opacity',0.9);
    }
    if (showLabels || isTraced) {
      let txt = lg.append('text')
        .attr('font-size', isTraced?'13px':'10.5px')
        .attr('font-family',"'JetBrains Mono',monospace")
        .attr('fill', isTraced?'#e03e5c':'rgba(0,0,0,0.55)')
        .attr('font-weight', isTraced?'700':'400')
        .text(leaf.name)
        .attr('class', showLabels?'':(isTraced?'':'hidden-label'));
      if (currentLayout==='radial' && leaf.angle!==undefined) {
        let deg=leaf.angle*180/Math.PI, flip=deg>90&&deg<270;
        txt.attr('transform',`rotate(${flip?deg-180:deg})`).attr('x',flip?-(isTraced?12:8):(isTraced?12:8))
          .attr('y',3.5).attr('text-anchor',flip?'end':'start');
      } else { txt.attr('x',isTraced?12:8).attr('y',3.5); }
    }
  });
  if(!showLabels) g.selectAll('.hidden-label').remove();

  // Initial transform
  if(currentLayout==='rect'){
    let s=Math.min(W/(dims.width+100),H/(dims.height+100),1);
    svg.call(zoom.transform,d3.zoomIdentity.translate(40,20).scale(s));
  } else {
    let s=Math.min(W,H)/(dims.width+200);
    svg.call(zoom.transform,d3.zoomIdentity.translate(W/2,H/2).scale(s));
  }
  document.getElementById('info-tip-count').textContent=allLeaves.length+' tips';
  updateLegend();
  updateScaleBar();
}

// ============ SCALE BAR (fixed overlay, zoom-aware) ============
function updateScaleBar(zoomK) {
  let overlay = document.getElementById('scale-bar-overlay');
  if (!lastDims || !lastDims.xScale || lastDims.maxDist <= 0) {
    overlay.style.display = 'none'; return;
  }
  // Get current zoom scale
  if (zoomK === undefined) {
    try { let t = d3.zoomTransform(svg.node()); zoomK = t.k; } catch(e) { zoomK = 1; }
  }
  let xScale = lastDims.xScale * zoomK; // pixels-per-unit at current zoom
  let maxDist = lastDims.maxDist;

  // Target bar width: 80-150 screen pixels
  let targetPx = 100;
  let targetVal = targetPx / xScale; // distance value for ~100px
  // Snap to a nice round number
  let mag = Math.pow(10, Math.floor(Math.log10(targetVal)));
  let candidates = [1,2,2.5,5,10].map(m => m * mag);
  let scaleVal = candidates.reduce((best, c) =>
    Math.abs(c * xScale - targetPx) < Math.abs(best * xScale - targetPx) ? c : best
  );
  let barPx = Math.round(scaleVal * xScale);
  if (barPx < 30) barPx = 30;
  if (barPx > 250) barPx = 250;

  let label = scaleVal >= 0.001 ? parseFloat(scaleVal.toPrecision(3)).toString() : scaleVal.toExponential(1);

  overlay.innerHTML = '<div class="sb-bar"><div class="sb-line" style="width:'+barPx+'px"></div></div>'
    + '<div class="sb-label">'+label+' '+DIST_UNIT+'</div>';
  overlay.style.display = 'block';
}

// ============ NODE → SCREEN COORDINATES ============
function leafToScreen(leaf) {
  // Convert SVG tree-space coords to screen (viewport) coords
  let svgEl = svg.node();
  let pt = svgEl.createSVGPoint();
  pt.x = leaf.x; pt.y = leaf.y;
  // Apply the current zoom/pan transform
  let ctm = g.node().getCTM();
  if (ctm) { let sp = pt.matrixTransform(ctm); return { x: sp.x, y: sp.y }; }
  return { x: 0, y: 0 };
}

// ============ TOOLTIP (hover) ============
function showTooltip(event, leaf) {
  let meta = META_DATA[leaf.name]||{};
  let tt = document.getElementById('tooltip');
  let html = '<div class="tt-name">'+leaf.name+'</div>';
  if(validTracedGenomes.includes(leaf.name)) html+='<div style="color:#e03e5c;font-size:11px;font-weight:700;margin-bottom:6px;">&#9733; TRACED GENOME</div>';
  TYPING_COLS.forEach(col => {
    html += '<div class="tt-row"><span class="tt-key">'+col+'</span><span class="tt-val">'+(meta[col]||'NA')+'</span></div>';
  });
  let keyGenes=['blaKPC','blaNDM','blaVIM','blaIMP','blaCTX-M-65','iucA','rmpA','rmpA2','ybtS'];
  let present=keyGenes.filter(g=>meta[g]==='1');
  if(present.length>0){
    html+='<div class="tt-section">Key Genes Detected</div>';
    html+='<div class="tt-genes">'+present.map(g=>'<span class="tt-gene present">'+g+'</span>').join('')+'</div>';
  }
  tt.innerHTML=html; tt.classList.add('visible');

  let ttW = tt.offsetWidth || 400, ttH = tt.offsetHeight || 200;
  let winW = window.innerWidth, winH = window.innerHeight;
  let sidebarW = document.getElementById('sidebar').offsetWidth;

  // Get node screen position → place tooltip to the right of the label
  let scr = leafToScreen(leaf);
  let gap = 160; // approximate label width in pixels at current zoom
  let x = scr.x + gap;
  let y = scr.y - ttH / 2;

  // If it overflows right, flip to left of node
  if (x + ttW > winW - 10) { x = scr.x - ttW - 20; }
  // If it overflows left into sidebar, fallback to cursor position
  if (x < sidebarW + 5) { x = event.clientX + 15; }
  // Clamp vertical
  if (y < 10) y = 10;
  if (y + ttH > winH - 10) y = winH - ttH - 10;

  tt.style.left = x + 'px'; tt.style.top = y + 'px';
}
function hideTooltip(){ document.getElementById('tooltip').classList.remove('visible'); }

// ============ DETAIL PANEL ============
let _detailOpacity = 0.95;
function setDetailOpacity(val) {
  _detailOpacity = val / 100;
  document.getElementById('detail-panel').style.opacity = _detailOpacity;
}

function showDetail(leaf, clickEvent) {
  let meta=META_DATA[leaf.name]||{};
  let dp=document.getElementById('detail-panel');
  let html='';
  if(validTracedGenomes.includes(leaf.name)) html+='<div class="dp-trace-badge">&#9733; TRACED GENOME</div>';
  html+='<div class="dp-name">'+leaf.name+'</div>';
  html+='<div class="dp-section">Typing Information</div>';
  TYPING_COLS.forEach(col=>{
    html+='<div class="dp-row"><span class="dp-key">'+col+'</span><span class="dp-val">'+(meta[col]||'NA')+'</span></div>';
  });
  if(VIR_GENES.length>0){
    html+='<div class="dp-section">Virulence Genes</div><div class="gene-grid">';
    VIR_GENES.forEach(g=>{ let p=meta[g]==='1'; html+='<div class="gene-chip '+(p?'present':'absent')+'">'+g+' '+(p?'&#10003;':'&#10007;')+'</div>'; });
    html+='</div>';
  }
  if(RES_GENES.length>0){
    html+='<div class="dp-section">Resistance Genes</div><div class="gene-grid">';
    RES_GENES.forEach(g=>{ let p=meta[g]==='1'; html+='<div class="gene-chip '+(p?'present':'absent')+'">'+g+' '+(p?'&#10003;':'&#10007;')+'</div>'; });
    html+='</div>';
  }
  document.getElementById('detail-content').innerHTML=html;
  dp.style.opacity = _detailOpacity;
  document.getElementById('dp-opacity-slider').value = Math.round(_detailOpacity * 100);
  dp.classList.add('visible');

  // Position: anchor beside the node
  let dpW = 370, winW = window.innerWidth, winH = window.innerHeight;
  let sidebarW = document.getElementById('sidebar').offsetWidth;
  let scr = leafToScreen(leaf);
  let dpH = dp.offsetHeight || 500;

  if (currentLayout === 'rect') {
    // Rectangular: place to the right of the tip label
    let labelGap = 180;
    let x = scr.x + labelGap;
    // If overflows right edge, try left of node
    if (x + dpW > winW - 10) { x = scr.x - dpW - 30; }
    // If overflows into sidebar, clamp to right edge
    if (x < sidebarW + 5) { x = winW - dpW - 12; }
    dp.style.left = x + 'px';
    dp.style.right = 'auto';
    // Vertical: align center with node, clamp
    let top = scr.y - dpH / 3;
    top = Math.max(12, Math.min(top, winH - dpH - 12));
    dp.style.top = top + 'px';
  } else {
    // Radial / fallback: use click position, opposite side
    let mx = clickEvent ? clickEvent.clientX : scr.x;
    let treeMid = sidebarW + (winW - sidebarW) / 2;
    if (mx > treeMid) {
      dp.style.left = (sidebarW + 12) + 'px'; dp.style.right = 'auto';
    } else {
      dp.style.left = (winW - dpW - 12) + 'px'; dp.style.right = 'auto';
    }
    let my = clickEvent ? clickEvent.clientY : scr.y;
    let top = Math.max(12, Math.min(my - dpH/3, winH - dpH - 12));
    dp.style.top = top + 'px';
  }
}
function closeDetail(){ document.getElementById('detail-panel').classList.remove('visible'); }

// ============ DRAG DETAIL PANEL ============
(function initDrag(){
  let dp, handle, dragging=false, ox, oy, sx, sy;
  function setup() {
    dp = document.getElementById('detail-panel');
    handle = document.getElementById('dp-drag-handle');
    if (!handle) return;
    handle.addEventListener('mousedown', onDown);
    handle.addEventListener('touchstart', onTouchDown, {passive:false});
  }
  function onDown(e) {
    // Ignore clicks on buttons / inputs inside toolbar
    if (e.target.tagName==='BUTTON'||e.target.tagName==='INPUT'||e.target.tagName==='LABEL') return;
    e.preventDefault();
    dragging=true;
    ox=e.clientX; oy=e.clientY;
    sx=dp.offsetLeft; sy=dp.offsetTop;
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onUp);
    dp.style.transition='none';
  }
  function onMove(e) {
    if(!dragging) return;
    let dx=e.clientX-ox, dy=e.clientY-oy;
    let nx=sx+dx, ny=sy+dy;
    // Clamp within viewport
    let winW=window.innerWidth, winH=window.innerHeight, dpW=dp.offsetWidth, dpH=dp.offsetHeight;
    nx=Math.max(0, Math.min(nx, winW-dpW));
    ny=Math.max(0, Math.min(ny, winH-dpH));
    dp.style.left=nx+'px'; dp.style.top=ny+'px'; dp.style.right='auto';
  }
  function onUp() {
    dragging=false;
    document.removeEventListener('mousemove', onMove);
    document.removeEventListener('mouseup', onUp);
    dp.style.transition='opacity 0.15s';
  }
  // Touch support
  function onTouchDown(e) {
    if (e.target.tagName==='BUTTON'||e.target.tagName==='INPUT'||e.target.tagName==='LABEL') return;
    e.preventDefault();
    let t=e.touches[0]; dragging=true;
    ox=t.clientX; oy=t.clientY;
    sx=dp.offsetLeft; sy=dp.offsetTop;
    document.addEventListener('touchmove', onTouchMove, {passive:false});
    document.addEventListener('touchend', onTouchEnd);
    dp.style.transition='none';
  }
  function onTouchMove(e) {
    if(!dragging) return; e.preventDefault();
    let t=e.touches[0], dx=t.clientX-ox, dy=t.clientY-oy;
    let nx=sx+dx, ny=sy+dy;
    let winW=window.innerWidth, winH=window.innerHeight, dpW=dp.offsetWidth, dpH=dp.offsetHeight;
    nx=Math.max(0, Math.min(nx, winW-dpW));
    ny=Math.max(0, Math.min(ny, winH-dpH));
    dp.style.left=nx+'px'; dp.style.top=ny+'px'; dp.style.right='auto';
  }
  function onTouchEnd() {
    dragging=false;
    document.removeEventListener('touchmove', onTouchMove);
    document.removeEventListener('touchend', onTouchEnd);
    dp.style.transition='opacity 0.15s';
  }
  if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',setup); else setup();
})();

// ============ LEGEND ============
function updateLegend() {
  let legend=document.getElementById('legend');
  if(currentColorField==='trace'){
    legend.innerHTML='<div class="legend-item"><div class="legend-dot" style="background:#e03e5c"></div>Traced genome</div>'
      +'<div class="legend-item"><div class="legend-dot" style="background:rgba(224,62,92,0.55)"></div>Traced clade</div>'
      +'<div class="legend-item"><div class="legend-dot" style="background:#b0b0c0"></div>Other</div>';
    return;
  }
  let vals={};
  allLeaves.forEach(l=>{let v=(META_DATA[l.name]||{})[currentColorField]||'NA';if(!vals[v])vals[v]={count:0,color:currentColorFn(l)};vals[v].count++;});
  let sorted=Object.entries(vals).sort((a,b)=>b[1].count-a[1].count);
  let mx=20;
  legend.innerHTML=sorted.slice(0,mx).map(([v,{count,color}])=>
    '<div class="legend-item" onclick="highlightValue(\''+v.replace(/'/g,"\\'")+'\')">'
    +'<div class="legend-dot" style="background:'+color+'"></div>'
    +'<span>'+v+'</span>'
    +'<span style="margin-left:auto;color:var(--text-muted);font-size:11px">'+count+'</span></div>'
  ).join('');
  if(sorted.length>mx) legend.innerHTML+='<div style="font-size:10px;color:var(--text-muted);padding:4px 0;">+' +(sorted.length-mx)+' more</div>';
}
function highlightValue(val) {
  g.selectAll('.tips g').each(function(){
    let leaf=d3.select(this).datum();if(!leaf)return;
    let v=(META_DATA[leaf.name]||{})[currentColorField]||'NA';
    d3.select(this).selectAll('circle').attr('opacity',v===val?1:0.1);
  });
  setTimeout(()=>{ g.selectAll('.tips g').selectAll('circle').attr('opacity',0.9); },3000);
}

// ============ CONTROLS ============
function setLayout(layout) {
  currentLayout=layout;
  document.getElementById('btn-rect').classList.toggle('active',layout==='rect');
  document.getElementById('btn-radial').classList.toggle('active',layout==='radial');
  render();
}
function focusOnTrace() {
  let targets=allLeaves.filter(l=>validTracedGenomes.includes(l.name));
  if(!targets.length)return;
  let W=document.getElementById('tree-container').clientWidth,H=document.getElementById('tree-container').clientHeight;
  if(targets.length===1){let t=targets[0];svg.transition().duration(750).call(zoom.transform,d3.zoomIdentity.translate(W/2-t.x*3,H/2-t.y*3).scale(3));}
  else{let xs=targets.map(t=>t.x),ys=targets.map(t=>t.y);let cx=(Math.min(...xs)+Math.max(...xs))/2,cy=(Math.min(...ys)+Math.max(...ys))/2;let dx=Math.max(...xs)-Math.min(...xs)+200,dy=Math.max(...ys)-Math.min(...ys)+200;let s=Math.min(W/dx,H/dy,4);svg.transition().duration(750).call(zoom.transform,d3.zoomIdentity.translate(W/2-cx*s,H/2-cy*s).scale(s));}
}
function resetView() {
  let W=document.getElementById('tree-container').clientWidth,H=document.getElementById('tree-container').clientHeight;
  if(currentLayout==='rect'){let d={width:allLeaves.length>0?Math.max(...allLeaves.map(l=>l.x))+300:1000,height:allLeaves.length*tipHeight+40};let s=Math.min(W/(d.width+100),H/(d.height+100),1);svg.transition().duration(500).call(zoom.transform,d3.zoomIdentity.translate(40,20).scale(s));}
  else{let r=computeRadialRadius();let s=Math.min(W,H)/(r*2+400);svg.transition().duration(500).call(zoom.transform,d3.zoomIdentity.translate(W/2,H/2).scale(s));}
}
function toggleLabels(){showLabels=!showLabels;render();}
function zoomIn(){svg.transition().duration(300).call(zoom.scaleBy,1.4);}
function zoomOut(){svg.transition().duration(300).call(zoom.scaleBy,0.7);}

// ============ SEARCH ============
let searchHighlighted = []; // track highlighted genomes on tree

document.getElementById('search-box').addEventListener('input', function() {
  let raw = this.value.trim();
  let results = document.getElementById('search-results');
  if (!raw) { results.innerHTML = ''; clearSearchHighlight(); return; }

  // Split by ; or whitespace to support multiple queries
  let queries = raw.split(/[;\s]+/).map(s => s.toLowerCase()).filter(Boolean);
  if (!queries.length) { results.innerHTML = ''; return; }

  // Collect unique matches across all queries
  let matchSet = new Map();
  queries.forEach(q => {
    allLeaves.filter(l => l.name.toLowerCase().includes(q)).forEach(l => {
      if (!matchSet.has(l.name)) matchSet.set(l.name, { leaf: l, matchedQ: q });
    });
  });
  let matches = [...matchSet.values()].slice(0, 15);

  // Build result list
  let html = '';
  if (matches.length > 1) {
    html += '<div class="search-result-item search-action" onclick="searchSelectAll()" '
         +  'style="color:var(--accent);font-weight:600;border-bottom:1px solid var(--border-light);padding-bottom:6px;margin-bottom:4px">'
         +  '&#10140; Locate all ' + matches.length + ' matches on tree</div>';
  }
  html += matches.map(({leaf: l, matchedQ: q}) => {
    let hi = l.name.replace(new RegExp('(' + q.replace(/[-\/\\^$*+?.()|[\]{}]/g, '\\$&') + ')', 'gi'), '<b style="color:var(--accent)">$1</b>');
    let star = validTracedGenomes.includes(l.name) ? ' <span style="color:var(--trace-color)">&#9733;</span>' : '';
    return '<div class="search-result-item" onclick="searchSelect(\'' + l.name + '\')">' + hi + star + '</div>';
  }).join('');

  results.innerHTML = html;

  // Store current search matches for "locate all"
  searchHighlighted = matches.map(m => m.leaf);
});

function searchSelect(name) {
  let leaf = allLeaves.find(l => l.name === name);
  if (!leaf) return;
  document.getElementById('search-results').innerHTML = '';
  document.getElementById('search-box').value = name;
  let W = document.getElementById('tree-container').clientWidth, H = document.getElementById('tree-container').clientHeight;
  svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity.translate(W/2 - leaf.x*4, H/2 - leaf.y*4).scale(4));
  showDetail(leaf, null);
  highlightSearchOnTree([leaf]);
}

function searchSelectAll() {
  if (!searchHighlighted.length) return;
  document.getElementById('search-results').innerHTML = '';

  let targets = searchHighlighted;
  let W = document.getElementById('tree-container').clientWidth, H = document.getElementById('tree-container').clientHeight;
  let xs = targets.map(t => t.x), ys = targets.map(t => t.y);
  let cx = (Math.min(...xs) + Math.max(...xs)) / 2;
  let cy = (Math.min(...ys) + Math.max(...ys)) / 2;
  let dx = Math.max(...xs) - Math.min(...xs) + 300;
  let dy = Math.max(...ys) - Math.min(...ys) + 300;
  let scale = Math.min(W / dx, H / dy, 4);
  svg.transition().duration(750).call(zoom.transform,
    d3.zoomIdentity.translate(W/2 - cx*scale, H/2 - cy*scale).scale(scale));

  highlightSearchOnTree(targets);
}

function highlightSearchOnTree(targets) {
  let nameSet = new Set(targets.map(t => t.name));
  // Add temporary rings to matched genomes
  g.selectAll('.search-ring').remove();
  targets.forEach(leaf => {
    g.append('circle')
      .attr('class', 'search-ring')
      .attr('cx', leaf.x).attr('cy', leaf.y)
      .attr('r', 12).attr('fill', 'none')
      .attr('stroke', 'var(--accent)').attr('stroke-width', 2.5)
      .attr('stroke-dasharray', '4,2')
      .attr('opacity', 1)
      .transition().duration(4000).attr('opacity', 0).remove();
  });
  // Dim other tips temporarily
  g.selectAll('.tips g').each(function() {
    let d = d3.select(this).datum();
    if (!d) return;
    d3.select(this).selectAll('circle').attr('opacity', nameSet.has(d.name) ? 1 : 0.15);
  });
  setTimeout(() => {
    g.selectAll('.tips g').selectAll('circle').attr('opacity', 0.9);
  }, 4000);
}

function clearSearchHighlight() {
  searchHighlighted = [];
  g.selectAll('.search-ring').remove();
  g.selectAll('.tips g').selectAll('circle').attr('opacity', 0.9);
}

// ============ COLOR SELECT ============
document.getElementById('color-select').addEventListener('change', function(){ currentColorField=this.value; render(); });

// ============ INIT ============
function init() {
  treeRoot = parseNewick(NEWICK_STR);
  assignIds(treeRoot); countLeaves(treeRoot); computeDepth(treeRoot);

  // P1 FIX: validate traced genomes exist in tree
  let leafNames = new Set(allLeaves.map(l => l.name));
  let missing = tracedGenomes.filter(g => !leafNames.has(g));
  validTracedGenomes = tracedGenomes.filter(g => leafNames.has(g));
  if (missing.length > 0) {
    console.warn('[mTracePatho] Traced genomes NOT found in tree:', missing);
  }

  // P0 FIX: pre-compute LCA once
  if (validTracedGenomes.length > 0) {
    traceLCA = findLCA(treeRoot, validTracedGenomes);
    traceLCADescIds = traceLCA ? new Set(getDescendants(traceLCA).map(n=>n.id)) : new Set();
  }
  isNodeTrace = validTracedGenomes.length > 1;

  let traceInfo=document.getElementById('trace-info');
  let traceHtml='<div class="label">Traced To'+(isNodeTrace?' (Node - Multiple Genomes)':' (Single Genome)')+'</div>'
    +'<div class="value">'+tracedGenomes.join('<br>')+'</div>';
  if (missing.length > 0) {
    traceHtml += '<div style="color:#d97706;font-size:11px;margin-top:6px;font-weight:600;">&#9888; ' + missing.length + ' genome(s) not found in tree: ' + missing.join(', ') + '</div>';
  }
  traceHtml += '<span class="type-badge '+(isNodeTrace?'node':'single')+'">'+(isNodeTrace?'Node Trace ('+validTracedGenomes.length+' genomes)':'Single Genome')+'</span>';
  traceInfo.innerHTML = traceHtml;

  let maxDist=getMaxDepth(treeRoot);
  document.getElementById('stats').innerHTML=
    '<div class="stat-row"><span class="stat-label">Species</span><span class="stat-value" style="font-style:italic">'+TAX_NAME+'</span></div>'
    +'<div class="stat-row"><span class="stat-label">Total Tips</span><span class="stat-value">'+allLeaves.length+'</span></div>'
    +'<div class="stat-row"><span class="stat-label">Internal Nodes</span><span class="stat-value">'+(allNodes.length-allLeaves.length)+'</span></div>'
    +'<div class="stat-row"><span class="stat-label">Max Root Dist</span><span class="stat-value">'+maxDist.toFixed(4)+'</span></div>'
    +'<div class="stat-row"><span class="stat-label">Traced Genomes</span><span class="stat-value">'+validTracedGenomes.length+(missing.length?' <span style="color:#d97706">('+missing.length+' missing)</span>':'')+'</span></div>'
    +'<div class="stat-row"><span class="stat-label">Trace Type</span><span class="stat-value">'+(isNodeTrace?'Node':'Single')+'</span></div>';

  render();
  document.getElementById('loading').style.display='none';

  // Auto-focus on traced genomes at startup (slightly zoomed out for context)
  if (validTracedGenomes.length > 0) {
    setTimeout(function(){
      let targets = allLeaves.filter(l => validTracedGenomes.includes(l.name));
      if (!targets.length) return;
      let W = document.getElementById('tree-container').clientWidth;
      let H = document.getElementById('tree-container').clientHeight;
      let ZOOM_REDUCE = 0.8; // 20% smaller than tight focus
      if (targets.length === 1) {
        let t = targets[0];
        let s = 3 * ZOOM_REDUCE;
        svg.transition().duration(900).call(zoom.transform,
          d3.zoomIdentity.translate(W/2 - t.x*s, H/2 - t.y*s).scale(s));
      } else {
        let xs = targets.map(t=>t.x), ys = targets.map(t=>t.y);
        let cx = (Math.min(...xs)+Math.max(...xs))/2;
        let cy = (Math.min(...ys)+Math.max(...ys))/2;
        let dx = Math.max(...xs)-Math.min(...xs)+300;
        let dy = Math.max(...ys)-Math.min(...ys)+300;
        let s = Math.min(W/dx, H/dy, 4) * ZOOM_REDUCE;
        svg.transition().duration(900).call(zoom.transform,
          d3.zoomIdentity.translate(W/2 - cx*s, H/2 - cy*s).scale(s));
      }
    }, 200);
  }
}
window.addEventListener('load',()=>setTimeout(init,100));

// P2 FIX: debounced resize
let _resizeTimer;
window.addEventListener('resize',()=>{
  clearTimeout(_resizeTimer);
  _resizeTimer = setTimeout(()=>{ if(treeRoot) render(); }, 200);
});
</script>
</body>
</html>'''
    return html



def generate_viewer_from_data(
        tree_path,
        meta_path,
        trace_genomes_str,
        taxname,
        output_html,
        dist_unit='Jaccard distance',
        log_fn=None):
    """Generate the HTML viewer directly from Python code.

    This is the programmatic interface expected by ``mTracePathoRun.py``.
    It preserves the original command-line ``main()`` interface while allowing
    the tracing pipeline to generate one viewer per sample without launching a
    subprocess.

    Returns
    -------
    bool
        ``True`` when both HTML and companion JSON files are generated;
        ``False`` when validation or generation fails.
    """
    def log(message):
        if log_fn is not None:
            log_fn(str(message))
        else:
            print(str(message))

    try:
        for path, label in [
            (tree_path, 'Tree'),
            (meta_path, 'Metadata'),
        ]:
            if not path or not os.path.isfile(path):
                log(f"[WARN] {label} file not found: {path}")
                return False

        trace_genomes_str = str(trace_genomes_str or '').strip()
        if not trace_genomes_str:
            log('[WARN] Empty tracing result; viewer was not generated')
            return False

        output_html = os.path.abspath(output_html)
        output_dir = os.path.dirname(output_html)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        tree_nwk = read_tree(tree_path)
        if not tree_nwk:
            log(f"[WARN] Tree file is empty: {tree_path}")
            return False

        meta, cols = read_meta(meta_path)
        cols = cols or []
        typing_cols, gene_cols = classify_columns(cols)
        vir_genes, res_genes = classify_genes(gene_cols)

        # Keep metadata records represented in the Newick tree.  Use a general
        # token parser as a fallback so non-GCF identifiers also remain usable.
        tree_genomes = set(re.findall(r'GCF_\d+\.\d+', tree_nwk))
        if tree_genomes:
            filtered_meta = {k: v for k, v in meta.items() if k in tree_genomes}
        else:
            filtered_meta = dict(meta)

        meta_json = json.dumps(filtered_meta, separators=(',', ':'))
        html = generate_html(
            tree_nwk,
            meta_json,
            trace_genomes_str,
            taxname,
            typing_cols,
            vir_genes,
            res_genes,
            dist_unit,
        )

        with open(output_html, 'w', encoding='utf-8') as handle:
            handle.write(html)

        traced_genomes = [
            genome.strip() for genome in trace_genomes_str.split(';')
            if genome.strip()
        ]
        json_output = os.path.splitext(output_html)[0] + '.data.json'
        interface_data = {
            '_comment': 'mTracePatho Viewer interface data — for frontend web integration',
            '_version': '1.0',
            'params': {
                'taxname': taxname,
                'dist_unit': dist_unit,
            },
            'data': {
                'NEWICK_STR': tree_nwk,
                'META_DATA': filtered_meta,
                'TRACE_RAW': trace_genomes_str,
                'TAX_NAME': format_taxname(taxname),
                'TYPING_COLS': typing_cols,
                'VIR_GENES': vir_genes,
                'RES_GENES': res_genes,
                'DIST_UNIT': dist_unit,
            },
            'stats': {
                'tree_tips': len(tree_genomes) if tree_genomes else len(filtered_meta),
                'meta_matched': len(filtered_meta),
                'traced_genomes': traced_genomes,
                'trace_type': 'node' if len(traced_genomes) > 1 else 'single',
                'typing_col_count': len(typing_cols),
                'vir_gene_count': len(vir_genes),
                'res_gene_count': len(res_genes),
            },
        }
        with open(json_output, 'w', encoding='utf-8') as handle:
            json.dump(interface_data, handle, ensure_ascii=False, indent=2)

        html_kb = os.path.getsize(output_html) / 1024
        json_kb = os.path.getsize(json_output) / 1024
        log(f"Viewer HTML generated: {output_html} ({html_kb:.0f} KB)")
        log(f"Viewer JSON generated: {json_output} ({json_kb:.0f} KB)")
        return True

    except Exception as exc:
        log(f"[WARN] Viewer generation failed: {type(exc).__name__}: {exc}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='mTracePatho Phylogenetic Tracing Viewer Generator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Example:
  python generate_trace_viewer.py \\
      --tree Jaccard_fastme_input.nwk \\
      --meta Klebsiella_pneumoniae-gene.tsv \\
      --trace trace_result.tsv \\
      --taxname Klebsiella_pneumoniae \\
      --dist-unit "Jaccard distance" \\
      --output mTracePatho_viewer.html
        '''
    )
    parser.add_argument('--tree', required=True, help='Newick tree file (.nwk)')
    parser.add_argument('--meta', required=True, help='Genome metadata TSV file')
    parser.add_argument('--trace', required=True, help='Trace result TSV (single line, genome IDs separated by ;)')
    parser.add_argument('--taxname', required=True, help='Species name, e.g. Klebsiella_pneumoniae')
    parser.add_argument('--output', default='mTracePatho_viewer.html', help='Output HTML file (default: mTracePatho_viewer.html)')
    parser.add_argument('--dist-unit', default='distance',
        help='Scale bar unit label, e.g. "substitutions/site", "Jaccard distance"')

    args = parser.parse_args()

    # Validate inputs
    for f, label in [(args.tree, 'Tree'), (args.meta, 'Meta'), (args.trace, 'Trace')]:
        if not os.path.isfile(f):
            print(f"[ERROR] {label} file not found: {f}", file=sys.stderr)
            sys.exit(1)

    print(f"[INFO] Reading tree: {args.tree}")
    tree_nwk = read_tree(args.tree)
    print(f"  Tree string length: {len(tree_nwk)} chars")

    print(f"[INFO] Reading metadata: {args.meta}")
    meta, cols = read_meta(args.meta)
    typing_cols, gene_cols = classify_columns(cols)
    vir_genes, res_genes = classify_genes(gene_cols)
    print(f"  Metadata entries: {len(meta)}")
    print(f"  Typing columns: {typing_cols}")
    print(f"  Gene columns: {len(gene_cols)} (vir: {len(vir_genes)}, res: {len(res_genes)})")

    print(f"[INFO] Reading trace result: {args.trace}")
    trace_result = read_trace(args.trace)
    traced_genomes = [g.strip() for g in trace_result.split(';') if g.strip()]
    print(f"  Traced genomes: {traced_genomes}")
    print(f"  Trace type: {'Node (multiple)' if len(traced_genomes) > 1 else 'Single genome'}")

    # Filter meta to only include genomes in tree
    tree_genomes = set(re.findall(r'GCF_\d+\.\d+', tree_nwk))
    filtered_meta = {k: v for k, v in meta.items() if k in tree_genomes}
    print(f"  Genomes in tree: {len(tree_genomes)}")
    print(f"  Matched meta: {len(filtered_meta)}")

    meta_json = json.dumps(filtered_meta, separators=(',', ':'))

    display_name = format_taxname(args.taxname)

    # ── Generate HTML ──
    print(f"[INFO] Generating HTML...")
    html = generate_html(tree_nwk, meta_json, trace_result, args.taxname,
                         typing_cols, vir_genes, res_genes, args.dist_unit)

    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)
    
    size_kb = os.path.getsize(args.output) / 1024
    print(f"[DONE] HTML output: {args.output} ({size_kb:.0f} KB)")

    # ── Generate JSON interface data ──
    json_output = os.path.splitext(args.output)[0] + '.data.json'

    interface_data = {
        "_comment": "mTracePatho Viewer interface data — for frontend web integration",
        "_version": "1.0",

        # Parameters
        "params": {
            "taxname":   args.taxname,
            "dist_unit": args.dist_unit,
        },

        # Processed data (ready for frontend JS consumption)
        "data": {
            "NEWICK_STR":  tree_nwk,
            "META_DATA":   filtered_meta,
            "TRACE_RAW":   trace_result,
            "TAX_NAME":    display_name,
            "TYPING_COLS": typing_cols,
            "VIR_GENES":   vir_genes,
            "RES_GENES":   res_genes,
            "DIST_UNIT":   args.dist_unit,
        },

        # Summary statistics
        "stats": {
            "tree_tips":        len(tree_genomes),
            "meta_matched":     len(filtered_meta),
            "traced_genomes":   traced_genomes,
            "trace_type":       "node" if len(traced_genomes) > 1 else "single",
            "typing_col_count": len(typing_cols),
            "vir_gene_count":   len(vir_genes),
            "res_gene_count":   len(res_genes),
        },
    }

    with open(json_output, 'w', encoding='utf-8') as f:
        json.dump(interface_data, f, ensure_ascii=False, indent=2)
    
    json_kb = os.path.getsize(json_output) / 1024
    print(f"[DONE] JSON output: {json_output} ({json_kb:.0f} KB)")
    print(f"  Open HTML in browser, or load JSON in your frontend app.")


if __name__ == '__main__':
    main()
