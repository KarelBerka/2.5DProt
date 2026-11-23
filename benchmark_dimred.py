#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
Comprehensive Benchmark of Dimensionality Reduction Algorithms in 2DProts.
Evaluates 12 layout algorithms across globular enzymes (Cytochromes P450, EC 1-7),
GPCRs, SLC & ABC transporters, Beta-barrels, 14-3-3, Hemoglobin, KcsA, Aquaporin,
and Extra-large complexes (F1-ATPase, GroEL/GroES) — both as single domains
and as whole multi-domain / multi-chain biological complexes.
"""

import os
import sys
import json

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
import time
import math
import gzip
import shutil
import logging
import warnings
from argparse import Namespace
from collections import OrderedDict

warnings.filterwarnings("ignore")

import numpy as np
import requests
import gemmi

import utils_dimred
import utils_sses_types
import utils_cluster
import utils_matrix2layout
import utils_angleslayout
import utils_generate_svg
import layout_loader
import superimpose2d

logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "benchmark_dimred_results")
CACHE_DIR = os.path.join(OUT_DIR, "cache")
SVG_DIR = os.path.join(OUT_DIR, "svgs")
JSON_DIR = os.path.join(OUT_DIR, "layouts")

os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(SVG_DIR, exist_ok=True)
os.makedirs(JSON_DIR, exist_ok=True)

METHODS = [
    ("evolution", "Evolution (Původní 2DProts)", "Stochastická mutace (35 000 kroků)"),
    ("pca", "PCA (Ortogonální 3D→2D)", "Lineární projekce 3D středů na 2 hlavní osy"),
    ("pcoa", "PCoA / Classical MDS", "Spektrální rozklad centrované Gramovy matice"),
    ("mds", "Metric MDS (SMACOF)", "Iterativní majorizace globálního metrického stressu"),
    ("nmmds", "Non-Metric MDS", "Monotónní regrese pořadí vzdáleností (Kruskal)"),
    ("sammon", "Sammon / Local wMDS", "Kvazi-Newtonova (L-BFGS-B) optimalizace lokální chyby <40 Å"),
    ("isomap", "Isomap (Geodetický)", "Nejkratší cesty po 3D kontaktním grafu + PCoA"),
    ("spectral", "Spectral (Laplacian)", "Vlastní vektory normalizovaného Laplaciánu grafu"),
    ("lle", "LLE (Locally Linear)", "Lokální lineární barycentrické rekonstrukční váhy"),
    ("kamada_kawai", "Kamada-Kawai (Pružinový)", "Fyzikální model soustavy pružin na grafu vzdáleností"),
    ("tsne", "t-SNE", "Student-t stochastické vnoření sousedů (perplexity)"),
    ("umap", "UMAP", "Uniform Manifold Approximation (fuzzy simpliciální komplex)"),
]

# Benchmark dataset specification:
# Each entry is either type="domain" (single CATH domain) or type="complex"
# (multi-domain / multi-chain biological assembly built from RCSB mmCIF + 2DProts domain templates).
BENCHMARK_TARGETS = [
    # 1. Cytochromy P450
    {
        "id": "1tqnA00",
        "name": "CYP3A4 (Cytochrom P450 3A4)",
        "category": "1. Globulární enzymy – Cytochromy P450",
        "kind": "domain",
        "pdb": "1tqn",
        "domains": [("A", "00", "1tqnA00")],
        "desc": "Lidský mikrosomální cytochrom P450 3A4 (EC 1.14.14.1) – 17 α-helixů (I-helix nad hemem) a 4 β-listy.",
    },
    {
        "id": "2hi4A00",
        "name": "CYP1A2 (Cytochrom P450 1A2)",
        "category": "1. Globulární enzymy – Cytochromy P450",
        "kind": "domain",
        "pdb": "2hi4",
        "domains": [("A", "00", "2hi4A00")],
        "desc": "Lidský cytochrom P450 1A2 (EC 1.14.14.1) – kompaktní P450 fold s úzkou aktivní dutinou.",
    },

    # 2. Globulární enzymy – EC 1 až EC 7
    {
        "id": "1adsA00",
        "name": "EC 1: Aldózareduktáza (TIM barel)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1ads",
        "domains": [("A", "00", "1adsA00")],
        "desc": "EC 1.1.1.21 (Oxidoreduktázy) – kanonický (β/α)8 TIM barel s 8 vnitřními β-vlákny a 8 vnějšími α-helixy.",
    },
    {
        "id": "1ldmA01",
        "name": "EC 1: Laktátdehydrogenáza (Rossmannova doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1ldm",
        "domains": [("A", "01", "1ldmA01")],
        "desc": "EC 1.1.1.27 (Oxidoreduktázy) – NAD-vazebná Rossmannova doména (3-vrstvý α/β/α sendvič).",
    },
    {
        "id": "1atpE01",
        "name": "EC 2: Proteinkináza A – N-terminální lalok (doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1atp",
        "domains": [("E", "01", "1atpE01")],
        "desc": "EC 2.7.11.11 (Transferázy) – malý N-terminální antiparalelní β-lalok cAMP-dependentní proteinkinázy.",
    },
    {
        "id": "1atpE02",
        "name": "EC 2: Proteinkináza A – C-terminální lalok (doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1atp",
        "domains": [("E", "02", "1atpE02")],
        "desc": "EC 2.7.11.11 (Transferázy) – velký C-terminální helikální katalytický lalok proteinkinázy A.",
    },
    {
        "id": "1atp_complex",
        "name": "EC 2: Proteinkináza A – Celý dvoudoménový komplex",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "complex",
        "pdb": "1atp",
        "use_au": True,
        "domains": [("E", "01", "1atpE01"), ("E", "02", "1atpE02")],
        "desc": "EC 2.7.11.11 – kompletní katalytická podjednotka cAMP-dependentní proteinkinázy (N-lalok + C-lalok).",
    },
    {
        "id": "1lz1A00",
        "name": "EC 3: Lidský lysozym (doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1lz1",
        "domains": [("A", "00", "1lz1A00")],
        "desc": "EC 3.2.1.17 (Hydrolázy) – klasický globulární α+β enzym štěpící peptidoglykan.",
    },
    {
        "id": "2ca2A00",
        "name": "EC 4: Karbonanhydráza II (doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "2ca2",
        "domains": [("A", "00", "2ca2A00")],
        "desc": "EC 4.2.1.1 (Lyázy) – lidská karbonanhydráza II s rozsáhlým centrálním 10-vláknovým zkrouceným β-listem.",
    },
    {
        "id": "1timA00",
        "name": "EC 5: Triózafosfátizomeráza – Monomer (doména)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1tim",
        "domains": [("A", "00", "1timA00")],
        "desc": "EC 5.3.1.1 (Izomerázy) – samostatná podjednotka triózafosfátizomerázy (TIM barel).",
    },
    {
        "id": "1tim_dimer",
        "name": "EC 5: Triózafosfátizomeráza – Homodimer (komplex)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "complex",
        "pdb": "1tim",
        "use_au": True,
        "domains": [("A", "00", "1timA00"), ("B", "00", "1timA00")],
        "desc": "EC 5.3.1.1 (Izomerázy) – kompletní biologický homodimer dvou interagujících TIM barelů (řetězce A + B).",
    },
    {
        "id": "1aszA01",
        "name": "EC 6: Aspartyl-tRNA syntetáza – Antikodonová doména",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1asz",
        "domains": [("A", "01", "1aszA01")],
        "desc": "EC 6.1.1.12 (Ligázy) – N-terminální OB-fold β-barel doména rozpoznávající antikodon tRNA.",
    },
    {
        "id": "1aszA02",
        "name": "EC 6: Aspartyl-tRNA syntetáza – Katalytická doména",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "domain",
        "pdb": "1asz",
        "domains": [("A", "02", "1aszA02")],
        "desc": "EC 6.1.1.12 (Ligázy) – C-terminální katalytická doména třídy II aminoacyl-tRNA syntetázy.",
    },
    {
        "id": "1asz_dimer",
        "name": "EC 6: Aspartyl-tRNA syntetáza – Homodimer (4-doménový komplex)",
        "category": "2. Globulární enzymy – Zástupci EC 1–7",
        "kind": "complex",
        "pdb": "1asz",
        "use_au": True,
        "domains": [
            ("A", "01", "1aszA01"), ("A", "02", "1aszA02"),
            ("B", "01", "1aszA01"), ("B", "02", "1aszA02"),
        ],
        "desc": "EC 6.1.1.12 (Ligázy) – kompletní homodimerní komplex (2× antikodonová + 2× katalytická doména).",
    },

    # 3. GPCR
    {
        "id": "1f88A00",
        "name": "GPCR: Bovinní rodopsin – 7TM monomer (doména)",
        "category": "3. GPCR (Receptory spřažené s G-proteiny)",
        "kind": "domain",
        "pdb": "1f88",
        "domains": [("A", "00", "1f88A00")],
        "desc": "Klasický GPCR třídy A – 7 transmembránových α-helixů (TM1–TM7) + cytoplazmatický helix H8.",
    },
    {
        "id": "1f88_dimer",
        "name": "GPCR: Bovinní rodopsin – Dimer (komplex A+B)",
        "category": "3. GPCR (Receptory spřažené s G-proteiny)",
        "kind": "complex",
        "pdb": "1f88",
        "use_au": True,
        "domains": [("A", "00", "1f88A00"), ("B", "00", "1f88A00")],
        "desc": "Dvojice 7TM podjednotek bovinního rodopsinu v asymetrické jednotce (řetězce A + B).",
    },
    {
        "id": "2rh1A02",
        "name": "GPCR: β2-adrenergní receptor – 7TM doména",
        "category": "3. GPCR (Receptory spřažené s G-proteiny)",
        "kind": "domain",
        "pdb": "2rh1",
        "domains": [("A", "02", "2rh1A02")],
        "desc": "Lidský β2-adrenergní receptor – transmembránové 7TM jádro receptoru.",
    },
    {
        "id": "2rh1_complex",
        "name": "GPCR: β2-adrenergní receptor + T4-lysozym (celý komplex)",
        "category": "3. GPCR (Receptory spřažené s G-proteiny)",
        "kind": "complex",
        "pdb": "2rh1",
        "use_au": True,
        "domains": [("A", "01", "2rh1A01"), ("A", "02", "2rh1A02")],
        "desc": "Dvoudoménový chimérický komplex lidského β2-adrenergního receptoru (7TM) a fúzního T4-lysozymu.",
    },

    # 4. Transportéry – SLC a ABC
    {
        "id": "1okcA00",
        "name": "SLC: Mitochondriální ADP/ATP přenašeč (SLC25A4)",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "domain",
        "pdb": "1okc",
        "domains": [("A", "00", "1okcA00")],
        "desc": "SLC25A4 (ANT1) – 6TM transmembránový kardinální transportér s trojčetnou vnitřní pseudo-symetrií.",
    },
    {
        "id": "1pv6A01",
        "name": "SLC / MFS: Laktózová permeáza LacY – N-doména (6TM)",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "domain",
        "pdb": "1pv6",
        "domains": [("A", "01", "1pv6A01")],
        "desc": "N-terminální 6-helixový svazek (TM1–TM6) MFS symportéru LacY.",
    },
    {
        "id": "1pv6A02",
        "name": "SLC / MFS: Laktózová permeáza LacY – C-doména (6TM)",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "domain",
        "pdb": "1pv6",
        "domains": [("A", "02", "1pv6A02")],
        "desc": "C-terminální 6-helixový svazek (TM7–TM12) MFS symportéru LacY.",
    },
    {
        "id": "1pv6_complex",
        "name": "SLC / MFS: Laktózová permeáza LacY – Celý 12TM komplex (N+C)",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "complex",
        "pdb": "1pv6",
        "use_au": True,
        "domains": [("A", "01", "1pv6A01"), ("A", "02", "1pv6A02")],
        "desc": "Kompletní 12TM MFS transportér LacY složený z N-terminální a C-terminální 6TM domény.",
    },
    {
        "id": "2hydA01",
        "name": "ABC (EC 7): Transportér Sav1866 – Transmembránová TMD doména",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "domain",
        "pdb": "2hyd",
        "domains": [("A", "01", "2hydA01")],
        "desc": "EC 7.6.2.- (ABC translokáza) – samostatná transmembránová doména (TMD) s dlouhými vnitrobuněčnými helixy.",
    },
    {
        "id": "2hydA02",
        "name": "ABC (EC 7): Transportér Sav1866 – Nukleotid-vázající NBD doména",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "domain",
        "pdb": "2hyd",
        "domains": [("A", "02", "2hydA02")],
        "desc": "EC 7.6.2.- (ABC translokáza) – ATPázová nukleotid-vázající doména (NBD).",
    },
    {
        "id": "2hyd_dimer",
        "name": "ABC (EC 7): Transportér Sav1866 – Celý homodimer (2×TMD + 2×NBD)",
        "category": "4. Membránové transportéry – SLC a ABC (EC 7)",
        "kind": "complex",
        "pdb": "2hyd",
        "use_au": True,
        "domains": [
            ("A", "01", "2hydA01"), ("A", "02", "2hydA02"),
            ("B", "01", "2hydA01"), ("B", "02", "2hydA02"),
        ],
        "desc": "Kompletní 4-doménový ABC exportér Sav1866 s propletenými (domain-swapped) transmembránovými helixy.",
    },

    # 5. Beta-barely
    {
        "id": "1bxwA00",
        "name": "Beta-barel (TM): OmpA (8-vláknový transmembránový barel)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "1bxw",
        "domains": [("A", "00", "1bxwA00")],
        "desc": "Vnější membránový protein A (E. coli) – kompaktní uzavřený 8-vláknový antiparalelní transmembránový β-barel.",
    },
    {
        "id": "1emaA00",
        "name": "Beta-barel (solubilní): GFP (11-vláknový solubilní barel)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "1ema",
        "domains": [("A", "00", "1emaA00")],
        "desc": "Zelený fluorescenční protein (GFP) – 11-vláknový solubilní β-barel obklopující centrální chromoforový helix.",
    },
    {
        "id": "1rbpA00",
        "name": "Beta-barel (solubilní lipokalin): Retinol-vázající protein RBP4 s retinolem (RTL)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "1rbp",
        "domains": [("A", "00", "1rbpA00")],
        "desc": "Solubilní krevní lipokalin (RBP4) – 8-vláknový netransmembránový antiparalelní β-barel uzavírající hydrofobní molekulu vitaminu A (retinol, RTL) uvnitř dutiny barelu.",
    },
    {
        "id": "1cbsA00",
        "name": "Beta-barel (solubilní iLBP): Protein vázající kyselinu retinovou CRABP-II (s REA)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "1cbs",
        "domains": [("A", "00", "1cbsA00")],
        "desc": "Vnitrobuněčný netransmembránový 10-vláknový β-barel (β-clam fold) s navázanou kyselinou all-trans-retinovou (REA) uvnitř dutiny.",
    },
    {
        "id": "1stpA00",
        "name": "Beta-barel (solubilní): Streptavidin – Monomer s biotinem (BTN)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "1stp",
        "domains": [("A", "00", "1stpA00")],
        "desc": "Solubilní 8-vláknový antiparalelní β-barel streptavidinu s vysokoafinitně vázaným biotinem (BTN, vitamin B7) v kapse barelu.",
    },
    {
        "id": "1stp_tetramer",
        "name": "Beta-barel (solubilní komplex): Streptavidin – Homotetramer (4×8 β-vláken + 4×BTN)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "complex",
        "pdb": "1stp",
        "use_au": False,
        "domains": [("A", "00", "1stpA00"), ("B", "00", "1stpA00"), ("C", "00", "1stpA00"), ("D", "00", "1stpA00")],
        "desc": "Kompletní solubilní homotetramer streptavidinu (D2 symetrie, 4× 8-vláknový netransmembránový β-barel) se 4 molekulami biotinu (BTN).",
    },
    {
        "id": "2omfA00",
        "name": "Beta-barel (TM): Porin OmpF – Monomer (16-vláknový barel)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "domain",
        "pdb": "2omf",
        "domains": [("A", "00", "2omfA00")],
        "desc": "Samostatná podjednotka bakteriálního porinu OmpF tvořící široký 16-vláknový transmembránový β-barel.",
    },
    {
        "id": "2omf_trimer",
        "name": "Beta-barel (TM komplex): Porin OmpF – Biologický trimer (3×16 = 48 β-vláken)",
        "category": "5. Beta-barely (Transmembránové i solubilní, monomery i trimer)",
        "kind": "complex",
        "pdb": "2omf",
        "use_au": False,
        "domains": [("A", "00", "2omfA00"), ("B", "00", "2omfA00"), ("C", "00", "2omfA00")],
        "desc": "Kompletní biologický trimer porinu OmpF se třemi 16-vláknovými β-barely uspořádanými kolem trojčetné osy.",
    },

    # 6. 14-3-3 protein
    {
        "id": "1qjbA00",
        "name": "14-3-3ζ protein – Monomer (doména)",
        "category": "6. Regulační a vazebné proteiny – 14-3-3",
        "kind": "domain",
        "pdb": "1qjb",
        "domains": [("A", "00", "1qjbA00")],
        "desc": "Jedna podjednotka lidského 14-3-3ζ proteinu (9 antiparalelních α-helixů tvořících amfipatický žlábek).",
    },
    {
        "id": "1qjb_dimer",
        "name": "14-3-3ζ protein – Homodimer (komplex A+B)",
        "category": "6. Regulační a vazebné proteiny – 14-3-3",
        "kind": "complex",
        "pdb": "1qjb",
        "use_au": True,
        "domains": [("A", "00", "1qjbA00"), ("B", "00", "1qjbA00")],
        "desc": "Biologický homodimer 14-3-3ζ ve tvaru W / podkovy se dvěma vazebnými žlábky pro fosfopeptidy.",
    },

    # 7. Hemoglobin
    {
        "id": "2hhbA00",
        "name": "Hemoglobin – α-podjednotka (doména)",
        "category": "7. Globiny – Hemoglobin (Domény a α2β2 tetramer)",
        "kind": "domain",
        "pdb": "2hhb",
        "domains": [("A", "00", "2hhbA00")],
        "desc": "Lidský hemoglobin, α-řetězec – klasický myoglobinový/globinový 3-vrstvý ortogonální helikální svazek.",
    },
    {
        "id": "2hhbB00",
        "name": "Hemoglobin – β-podjednotka (doména)",
        "category": "7. Globiny – Hemoglobin (Domény a α2β2 tetramer)",
        "kind": "domain",
        "pdb": "2hhb",
        "domains": [("B", "00", "2hhbB00")],
        "desc": "Lidský hemoglobin, β-řetězec – globinový fold s hemovou kapsou.",
    },
    {
        "id": "2hhb_tetramer",
        "name": "Hemoglobin – Celý α2β2 tetramer (komplex A+B+C+D)",
        "category": "7. Globiny – Hemoglobin (Domény a α2β2 tetramer)",
        "kind": "complex",
        "pdb": "2hhb",
        "use_au": True,
        "domains": [
            ("A", "00", "2hhbA00"), ("B", "00", "2hhbB00"),
            ("C", "00", "2hhbA00"), ("D", "00", "2hhbB00"),
        ],
        "desc": "Kompletní kvartérní tetramer lidského deoxyhemoglobinu (α1, β1, α2, β2).",
    },

    # 8. KcsA K+ kanál
    {
        "id": "1bl8A00",
        "name": "Draslíkový kanál KcsA – Jedna podjednotka (doména)",
        "category": "8. Iontové kanály – KcsA (Doména a tetramerní pór)",
        "kind": "domain",
        "pdb": "1bl8",
        "domains": [("A", "00", "1bl8A00")],
        "desc": "Samostatná podjednotka K+ kanálu KcsA (vnější TM1 helix, krátký pore-helix a vnitřní TM2 helix).",
    },
    {
        "id": "1bl8_tetramer",
        "name": "Draslíkový kanál KcsA – Celý biologický tetramer (komplex A+B+C+D)",
        "category": "8. Iontové kanály – KcsA (Doména a tetramerní pór)",
        "kind": "complex",
        "pdb": "1bl8",
        "use_au": True,
        "domains": [
            ("A", "00", "1bl8A00"), ("B", "00", "1bl8A00"),
            ("C", "00", "1bl8A00"), ("D", "00", "1bl8A00"),
        ],
        "desc": "Kompletní symetrický tetramer K+ kanálu KcsA obklopující centrální selektivitní filtr a vodní dutinu.",
    },

    # 9. Akvaporin AQP1
    {
        "id": "1j4nA00",
        "name": "Akvaporin AQP1 – Monomer (doména)",
        "category": "9. Vodní kanály – Akvaporin AQP1 (Doména a tetramer)",
        "kind": "domain",
        "pdb": "1j4n",
        "domains": [("A", "00", "1j4nA00")],
        "desc": "Monomerní vodní kanál AQP1 – 6 transmembránových helixů a 2 poloviční NPA helixy (přesýpací hodiny).",
    },
    {
        "id": "1j4n_tetramer",
        "name": "Akvaporin AQP1 – Celý biologický tetramer (4× kanál)",
        "category": "9. Vodní kanály – Akvaporin AQP1 (Doména a tetramer)",
        "kind": "complex",
        "pdb": "1j4n",
        "use_au": False,
        "domains": [
            ("A", "00", "1j4nA00"), ("B", "00", "1j4nA00"),
            ("C", "00", "1j4nA00"), ("D", "00", "1j4nA00"),
        ],
        "desc": "Biologický tetramer akvaporinu AQP1 v membráně (4 vodní kanály kolem čtyřčetné osy symetrie).",
    },

    # 10. Extra velké systémy
    {
        "id": "1bmfA02",
        "name": "Extra velký systém (EC 7): F1-ATPáza – Katalytická α/β doména",
        "category": "10. Extra velké multipodjednotkové komplexy (F1-ATPáza, GroEL)",
        "kind": "domain",
        "pdb": "1bmf",
        "domains": [("A", "02", "1bmfA02")],
        "desc": "EC 7.1.2.2 – nukleotid-vázající RecA-like α/β doména jedné podjednotky F1-ATPázy.",
    },
    {
        "id": "1bmf_complex",
        "name": "Extra velký systém (EC 7): F1-ATPáza – Celý motor α3β3γ (19 domén!)",
        "category": "10. Extra velké multipodjednotkové komplexy (F1-ATPáza, GroEL)",
        "kind": "complex",
        "pdb": "1bmf",
        "use_au": True,
        "domains": [
            ("A", "01", "1bmfA01"), ("A", "02", "1bmfA02"), ("A", "03", "1bmfA03"),
            ("B", "01", "1bmfA01"), ("B", "02", "1bmfA02"), ("B", "03", "1bmfA03"),
            ("C", "01", "1bmfA01"), ("C", "02", "1bmfA02"), ("C", "03", "1bmfA03"),
            ("D", "01", "1bmfD01"), ("D", "02", "1bmfD02"), ("D", "03", "1bmfD03"),
            ("E", "01", "1bmfD01"), ("E", "02", "1bmfD02"), ("E", "03", "1bmfD03"),
            ("F", "01", "1bmfD01"), ("F", "02", "1bmfD02"), ("F", "03", "1bmfD03"),
            ("G", "00", "1bmfG00"),
        ],
        "desc": "EC 7.1.2.2 – kompletní mitochondriální F1-ATPáza (hexamer α3β3 + centrální rotační hřídel γ, 19 domén, >150 SSE).",
    },
    {
        "id": "1aonA01",
        "name": "Extra velký systém: Chaperonin GroEL – Ekvatoriální ATPázová doména",
        "category": "10. Extra velké multipodjednotkové komplexy (F1-ATPáza, GroEL)",
        "kind": "domain",
        "pdb": "1aon",
        "domains": [("A", "01", "1aonA01")],
        "desc": "Ekvatoriální doména podjednotky chaperoninu GroEL zajišťující vazbu ATP a kontakt mezi kruhy.",
    },
    {
        "id": "1aon_ring",
        "name": "Extra velký systém: Chaperonin GroEL – Heptamerní kruh (7 podjednotek)",
        "category": "10. Extra velké multipodjednotkové komplexy (F1-ATPáza, GroEL)",
        "kind": "complex",
        "pdb": "1aon",
        "use_au": True,
        "domains": [
            ("A", "01", "1aonA01"), ("B", "01", "1aonA01"), ("C", "01", "1aonA01"),
            ("D", "01", "1aonA01"), ("E", "01", "1aonA01"), ("F", "01", "1aonA01"),
            ("G", "01", "1aonA01"),
        ],
        "desc": "Sedmičetný symetrický cis-kruh ekvatoriálních domén chaperoninu GroEL (řetězce A–G, ~100 SSE).",
    },

    # 11. Na+/K+-ATPáza, Kompletní protonová pumpa a Komplexy dýchacího řetězce (s ligandy)
    {
        "id": "2zxeA02",
        "name": "Na+/K+-ATPáza (EC 7): 10TM transmembránová doména α-podjednotky (s K+ a CLR)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "domain",
        "pdb": "2zxe",
        "domains": [("A", "02", "2zxeA02")],
        "desc": "EC 7.2.2.13 – 10TM transmembránové jádro α-podjednotky sodno-draselné pumpy se 2 vazebnými místy pro K+ v membráně, regulačním K+ a cholesterolem (CLR).",
    },
    {
        "id": "2zxe_complex",
        "name": "Na+/K+-ATPáza (EC 7): Kompletní αβ-heterokomplex (6 domén, K+, MgF4, Mg2+, CLR)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "complex",
        "pdb": "2zxe",
        "use_au": True,
        "domains": [
            ("A", "01", "2zxeA01"), ("A", "02", "2zxeA02"),
            ("A", "03", "2zxeA03"), ("A", "04", "2zxeA04"),
            ("B", "01", "2zxeB01"), ("B", "02", "2zxeB02"),
        ],
        "desc": "EC 7.2.2.13 – kompletní sodno-draselná pumpa (α-podjednotka: A, N, P a 10TM doména + β-podjednotka) s vázanými ionty K+, Mg2+, fosfátovým analogem MgF4, cholesterolem a glykany.",
    },
    {
        "id": "1c3wA00",
        "name": "Protonová pumpa: Bakteriorhodopsin – kompletní 7TM pumpa s retinalem (RET)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "domain",
        "pdb": "1c3w",
        "domains": [("A", "00", "1c3wA00")],
        "desc": "Světlem poháněná 7TM protonová pumpa (Halobacterium salinarum) s kovalentně vázaným retinalem (RET), skvalenem (SQU) a membránovými lipidy (LI1).",
    },
    {
        "id": "2xnd_complex",
        "name": "Kompletní protonová pumpa (EC 7): F1c8-ATP syntáza – α3β3γδε + TM rotor c8 (31 domén!)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "complex",
        "pdb": "2xnd",
        "use_au": True,
        "domains": [
            ("A", "01", "2xndA01"), ("A", "02", "2xndA02"), ("A", "03", "2xndA03"),
            ("B", "01", "2xndA01"), ("B", "02", "2xndA02"), ("B", "03", "2xndA03"),
            ("C", "01", "2xndA01"), ("C", "02", "2xndA02"), ("C", "03", "2xndA03"),
            ("D", "01", "2xndD01"), ("D", "02", "2xndD02"), ("D", "03", "2xndD03"),
            ("E", "01", "2xndD01"), ("E", "02", "2xndD02"), ("E", "03", "2xndD03"),
            ("F", "01", "2xndD01"), ("F", "02", "2xndD02"), ("F", "03", "2xndD03"),
            ("G", "01", "2xndG01"), ("G", "02", "2xndG02"),
            ("H", "01", "2xndH01"), ("H", "02", "2xndH02"),
            ("I", "00", "2xndI00"),
            ("J", "00", "2xndJ00"), ("K", "00", "2xndJ00"), ("L", "00", "2xndJ00"), ("M", "00", "2xndJ00"),
            ("N", "00", "2xndJ00"), ("O", "00", "2xndJ00"), ("P", "00", "2xndJ00"), ("Q", "00", "2xndJ00"),
        ],
        "desc": "EC 7.1.2.2 – kompletní mitochondriální F1c8 protonová pumpa / ATP syntáza: katalytický hexamer α3β3 (s 5× ANP a 5× Mg2+), centrální rotor γδε a transmembránový protonový turbínový kruh c8 (řetězce J–Q).",
    },
    {
        "id": "1zoy_complex",
        "name": "Dýchací řetězec – Komplex II: Sukcinátdehydrogenáza (komplex A+B+C+D, FAD, 3×FeS, UQ1, HEM)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "complex",
        "pdb": "1zoy",
        "use_au": True,
        "domains": [
            ("A", "01", "1zoyA01"), ("A", "02", "1zoyA02"),
            ("A", "03", "1zoyA03"), ("A", "04", "1zoyA04"),
            ("B", "01", "1zoyB01"), ("B", "02", "1zoyB02"),
            ("C", "01", "1zoyC01"), ("C", "02", "1zoyC02"),
            ("D", "00", "1zoyD00"),
        ],
        "desc": "EC 1.3.5.1 – kompletní mitochondriální Komplex II dýchacího řetězce (9 domén) s kompletní kaskádou přenosu elektronů: FAD -> [2Fe-2S] (FES) -> [4Fe-4S] (SF4) -> [3Fe-4S] (F3S) -> ubichinon (UQ1) a transmembránový hem b (HEM).",
    },
    {
        "id": "1occA00",
        "name": "Dýchací řetězec – Komplex IV: Cytochrom c oxidáza – podjednotka I (12TM jádro, 2×HEA, CuB, Mg)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "domain",
        "pdb": "1occ",
        "domains": [("A", "00", "1occA00")],
        "desc": "EC 7.1.1.9 – katalytická 12TM podjednotka I cytochrom c oxidázy s nízkospinovým hemem a (HEA), vysokospinovým hemem a3 (HEA), měďnatým centrem CuB (CU) a iontem Mg2+.",
    },
    {
        "id": "1occ_complex",
        "name": "Dýchací řetězec – Komplex IV: Cytochrom c oxidáza – kompletní 13-podjednotkový komplex (A–M, 15 domén)",
        "category": "11. Bioenergetické pumpy a komplexy dýchacího řetězce (Na+/K+-ATPáza, Protonové pumpy, Komplex II a IV)",
        "kind": "complex",
        "pdb": "1occ",
        "use_au": True,
        "domains": [
            ("A", "00", "1occA00"),
            ("B", "01", "1occB01"), ("B", "02", "1occB02"),
            ("C", "01", "1occC01"), ("C", "02", "1occC02"),
            ("D", "00", "1occD00"), ("E", "00", "1occE00"), ("F", "00", "1occF00"),
            ("G", "00", "1occG00"), ("H", "00", "1occH00"), ("I", "00", "1occI00"),
            ("J", "00", "1occJ00"), ("K", "00", "1occK00"), ("L", "00", "1occL00"), ("M", "00", "1occM00"),
        ],
        "desc": "EC 7.1.1.9 – kompletní 13-podjednotkový mitochondriální Komplex IV (řetězce A–M, 15 domén) s binukleárním centrem CuA (2×CU), hemem a, hemem a3 (2×HEA), CuB (CU), Mg2+ a Zn2+.",
    },
]


class MockSSESStats:
    """Provides exact prime/sec/third SSE classification from a 2DProts h_type & size dict."""

    def __init__(self, h_type, sizes):
        self.label_to_type = dict(h_type)
        self.sizes = dict(sizes)

    def most_probable(self, sses_set):
        best_sse = None
        best_sz = -1
        for sse in sses_set:
            sz = self.sizes.get(sse, 0)
            if sz > best_sz:
                best_sz = sz
                best_sse = sse
        return best_sse if best_sse is not None else next(iter(sses_set))

    def determine_types(self, label, size=None):
        if label in self.label_to_type:
            return self.label_to_type[label]
        sz_map = size if size is not None else self.sizes
        if label in sz_map and utils_sses_types.is_nontrivial_sse(label, sz_map[label]):
            return 2
        return 3

    def is_prime(self, label):
        return self.label_to_type.get(label, 3) == 1

    def is_sec(self, label, size=None):
        return self.determine_types(label, size) == 2

    def is_third(self, label):
        return self.label_to_type.get(label, 3) == 3

    def get_avg_size(self, label):
        return self.sizes.get(label, 5)


def fetch_domain_json(dom_id):
    """Fetch and cache domain layout.json from 2dprots.ncbr.muni.cz."""
    cache_path = os.path.join(CACHE_DIR, f"domain_{dom_id}.json")
    if os.path.exists(cache_path):
        with open(cache_path, "r", encoding="utf-8") as f:
            return json.load(f)
    url = f"https://2dprots.ncbr.muni.cz/files/domain/{dom_id}/layout.json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    return data


def fetch_pdb_structure(pdb_id, use_au=True):
    """Fetch and cache mmCIF from RCSB and return the gemmi Model (AU or Assembly 1)."""
    cache_path = os.path.join(CACHE_DIR, f"pdb_{pdb_id}.cif.gz")
    if not os.path.exists(cache_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
        r = requests.get(url, timeout=45)
        r.raise_for_status()
        with open(cache_path, "wb") as f:
            f.write(r.content)
    with gzip.open(cache_path, "rb") as f:
        cif_text = f.read().decode("utf-8", errors="replace")
    doc = gemmi.cif.read_string(cif_text)
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    if not use_au and len(st.assemblies) > 0:
        return gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.Short)
    return st[0]


def get_polymer_chains(model):
    """Return dict mapping chain_name -> gemmi.Chain (or ResidueSpan) with >15 polymer residues."""
    out = {}
    for ch in model:
        poly = ch.get_polymer()
        if len(poly) >= 15 and ch.name not in out:
            out[ch.name] = poly
        elif len(ch) >= 15 and ch.name not in out:
            # In generated assemblies, setup_entities might not mark polymer subchains
            ca_count = sum(1 for res in ch if "CA" in res)
            if ca_count >= 15:
                out[ch.name] = ch
    return out


def extract_ca_map(chain_obj):
    """Return dict mapping seqid.num -> np.array([x, y, z]) of CA atom."""
    ca_map = {}
    for res in chain_obj:
        if "CA" in res:
            pos = res["CA"][0].pos
            ca_map[res.seqid.num] = np.array([pos.x, pos.y, pos.z], dtype=float)
    return ca_map


def transform_domain_3d_to_chain(dom_data, ca_map):
    """
    Map the domain's 3dlayout vectors into the true 3D coordinate frame of `ca_map`
    using rigid Kabsch/superposition on SSE endpoints/midpoints + direct CA line segments.
    """
    hr = dom_data["helices_residues"]
    l3d = dom_data["3dlayout"]
    src_pts = []
    dst_pts = []
    direct_vectors = {}

    # Detect if domain `helices_residues` uses 1-based sequential polymer numbering
    # while the CIF `ca_map` uses offset author numbering (e.g. 1bl8 starts at res 23 while hr starts at res 4)
    res_offset = 0
    if ca_map and hr:
        min_hr = min(min(pair) for pair in hr.values())
        min_ca = min(ca_map.keys())
        if min_hr <= 10 and min_ca > 15 and (min_ca - 1 + max(max(pair) for pair in hr.values())) in ca_map:
            res_offset = min_ca - 1

    for sse, res_pair in hr.items():
        if sse not in l3d:
            continue
        r1, r2 = min(res_pair) + res_offset, max(res_pair) + res_offset
        coords = [ca_map[r] for r in range(r1, r2 + 1) if r in ca_map]
        if len(coords) >= 2:
            pts = np.array(coords)
            # Fit 3D line segment to CA atoms, oriented from r1 to r2
            mean = pts.mean(axis=0)
            centered = pts - mean
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            direction = vt[0]
            if np.dot(pts[-1] - pts[0], direction) < 0:
                direction = -direction
            proj = centered @ direction
            p_start = mean + proj[0] * direction
            p_end = mean + proj[-1] * direction
            direct_vectors[sse] = [p_start.tolist(), p_end.tolist()]
            src_pts.append(l3d[sse][0])
            src_pts.append(l3d[sse][1])
            dst_pts.append(p_start.tolist())
            dst_pts.append(p_end.tolist())

    # If any SSEs missed CA atoms in the crystal structure, transform them from consensus via Kabsch
    if len(src_pts) >= 4:
        sup = gemmi.superpose_positions(
            [gemmi.Position(*p) for p in dst_pts],
            [gemmi.Position(*p) for p in src_pts],
        )
        for sse, vec in l3d.items():
            if sse not in direct_vectors:
                p0 = sup.transform.apply(gemmi.Position(*vec[0]))
                p1 = sup.transform.apply(gemmi.Position(*vec[1]))
                direct_vectors[sse] = [[p0.x, p0.y, p0.z], [p1.x, p1.y, p1.z]]
    else:
        for sse, vec in l3d.items():
            if sse not in direct_vectors:
                direct_vectors[sse] = [list(vec[0]), list(vec[1])]

    return direct_vectors


def reconstruct_sheet_triples(parts, l3d, rename_fn=None):
    """
    Reconstruct the true 3D H-bonded neighbor pairs [s1, s2, orient] for a beta-sheet/barrel
    cluster `parts` (since 2dprots layout.json stores `beta_connectivity` strings in BFS order).
    Builds a minimum-weight degree<=2 spanning path (and closes the cycle if it's a barrel).
    """
    if rename_fn is None:
        rename_fn = lambda x: x
    k = len(parts)
    if k < 2:
        return []
    import utils
    sub_vecs = {s: l3d[s] for s in parts}
    dmat = utils.get_distance_matrix4n(sub_vecs)

    # All candidate edges sorted by 3D line-to-line distance
    edges = []
    for i in range(k):
        for j in range(i + 1, k):
            s1, s2 = parts[i], parts[j]
            edges.append((float(dmat[s1][s2]), i, j))
    edges.sort(key=lambda x: x[0])

    deg = [0] * k
    parent = list(range(k))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    chosen = []
    for dist_val, i, j in edges:
        if deg[i] < 2 and deg[j] < 2 and find(i) != find(j):
            chosen.append((i, j, dist_val))
            deg[i] += 1
            deg[j] += 1
            parent[find(i)] = find(j)
            if len(chosen) == k - 1:
                break

    # If k >= 5 (potential beta-barrel), check if the two degree-1 endpoints are close in 3D (< 6.5 Å)
    if k >= 5:
        endpoints = [idx for idx in range(k) if deg[idx] == 1]
        if len(endpoints) == 2:
            i, j = endpoints[0], endpoints[1]
            d_close = float(dmat[parts[i]][parts[j]])
            if d_close < 6.5:
                chosen.append((i, j, d_close))

    triples = []
    for i, j, _ in chosen:
        s1, s2 = parts[i], parts[j]
        v1 = np.array(l3d[s1][1]) - np.array(l3d[s1][0])
        v2 = np.array(l3d[s2][1]) - np.array(l3d[s2][0])
        orient = 1 if float(np.dot(v1, v2)) >= 0 else -1
        triples.append([rename_fn(s1), rename_fn(s2), orient])
    return triples


def build_target_annotation(target):
    """
    Build an Annotation JSON file, Annotation instance, MockSSESStats, and reference
    start_layout / start_angles for either a single domain or a multi-domain complex.
    """
    tid = target["id"]
    ann_path = os.path.join(CACHE_DIR, f"annotation_{tid}.json")

    if target["kind"] == "domain":
        chain_id, dom_id, dom_code = target["domains"][0]
        dom_data = fetch_domain_json(dom_code)
        l3d = dom_data["3dlayout"]
        hr = dom_data["helices_residues"]
        sizes = dom_data["size"]
        h_type = dom_data["h_type"]
        colors = dom_data.get("color", {})
        bconn_clusters = dom_data.get("beta_connectivity", [])
        ref_layout = {k: tuple(v) for k, v in dom_data.get("layout", {}).items()}
        ref_angles = dict(dom_data.get("angles", {}))

        sses_list = []
        for sse in l3d:
            if sse not in hr or sse not in sizes:
                continue
            r_start, r_end = min(hr[sse]), max(hr[sse])
            # Ensure r_end - r_start matches sizes[sse]
            r_end = r_start + int(sizes[sse])
            sses_list.append({
                "label": sse,
                "start": r_start,
                "end": r_end,
                "start_vector": l3d[sse][0],
                "end_vector": l3d[sse][1],
                "chain_id": chain_id,
                "rainbow": colors.get(sse, "grey"),
            })

        # Reconstruct pairwise beta_connectivity triples from cluster strings
        bconn_triples = []
        for cl in bconn_clusters:
            parts = [p for p in cl.split(":") if p in l3d]
            bconn_triples.extend(reconstruct_sheet_triples(parts, l3d))

        ann_json = {
            tid: {
                "secondary_structure_elements": sses_list,
                "beta_connectivity": bconn_triples,
            }
        }
        with open(ann_path, "w", encoding="utf-8") as f:
            json.dump(ann_json, f)

        ann = layout_loader.Annotation(ann_path, chain_id, dom_id)
        stats = MockSSESStats(h_type, sizes)
        return ann, stats, ref_layout, ref_angles, chain_id, dom_id

    else:
        # Whole complex: map each (chain, dom, template_dom_code) onto the 3D assembly
        pdb_id = target["pdb"]
        use_au = target.get("use_au", True)
        model = fetch_pdb_structure(pdb_id, use_au=use_au)
        poly_chains = get_polymer_chains(model)

        sses_list = []
        bconn_triples = []
        h_type_all = {}
        sizes_all = {}
        ndom = len(target["domains"])

        # Distinct subunit/domain spectrum colors ("s0" .. "s850")
        for d_idx, (ch_id, dom_id, dom_code) in enumerate(target["domains"]):
            dom_data = fetch_domain_json(dom_code)
            if ch_id in poly_chains:
                ca_map = extract_ca_map(poly_chains[ch_id])
                l3d = transform_domain_3d_to_chain(dom_data, ca_map)
            else:
                l3d = dom_data["3dlayout"]

            hr = dom_data["helices_residues"]
            sizes = dom_data["size"]
            h_type = dom_data["h_type"]
            bconn_clusters = dom_data.get("beta_connectivity", [])

            prefix = f"{d_idx + 1:02d}" if ndom >= 10 else f"{d_idx + 1}"
            dom_color = f"s{int(820.0 * d_idx / max(1, ndom))}"

            def rename_sse(s):
                if s.startswith("_"):
                    return f"{s[:2]}{prefix}_{s[2:]}"
                return f"{s[0]}{prefix}_{s[1:]}"

            for sse in l3d:
                if sse not in hr or sse not in sizes:
                    continue
                # For large complexes (>6 domains), keep only prime & secondary SSEs (skip tiny 3rd-level loops)
                ht_val = h_type.get(sse, 3)
                if ndom >= 4 and ht_val == 3:
                    continue
                new_label = rename_sse(sse)
                r_start, r_end = min(hr[sse]), max(hr[sse])
                r_end = r_start + int(sizes[sse])
                sses_list.append({
                    "label": new_label,
                    "start": r_start + d_idx * 2000,
                    "end": r_end + d_idx * 2000,
                    "start_vector": l3d[sse][0],
                    "end_vector": l3d[sse][1],
                    "chain_id": ch_id,
                    "rainbow": dom_color,
                })
                h_type_all[new_label] = 1 if ht_val == 1 else 2
                sizes_all[new_label] = int(sizes[sse])

            for cl in bconn_clusters:
                parts = [p for p in cl.split(":") if p in l3d and (ndom < 4 or h_type.get(p, 3) != 3)]
                bconn_triples.extend(reconstruct_sheet_triples(parts, l3d, rename_fn=rename_sse))

        ann_json = {
            tid: {
                "secondary_structure_elements": sses_list,
                "beta_connectivity": bconn_triples,
            }
        }
        with open(ann_path, "w", encoding="utf-8") as f:
            json.dump(ann_json, f)

        ann = layout_loader.Annotation(ann_path, "ALL", "ALL")
        stats = MockSSESStats(h_type_all, sizes_all)

        # Build a canonical PCoA reference frame for the complex so all 12 methods share one global orientation
        dmat_dict = ann.get_distances()
        sses_sorted = sorted(ann.sses)
        n_s = len(sses_sorted)
        D_s = np.zeros((n_s, n_s), dtype=float)
        for i in range(n_s):
            s_i = sses_sorted[i]
            for j in range(i + 1, n_s):
                s_j = sses_sorted[j]
                d_val = dmat_dict.get(s_i, {}).get(s_j, 20.0)
                D_s[i, j] = D_s[j, i] = d_val
        pcoa_coords = utils_dimred.embed_distance_matrix(D_s, method="pcoa")
        ref_layout = {
            sses_sorted[i]: (float(pcoa_coords[i, 0]), float(pcoa_coords[i, 1]))
            for i in range(n_s)
        }
        return ann, stats, ref_layout, {}, "ALL", "ALL"


def evaluate_layout_metrics(ann, l_obj):
    """
    Compute quantitative quality metrics comparing the 2D layout against 3D geometry:
    1. local_err: 2DProts size-weighted local distance error (< 40 Å)
    2. stress_1: Kruskal normalized global metric stress-1
    3. knn_pres: k-NN neighborhood preservation (k=4)
    4. conflicts: count of intersecting / too-close SSE pairs
    """
    shown_sses = [
        s for s in l_obj.layout
        if s in ann.res_sizes and l_obj.h_type.get(s, 3) in (1, 2)
    ]
    n = len(shown_sses)
    if n < 3:
        return {"local_err": 0.0, "stress_1": 0.0, "knn_pres": 1.0, "conflicts": 0, "n_sses": n}

    d3d_dict = ann.get_distances()
    D3 = np.zeros((n, n), dtype=float)
    D2 = np.zeros((n, n), dtype=float)
    weights = np.array([max(1.0, float(ann.res_sizes.get(s, 5))) for s in shown_sses])

    for i in range(n):
        si = shown_sses[i]
        pi = l_obj.layout[si]
        for j in range(i + 1, n):
            sj = shown_sses[j]
            pj = l_obj.layout[sj]
            d3 = float(d3d_dict[si][sj])
            d2 = math.hypot(pi[0] - pj[0], pi[1] - pj[1])
            D3[i, j] = D3[j, i] = d3
            D2[i, j] = D2[j, i] = d2

    # 1. 2DProts Local Weighted Error (< 40 Å) per SSE
    local_Cutoff = 40.0
    mask = (D3 < local_Cutoff) & (np.triu(np.ones((n, n), dtype=bool), k=1))
    if np.any(mask):
        w_ij = np.outer(weights, weights)[mask]
        local_err = float(np.sum(np.abs(D2[mask] - D3[mask]) * w_ij) / np.sum(w_ij))
    else:
        triu = np.triu_indices(n, k=1)
        local_err = float(np.mean(np.abs(D2[triu] - D3[triu])))

    # 2. Kruskal Stress-1
    triu = np.triu_indices(n, k=1)
    denom = float(np.sum(D3[triu] ** 2))
    stress_1 = float(np.sqrt(np.sum((D2[triu] - D3[triu]) ** 2) / max(denom, 1e-9)))

    # 3. k-NN preservation (k = min(4, n-1))
    k = min(4, n - 1)
    overlap_sum = 0.0
    for i in range(n):
        d3_row = D3[i].copy()
        d3_row[i] = np.inf
        d2_row = D2[i].copy()
        d2_row[i] = np.inf
        nn3 = set(np.argsort(d3_row)[:k])
        nn2 = set(np.argsort(d2_row)[:k])
        overlap_sum += len(nn3 & nn2) / float(k)
    knn_pres = float(overlap_sum / n)

    # 4. Approximate 2D segment conflict count (centers closer than 1.2 Å from different clusters)
    cl_map = {}
    for cl in ann.get_beta_connectivity_clusters():
        for s in cl.split(":"):
            cl_map[s] = cl
    conflicts = 0
    for i in range(n):
        for j in range(i + 1, n):
            if cl_map.get(shown_sses[i]) != cl_map.get(shown_sses[j]) and D2[i, j] < 1.2:
                conflicts += 1

    return {
        "local_err": round(local_err, 2),
        "stress_1": round(stress_1, 3),
        "knn_pres": round(knn_pres * 100.0, 1),
        "conflicts": conflicts,
        "n_sses": n,
    }


def render_svg_from_json(json_path, svg_path):
    """Render standard 2DProts SVG from a layout.json file using utils_generate_svg.compute."""
    args = Namespace(
        only="",
        OPACITY=1.0,
        transparency=False,
        multi=False,
        orotation=None,
        irotation=None,
        otemplate=None,
        itemplate=None,
        oligands=None,
        iligands=None,
        groups=False,
        add_descriptions="0",
        ligands=True,
        dominik=False,
        fill=True,
        data=None,
        debug=False,
        layouts=[json_path],
        output=svg_path,
    )
    utils_generate_svg.compute(args)


_TARGET_LIGANDS_3D_CACHE = {}


def get_cached_target_ligands_3d(target, ann):
    tid = target["id"]
    if tid not in _TARGET_LIGANDS_3D_CACHE:
        import utils_benchmark_ligands
        _TARGET_LIGANDS_3D_CACHE[tid] = utils_benchmark_ligands.extract_target_ligands_3d(target, ann, CACHE_DIR)
    return _TARGET_LIGANDS_3D_CACHE[tid]


def run_single_target_method(target, method_id, ann, stats, ref_layout, ref_angles, chain_id, dom_id):
    """Run one (target, method_id) pair, save layout.json & SVG, and return metrics."""
    import utils_benchmark_ligands
    tid = target["id"]
    json_path = os.path.join(JSON_DIR, f"{tid}_{method_id}.json")
    svg_path = os.path.join(SVG_DIR, f"{tid}_{method_id}.svg")

    l_obj = layout_loader.Layout(sses_stats=stats, annotation=ann, method=method_id)
    l_obj.apply_membrane = False
    t0 = time.perf_counter()
    # Run 100% blind from scratch (no start_layout) for every algorithm
    l_obj.compute(chain_id, dom_id, start_layout={}, start_angles={})
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    # Rigidly align finished 2D layout via Procrustes to ref_layout so all 12 methods share orientation
    common = [k for k in l_obj.layout if k in ref_layout]
    if len(common) >= 2:
        A = np.array([[l_obj.layout[k][0], l_obj.layout[k][1]] for k in common], dtype=float).T
        B = np.array([[ref_layout[k][0], ref_layout[k][1]] for k in common], dtype=float).T
        R, t = superimpose2d.optimal_rotation_translation(A, B, allow_mirror=True)
        det_R = float(np.linalg.det(R))
        for k in list(l_obj.layout.keys()):
            pt = np.array([[l_obj.layout[k][0]], [l_obj.layout[k][1]]], dtype=float)
            pt_new = R @ pt + t
            l_obj.layout[k] = (float(pt_new[0, 0]), float(pt_new[1, 0]))
        for k in list(l_obj.a_layout.keys()):
            ang_rad = math.radians(float(l_obj.a_layout[k]))
            v = np.array([math.sin(ang_rad), math.cos(ang_rad)], dtype=float)
            v_rot = R @ v
            l_obj.a_layout[k] = math.degrees(math.atan2(float(v_rot[0]), float(v_rot[1]))) % 360.0

    # Prevent utils_generate_svg from overwriting the DR layout with sequence-unrolled membrane columns
    l_obj.membrane = {"is_transmembrane": False, "angle_step": 0}
    ligs_3d = get_cached_target_ligands_3d(target, ann)
    ligs_2d = utils_benchmark_ligands.project_ligands_to_2d(ligs_3d, ann, l_obj.layout, l_obj.a_layout, mem_info=None)
    l_obj.ligands_data = ligs_2d
    l_obj.save(chain_id, dom_id, json_path)

    render_svg_from_json(json_path, svg_path)
    metrics = evaluate_layout_metrics(ann, l_obj)
    metrics["time_ms"] = round(elapsed_ms, 1)
    metrics["n_ligands"] = len(ligs_2d)
    return metrics


def main():
    only_targets = sys.argv[1:] if len(sys.argv) > 1 else None
    summary_path = os.path.join(OUT_DIR, "results_summary.json")
    prev_results = {}
    if os.path.exists(summary_path):
        try:
            with open(summary_path, "r", encoding="utf-8") as f:
                prev_results = json.load(f)
        except Exception:
            prev_results = {}

    # Warm up UMAP / Numba JIT once before timing so the first protein isn't penalized by JIT compilation
    print("Warming up Numba/UMAP JIT...", flush=True)
    _d_warm = np.array([
        [0.0, 4.0, 7.0, 9.0, 11.0, 13.0],
        [4.0, 0.0, 4.5, 8.0, 10.0, 12.0],
        [7.0, 4.5, 0.0, 4.2, 7.5, 9.5],
        [9.0, 8.0, 4.2, 0.0, 4.0, 7.0],
        [11.0, 10.0, 7.5, 4.0, 0.0, 4.5],
        [13.0, 12.0, 9.5, 7.0, 4.5, 0.0],
    ])
    _ = utils_dimred.embed_distance_matrix(_d_warm, method="umap")

    all_results = dict(prev_results)

    for target in BENCHMARK_TARGETS:
        tid = target["id"]
        if only_targets and tid not in only_targets:
            continue
        if not only_targets and tid in all_results and all_results[tid]["methods"].get("umap", {}).get("time_ms", 99999) < 5000:
            print(f"Skipping already-computed [{tid}] {target['name']}", flush=True)
            continue
        print(f"\n=== [{tid}] {target['name']} ({target['kind']}) ===", flush=True)
        ann, stats, ref_layout, ref_angles, chain_id, dom_id = build_target_annotation(target)
        n_clusters = len(ann.get_beta_connectivity_clusters())
        print(f"    SSEs: {len(ann.sses)}, Clusters: {n_clusters}", flush=True)

        all_results[tid] = {
            "target": target,
            "n_sses": len(ann.sses),
            "n_clusters": n_clusters,
            "methods": {},
        }
        for method_id, method_title, _ in METHODS:
            m = run_single_target_method(
                target, method_id, ann, stats, ref_layout, ref_angles, chain_id, dom_id
            )
            all_results[tid]["methods"][method_id] = m
            print(
                f"    - {method_id:13s}: time={m['time_ms']:7.1f} ms | "
                f"local_err={m['local_err']:5.2f} Å | stress={m['stress_1']:.3f} | "
                f"kNN={m['knn_pres']:4.1f}% | conflicts={m['conflicts']}",
                flush=True,
            )

        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False)

    print(f"\nSaved summary to {summary_path} ({len(all_results)} targets)", flush=True)


if __name__ == "__main__":
    main()
