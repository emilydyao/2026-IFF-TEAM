# -*- coding: utf-8 -*-
"""
Created on 2026-04-08 16:47:40
@author: eyao

gen_sys.py  —  Build or extend MARTINI hair surface simulation inputs.

Two operating modes
---
1. Build from scratch (default):
     Constructs graphene slab + top/bottom MEA monolayers, then optionally
     adds odorant and water vapor.

2. Extend an existing minimized hair structure (--input_hair):
     Loads a pre-minimized GRA+TMEA+BMEA .gro file and appends odorant
     and/or water beads on top.  Use this for the 8-system campaign:
       4 odorants (DMS / LIN / SYR / PHE) × 2 modes (vapor / droplet)

Output folder is named automatically:
  --odorant DMS  --droplet False  →  DMS_single/
  --odorant DMS  --droplet True   →  DMS_droplet/
  --odorant None                  →  hair_clean/

Each folder contains:
  system.gro          — all-atom coordinate file
  system.top          — master GROMACS topology
  forcefield.itp      — MARTINI v2.2 nonbonded parameters
  graphene.itp        — graphene molecule topology
  lipidchain_top.itp  — top MEA lipid chain topology
  lipidchain_bot.itp  — bottom MEA lipid chain topology
  odorant.itp       — odorant topology (if applicable)
  sysinfo.json        — machine-readable run metadata

Usage examples
  # Build from scratch, DMS vapor:
  python gen_sys.py --odorant DMS

  # Extend minimized hair, all 8 systems:
  for ODO in DMS LIN SYR PHE; do
    python gen_sys.py --input_hair ./hair_clean_v2/em3.gro --odorant $ODO --droplet False
    python gen_sys.py --input_hair ./hair_clean/nvt1.gro --odorant $ODO --droplet True
  done

  # No odorant, build from scratch:
  python gen_sys.py --odorant None --o hair_clean
"""

import argparse
import json
import math
import os
import shutil
import numpy as np
from collections import defaultdict
from gromacs_helpers import *
from ff_data import *

# Graphene Sheet Topography

def define_grid_dimensions(target_lx, target_ly, b):
    a1x    = b * math.sqrt(3)
    a2y    = 1.5 * b
    nx     = max(1, round(target_lx / a1x))
    ny_raw = round(target_ly / a2y)
    if ny_raw % 2 != 0:
        ny_lo = max(2, ny_raw - 1)
        ny_hi = ny_raw + 1
        ny = ny_lo if abs(ny_lo * a2y - target_ly) <= abs(ny_hi * a2y - target_ly) else ny_hi
    else:
        ny = max(2, ny_raw)
    return nx, ny, nx * a1x, ny * a2y


def build_sheet(b, nx, ny, z=0.0, resname="GRA"):
    if ny % 2 != 0:
        raise ValueError(f"ny must be even for PBC compatibility (got {ny}).")
    a1x = b * math.sqrt(3)
    a2x = b * math.sqrt(3) / 2
    a2y = 3 * b / 2
    Lx  = nx * a1x
    Ly  = ny * a2y
    beads    = []
    atom_idx = 1
    for m in range(nx):
        for n in range(ny):
            ox = (m * a1x + n * a2x) % Lx
            oy =  n * a2y
            beads.append(dict(resnum=1, resname=resname, atomname="C1",
                               atomnum=atom_idx,
                               x=round(ox, 6), y=round(oy, 6), z=round(z, 6)))
            atom_idx += 1
            beads.append(dict(resnum=1, resname=resname, atomname="C1",
                               atomnum=atom_idx,
                               x=round(ox, 6), y=round(oy + b, 6), z=round(z, 6)))
            atom_idx += 1
    return beads, Lx, Ly


def place_graphene(sheet_beads, n_layers, z_start, interlayer_dz=INTERLAYER_DZ):
    b    = 0.1418
    Lx   = max(a['x'] for a in sheet_beads) + b * math.sqrt(3) / 2
    Ly   = max(a['y'] for a in sheet_beads) + b
    dx_B = b * math.sqrt(3) / 2
    dy_B = b / 2
    slab = []
    for k in range(n_layers):
        z_layer = round(z_start + k * interlayer_dz, 6)
        is_B    = (k % 2 == 1)
        for a in sheet_beads:
            new           = dict(a)
            new['z']      = z_layer
            new['resnum'] = k + 1
            if is_B:
                new['x'] = round((a['x'] + dx_B) % Lx, 6)
                new['y'] = round((a['y'] + dy_B) % Ly, 6)
            slab.append(new)
    return slab

# Graphene force field topology helper

def write_gra_bonds(beads, b, lx, ly, tol=0.05):
    bonds = []
    n = len(beads)
    for i in range(n):
        for j in range(i + 1, n):
            dx = beads[i]['x'] - beads[j]['x']
            dy = beads[i]['y'] - beads[j]['y']
            dx -= round(dx / lx) * lx
            dy -= round(dy / ly) * ly
            d   = math.sqrt(dx * dx + dy * dy)
            if d < b * (1.0 + tol):
                bonds.append((beads[i]['atomnum'], beads[j]['atomnum']))
    return bonds


def write_gra_angles(bonds):
    adj  = defaultdict(list)
    for (a, b_) in bonds:
        adj[a].append(b_)
        adj[b_].append(a)
    angles = []
    seen   = set()
    for center in adj:
        nbrs = sorted(adj[center])
        for p in range(len(nbrs)):
            for q in range(p + 1, len(nbrs)):
                j, k = nbrs[p], nbrs[q]
                key  = (j, center, k)
                if key not in seen:
                    seen.add(key)
                    angles.append(key)
    return angles

# Triangular grafting lattice

def rotate_xy(positions, phi, cx, cy, Lx, Ly):
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    rotated = []
    for (x, y) in positions:
        dx, dy = x - cx, y - cy
        xr = cx + dx * cos_phi - dy * sin_phi
        yr = cy + dx * sin_phi + dy * cos_phi
        rotated.append((round(xr % Lx, 6), round(yr % Ly, 6)))
    return rotated


def make_triangular_grid(Lx, Ly, grafting_distance, x_offset=0.0, y_offset=0.0):
    ax     = grafting_distance
    ay     = grafting_distance * math.sqrt(3) / 2
    nx     = int(Lx // ax)
    ny     = int(Ly // ay)
    ax_eff = Lx / nx
    ay_eff = Ly / ny
    positions = []
    for iy in range(ny):
        y_pos   = (y_offset + (iy + 0.5) * ay_eff) % Ly
        x_shift = (ax_eff / 2.0) if (iy % 2 == 1) else 0.0
        for ix in range(nx):
            x_pos = (x_offset + x_shift + (ix + 0.5) * ax_eff) % Lx
            positions.append((round(x_pos, 6), round(y_pos, 6)))
    density = len(positions) / (Lx * Ly)
    print(f"  Grid: nx={nx} ny={ny} | ax_eff={ax_eff:.4f} nm | "
          f"N={len(positions)} chains | density={density:.3f} chains/nm²")
    return positions


def place_monolayer(Lx, Ly, z_surface, grafting_distance,
                    l_bond=BOND_LJ, flip_z=False,
                    x_offset=0.0, y_offset=0.0):
    sign        = -1.0 if flip_z else +1.0
    mea_resname = ("B" if flip_z else "T") + MEA_RESNAME

    positions = make_triangular_grid(Lx, Ly, grafting_distance, x_offset, y_offset)
    phi       = math.atan(math.sqrt(3) / 2 / 4.5)
    positions = rotate_xy(positions, phi, Lx / 2, Ly / 2, Lx, Ly)

    chain_beads = []
    for mol_idx, (cx, cy) in enumerate(positions, 1):
        for local_idx, (name, mtype, chg, mass) in enumerate(LIPID_BEADS, 1):
            dz = sign * (local_idx - 1) * l_bond
            chain_beads.append(dict(
                resnum   = mol_idx,
                resname  = mea_resname,
                atomname = name,
                atomnum  = 0,
                x=cx, y=cy,
                z=round(z_surface + dz, 6),
                _mtype=mtype, _charge=chg, _mass=mass,
            ))
    return chain_beads, len(positions), mea_resname

# Solvent placement
def place_plane(n, z_center, z_half_width, rng, lx, ly,
                min_sep=BOND_LJ, existing=None, max_tries=500, label=None):
    if existing is None:
        existing = []
    placed = []

    for mol_num in range(1, n + 1):
        for attempt in range(1, max_tries + 1):
            x = rng.uniform(0.0, lx)
            y = rng.uniform(0.0, ly)
            z = rng.uniform(z_center - z_half_width, z_center + z_half_width)

            too_close = False
            for (px, py, pz) in existing:
                dx = abs(x - px); dx = min(dx, lx - dx)
                dy = abs(y - py); dy = min(dy, ly - dy)
                dz = abs(z - pz)
                if math.sqrt(dx*dx + dy*dy + dz*dz) < min_sep:
                    too_close = True
                    break

            if not too_close:
                placed.append((x, y, z))
                existing.append((x, y, z))
                print(f"  [{label}] placed {mol_num}/{n}")
                break
        else:
            print(f"  WARNING [{label}]: could not place molecule {mol_num}/{n} "
                  f"after {max_tries} attempts (box crowded for min_sep={min_sep} nm)")

    return placed


def compute_mols(lx, ly, lz, T=298.0):
    MOLECULES_PER_BEAD  = 4
    WATER_MOLE_FRACTION = 0.04
    k_B     = 1.380649e-23
    P_atm   = 101325.0
    P_water = WATER_MOLE_FRACTION * P_atm
    n_mol_per_m3 = P_water / (k_B * T)
    V_box_m3     = lx * ly * lz * 1e-27
    n_total  = max(2, round(n_mol_per_m3 * V_box_m3 / MOLECULES_PER_BEAD))
    n_top    = (n_total + 1) // 2
    n_bot    = n_total // 2
    return n_top, n_bot


def add_odorant(resname, lx, ly, lz, n_top, n_bot, z_top, z_bot,
                mol_offset=0, droplet=False, max_tries=100):
    if droplet:
        beads_list, r = load_droplet(resname, lx, ly, z_top)
        resnum   = mol_offset + 1
        prev_res = beads_list[0]['resnum']
        for a in beads_list:
            if a['resnum'] != prev_res:
                resnum  += 1
                prev_res = a['resnum']
            a['resnum'] = resnum
        n_mols = resnum - mol_offset
        print(f"  Droplet placed: {len(beads_list)} beads, "
              f"~{n_mols} molecules, radius ~{r:.3f} nm")
        return beads_list, n_mols, []   # no xyz list for droplet mode

    # Vapor mode
    MIN_ODORANT_SEP = 1.3   
    model_beads, r_eff = load_single_odorant(resname)
    min_sep = max(round(2.0 * r_eff + 0.1, 2), MIN_ODORANT_SEP)
    print(f"  Odorant min_sep: {min_sep:.3f} nm  (r_eff={r_eff:.3f} nm)")

    rng      = np.random.default_rng()
    existing = []
    top_xyz  = place_plane(n_top, z_top, 0.5, rng, lx, ly,
                           min_sep=min_sep, existing=existing,
                           max_tries=max_tries, label=f"{resname}-top")
    bot_xyz  = place_plane(n_bot, z_bot, 0.5, rng, lx, ly,
                           min_sep=min_sep, existing=existing,
                           max_tries=max_tries, label=f"{resname}-bot")

    odor_beads = []
    for mol_idx, (cx, cy, cz) in enumerate(top_xyz + bot_xyz, 1):
        resnum = mol_offset + mol_idx
        for b in model_beads:
            odor_beads.append(dict(
                resnum   = resnum,
                resname  = resname,
                atomname = b['atomname'],
                atomnum  = 0,
                x = round((cx + b['x']) % lx, 6),
                y = round((cy + b['y']) % ly, 6),
                z = round( cz + b['z'],        6),
                _mtype  = b.get('_mtype',  'C1'),
                _charge = b.get('_charge', 0.0),
                _mass   = b.get('_mass',   72.0),
            ))

    placed_xyz = top_xyz + bot_xyz
    print(f"  Odorant beads placed: {len(odor_beads)}  "
          f"({len(top_xyz)} top + {len(bot_xyz)} bot molecules)")
    return odor_beads, len(placed_xyz), placed_xyz


def add_water(lx, ly, lz, n_top, n_bot, z_top, z_bot, T=298.0,
              existing=None, max_tries=500):
    print(f"\nAdding water vapor (100% RH, T={T:.1f} K):")
    print(f"  Total MARTINI W beads: {n_top + n_bot}  ({n_top} top + {n_bot} bot)")

    if existing is None:
        existing = []

    rng     = np.random.default_rng()
    top_xyz = place_plane(n_top, z_top, 0.5, rng, lx, ly, min_sep=0.47,
                          existing=existing, max_tries=max_tries, label="W-top")
    bot_xyz = place_plane(n_bot, z_bot, 0.5, rng, lx, ly, min_sep=0.47,
                          existing=existing, max_tries=max_tries, label="W-bot")

    water_beads = []
    for mol_idx, (x, y, z) in enumerate(top_xyz + bot_xyz, 1):
        water_beads.append(dict(resnum=mol_idx, resname="W", atomname="W",
                                atomnum=0,
                                x=round(x, 6), y=round(y, 6), z=round(z, 6)))

    n_placed = len(water_beads)
    print(f"  Beads placed: {n_placed}  ({len(top_xyz)} top + {len(bot_xyz)} bot)")
    return water_beads, n_placed

# Force-field helper

def collect_bead_types(all_beads):
    _ATOMNAME_TO_TYPE = {
        "W":   "P4",
        "WF":  "BP4",
        "CGH": "CGH",   # ghost anchor — must stay distinct from C1
    }
    seen  = {}
    order = []
    for a in all_beads:
        if a['atomname'] == "CGH":
            mtype = "CGH"
        else:
            mtype = a.get('_mtype') or _ATOMNAME_TO_TYPE.get(a['atomname'], 'C1')
        if mtype not in seen:
            seen[mtype] = True
            order.append(mtype)
    return order

# Directory / run-name helper

def make_run_dir(odorant, is_droplet, override=None):
    if override and override not in ("system", ""):
        return override
    if odorant is None or odorant.lower() == "none":
        return "hair_clean"
    suffix = "droplet" if is_droplet else "single"
    return f"{odorant}_{suffix}"


def _z_extents_from_gro(beads):
    extents = {}
    for a in beads:
        rn = a['resname']
        z  = a['z']
        if rn not in extents:
            extents[rn] = [z, z]
        else:
            if z < extents[rn][0]: extents[rn][0] = z
            if z > extents[rn][1]: extents[rn][1] = z
    return {k: tuple(v) for k, v in extents.items()}


def _mol_counts_from_gro(beads):
    counts  = {}
    prev_key = None
    for a in beads:
        key = (a['resnum'], a['resname'])
        if key != prev_key:
            counts[a['resname']] = counts.get(a['resname'], 0) + 1
            prev_key = key
    return counts

ODORANT_RESNAMES = set(ODORANT_TEMPLATES.keys())   

def write_index_ndx(filename, all_beads, odorant=None):
    if odorant is not None:
        odo_resnames = {odorant}
    else:
        odo_resnames = ODORANT_RESNAMES

    groups = {"MEA": [], "GRA": [], "WATER": [], "ODO": [], "System": []}

    for a in all_beads:
        idx = a['atomnum']
        rn  = a['resname']
        groups["System"].append(idx)
        if rn in ("TMEA", "BMEA"):
            groups["MEA"].append(idx)
        elif rn == "GRA":
            groups["GRA"].append(idx)
        elif rn in ("W", "WF"):
            groups["WATER"].append(idx)
        elif rn in odo_resnames:
            groups["ODO"].append(idx)
        elif rn in ("TMEA"):
            groups["TMEA"].append(idx)
        elif rn in ("BMEA"):
            groups["BMEA"].append(idx)

    def _write_group(f, name, indices):
        f.write(f"[ {name} ]\n")
        for i, idx in enumerate(indices):
            f.write(f"{idx:6d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        if indices and len(indices) % 15 != 0:
            f.write("\n")
        f.write("\n")

    with open(filename, 'w') as f:
        f.write("; GROMACS index file\n")
        f.write("; Generated by gen_sys.py — do not edit manually\n\n")
        _write_group(f, "MEA",    groups["MEA"])
        _write_group(f, "GRA",    groups["GRA"])
        _write_group(f, "WATER",  groups["WATER"])
        _write_group(f, "TMEA",    groups["TMEA"])
        _write_group(f, "BMEA",    groups["BMEA"])
        if groups["ODO"]:
            _write_group(f, "ODO", groups["ODO"])
        _write_group(f, "System", groups["System"])

    # Report
    print(f"Wrote: {filename}")
    print(f"  MEA   : {len(groups['MEA'])} beads")
    print(f"  GRA   : {len(groups['GRA'])} beads")
    print(f"  WATER : {len(groups['WATER'])} beads")
    if groups["ODO"]:
        print(f"  ODO   : {len(groups['ODO'])} beads  ({odorant or '/'.join(sorted(odo_resnames))})")
    else:
        print(f"  ODO   : (no odorant beads — group omitted)")
    print(f"  System: {len(groups['System'])} beads")


def main():
    parser = argparse.ArgumentParser(
        description="Build or extend MARTINI hair surface GROMACS inputs."
    )
    parser.add_argument("--input_hair", default=None,
                        help="Path to a pre-minimized hair .gro file "
                             "(e.g. ./hair_clean/nvt1.gro). "
                             "If given, graphene+lipid beads are loaded from "
                             "this file and only odorant/water are appended. "
                             "Box dimensions, n_layers, n_top, n_bot are all "
                             "inferred from the file; --lx/--ly/--lz/--n are ignored.")
    parser.add_argument("-b", type=float, default=BOND_LJ,
                        help=f"C-C bond length for graphene ITP (nm, default {BOND_LJ})")
    parser.add_argument("--lx", type=float, default=35,
                        help="Target box x (nm) — ignored with --input_hair")
    parser.add_argument("--ly", type=float, default=35,
                        help="Target box y (nm) — ignored with --input_hair")
    parser.add_argument("--lz", type=float, default=80,
                        help="Box z (nm) — ignored with --input_hair")
    parser.add_argument("--n", type=int, default=4,
                        help="Number of graphene layers — ignored with --input_hair")
    parser.add_argument("--g_gap", type=float, default=0.7,
                        help="Gap between graphene surface and ghost beads (nm)")
    parser.add_argument("--s_gap", type=float, default=2.0,
                        help="Gap between lipid tips and solvent placement zone (nm)")
    parser.add_argument("--d_graft", type=float, default=0.65,
                        help="Grafting distance between lipid chains (nm) — "
                             "ignored with --input_hair")
    parser.add_argument("--add_water", type=lambda x: x.lower() != 'false', default=True,
                        help="Add water vapor beads (default True)")
    parser.add_argument("--antifreeze", type=lambda x: x.lower() != 'false', default=True,
                        help="Replace 10% of water beads with antifreeze beads (default True)")                    
    parser.add_argument("--odorant", default="DMS",
                        help="Odorant resname: DMS / LIN / SYR / PHE / None")
    parser.add_argument("--droplet", type=lambda x: x.lower() == 'true', default=False,
                        help="Place odorant as a pre-formed droplet (default False)")
    parser.add_argument("--tries", type=int, default=500,
                        help="Max placement attempts per solvent molecule")
    parser.add_argument("--o", default="system",
                        help="Override output folder name")
    args = parser.parse_args()

    b          = args.b
    s_gap      = args.s_gap
    g_gap      = args.g_gap
    d_graft    = args.d_graft
    odorant    = None if args.odorant.lower() == "none" else args.odorant
    is_droplet = args.droplet
    tries      = args.tries

    # Create output folder
    run_dir = make_run_dir(odorant, is_droplet, args.o)
    os.makedirs(run_dir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"  Run folder : {run_dir}/")
    print(f"{'='*60}\n")

    # PATH A — load pre-minimized hair structure
    if args.input_hair:
        if any([args.lx != 50, args.ly != 50, args.lz != 80, args.n != 4]):
            print("  NOTE: --lx / --ly / --lz / --n are ignored "
                  "when --input_hair is set; box and layer count "
                  "are read from the .gro file.\n")

        print(f"Loading existing hair structure: {args.input_hair}")
        title_in, hair_beads, box = parse_gro(args.input_hair)
        lx, ly, lz = box
        print(f"  Box: {lx:.4f} x {ly:.4f} x {lz:.4f} nm  "
              f"({len(hair_beads)} beads)")

        extents    = _z_extents_from_gro(hair_beads)
        mol_counts_in = _mol_counts_from_gro(hair_beads)

        # Validate expected residues
        for rn in ("GRA", "TMEA", "BMEA"):
            if rn not in extents:
                raise ValueError(
                    f"Expected residue '{rn}' not found in {args.input_hair}. "
                    f"Found: {list(extents.keys())}"
                )

        z_gra_bottom = extents["GRA"][0]
        z_gra_top    = extents["GRA"][1]
        z_top_term   = extents["TMEA"][1]   # highest bead of top monolayer
        z_bot_tip    = extents["BMEA"][0]   # lowest bead of bottom monolayer
        n_top        = mol_counts_in["TMEA"]
        n_bot        = mol_counts_in["BMEA"]
        # n_layers = number of distinct GRA resnums = one per layer
        n_layers     = mol_counts_in["GRA"]
        top_resname  = "TMEA"
        bot_resname  = "BMEA"

        print(f"  GRA z: {z_gra_bottom:.3f} – {z_gra_top:.3f} nm  "
              f"({n_layers} layers)")
        print(f"  TMEA tip z: {z_top_term:.3f} nm  ({n_top} chains)")
        print(f"  BMEA tip z: {z_bot_tip:.3f} nm  ({n_bot} chains)")

        nx, ny, _, _ = define_grid_dimensions(lx, ly, b)
        gra_beads_1layer, _, _ = build_sheet(b, nx, ny, z=0.0, resname=GRA_RESNAME)
        n_sheet    = len(gra_beads_1layer)
        gra_bonds  = write_gra_bonds(gra_beads_1layer, b, lx, ly)
        gra_angles = write_gra_angles(gra_bonds)

        base_beads = hair_beads

    # PATH B — build hair from scratch
    else:
        # Build Graphene Sheet
        nx, ny, lx, ly = define_grid_dimensions(args.lx, args.ly, b)
        print(f"  Target  : {args.lx:.4f} x {args.ly:.4f} nm")
        print(f"  Grid    : nx={nx}, ny={ny}")
        print(f"  Actual  : {lx:.5f} x {ly:.5f} nm")

        lz = args.lz

        gra_beads_1layer, lx, ly = build_sheet(b, nx, ny, z=0.0, resname=GRA_RESNAME)
        n_sheet    = len(gra_beads_1layer)
        gra_bonds  = write_gra_bonds(gra_beads_1layer, b, lx, ly)
        gra_angles = write_gra_angles(gra_bonds)

        n_layers       = args.n
        slab_thickness = (n_layers - 1) * INTERLAYER_DZ
        z_gra_bottom   = round(lz / 2.0 - slab_thickness / 2.0, 6)
        z_gra_top      = round(z_gra_bottom + slab_thickness,    6)
        print(f"\nStacking {n_layers}-layer graphene slab...")
        print(f"  z_bottom={z_gra_bottom:.4f} nm  z_top={z_gra_top:.4f} nm")

        gra_slab_beads = place_graphene(gra_beads_1layer, n_layers, z_gra_bottom,
                                        INTERLAYER_DZ)

        # Top MEA monolayer
        print("\nPlacing top lipid monolayer...")
        l_top, n_top, top_resname = place_monolayer(
            lx, ly, z_surface=z_gra_top + g_gap,
            grafting_distance=d_graft, flip_z=False,
        )
        z_top_term = z_gra_top + g_gap + (N_LIPID_BEADS - 1) * BOND_LJ
        print(f"  Ghost z={z_gra_top + g_gap:.4f}  tip z={z_top_term:.4f}  "
              f"chains={n_top}")

        # Bottom MEA monolayer
        print("\nPlacing bottom lipid monolayer...")
        l_bot, n_bot, bot_resname = place_monolayer(
            lx, ly, z_surface=z_gra_bottom - g_gap,
            grafting_distance=d_graft, flip_z=True,
            x_offset=d_graft / 2.0, y_offset=d_graft / 2.0,
        )
        z_bot_tip = z_gra_bottom - g_gap - (N_LIPID_BEADS - 1) * BOND_LJ
        print(f"  Ghost z={z_gra_bottom - g_gap:.4f}  tip z={z_bot_tip:.4f}  "
              f"chains={n_bot}")

        # Sanity-check box size
        R_LIST = 1.3  # nm — r_cut (1.2) + neighbour-list skin (0.1)
        assert lx > 2 * R_LIST, f"lx={lx:.3f} nm too small (need > {2*R_LIST} nm)"
        assert ly > 2 * R_LIST, f"ly={ly:.3f} nm too small"
        assert lz > 2 * R_LIST, f"lz={lz:.3f} nm too small"
        slab_z_span = (z_top_term - z_bot_tip) + 2 * s_gap
        assert lz > slab_z_span + 2 * R_LIST, (
            f"lz too small: need {slab_z_span + 2 * R_LIST:.2f} nm, "
            f"have {lz:.2f} nm"
        )

        base_beads = gra_slab_beads + l_top + l_bot

    # Solvent placement 
    z_solv_top = z_top_term + s_gap
    z_solv_bot = z_bot_tip  - s_gap

    odorant_beads = []
    n_odorants    = 0
    solv_xyz      = []   # accumulates placed centres for checking overlaps

    if odorant is not None:
        n_odorants_top, n_odorants_bot = compute_mols(lx, ly, lz)
        print(f"\nPlacing {odorant} odorant "
              f"({'droplet' if is_droplet else 'vapor'} mode)...")
        odorant_beads, n_odorants, solv_xyz = add_odorant(
                    odorant, lx, ly, lz,
                    n_odorants_top, n_odorants_bot,
                    z_top=z_solv_top, z_bot=z_solv_bot,
                    droplet=is_droplet, max_tries=tries)

    water_beads      = []
    antifreeze_beads = []
    n_water          = 0
    n_afwater        = 0

    if args.add_water:
        n_water_top, n_water_bot = compute_mols(lx, ly, lz)
        water_beads, n_water = add_water(
            lx, ly, lz, n_water_top, n_water_bot,
            z_top=z_solv_top, z_bot=z_solv_bot, T=298.0,
            existing=solv_xyz, max_tries=tries,
        )

        if args.antifreeze:
            rng_af     = np.random.default_rng()
            n_replace  = max(1, round(n_water * 0.10))
            af_indices = set(rng_af.choice(n_water, size=n_replace, replace=False).tolist())

            new_water_beads = []

            for i, bead in enumerate(water_beads):
                if i in af_indices:
                    af_bead             = dict(bead)
                    af_bead['resname']  = 'WF'
                    af_bead['atomname'] = 'WF'
                    af_bead['_mtype']   = 'BP4'
                    antifreeze_beads.append(af_bead)
                else:
                    new_water_beads.append(bead)

            water_beads = new_water_beads
            n_water     = len(water_beads)
            n_afwater   = len(antifreeze_beads)

    # Assemble all beads
    all_beads = base_beads + odorant_beads + water_beads + antifreeze_beads
    for k, a in enumerate(all_beads, 1):
        a['atomnum'] = k
    n_total = len(all_beads)

    # Write ITP files
    gra_itp = os.path.join(run_dir, "graphene.itp")
    write_gra_itp(gra_itp, GRA_RESNAME, gra_beads_1layer, gra_bonds, gra_angles, b)
    print(f"\nWrote: {gra_itp}  ({n_sheet} beads/layer, "
          f"{len(gra_bonds)} bonds, {len(gra_angles)} angles)")

    # Write system .gro
    gro_file = os.path.join(run_dir, "system.gro")
    title = (f"MARTINI hair | {run_dir} | {n_layers} GRA layers | "
             f"{n_top}+{n_bot} MEA chains | {n_water} W | "
             f"{n_odorants} {odorant or 'none'} | "
             f"{lx:.3f}x{ly:.3f}x{lz:.3f} nm")
    write_gro(gro_file, title, all_beads, (lx, ly, lz))
    print(f"Wrote: {gro_file}  ({n_total} beads total)")

    # Write forcefield.itp
    bead_types = collect_bead_types(all_beads)

    # Write system.top
    # Molecule order must be the same as the .gro assembly order above.
    mol_counts = {}
    mol_counts[GRA_RESNAME] = n_layers   # one GRA molecule per layer
    mol_counts[top_resname] = n_top      # TMEA chains
    mol_counts[bot_resname] = n_bot      # BMEA chains
    if n_odorants > 0 and odorant:
        mol_counts[odorant] = n_odorants
    if n_water > 0:
        mol_counts["W"] = n_water
    if n_afwater > 0:
        mol_counts["WF"] = n_afwater

    itp_files = ["graphene.itp", "lipidchain_top.itp", "lipidchain_bot.itp", "odorants.itp"]

    top_file = os.path.join(run_dir, "system.top")
    write_system_top(
        top_file, run_dir,
        mol_counts      = mol_counts,
        bead_types_used = bead_types,
        itp_files       = itp_files,
        box             = (lx, ly, lz),
        odorant_resname = odorant,
        is_droplet      = is_droplet,
        if_water        = args.add_water,
        if_antifreeze   = args.antifreeze
    )
    print(f"Wrote: {top_file}")

    # Summary
    print(f"\n{'='*60}")
    print(f"  Folder  : {run_dir}/")
    print(f"  Beads   : {n_total} total")
    print(f"  Box     : {lx:.3f} x {ly:.3f} x {lz:.3f} nm")
    print(f"  Graphene: {n_layers} layers  ({n_sheet} beads/layer)")
    print(f"  Lipids  : {n_top} (top) + {n_bot} (bot)")
    if not args.input_hair:
        print(f"  d_graft : {d_graft} nm")
    print(f"  Water   : {n_water} W beads")
    print(f"  Odorant : {n_odorants} {odorant or 'none'} "
          f"({'droplet' if is_droplet else 'vapor'})")
    print(f"\n  Files written:")
    written = [gro_file, top_file, gra_itp]
    for fname in written:
        print(f"    {fname}")
    print(f"{'='*60}\n")

    print("[ molecules ] block:")
    for rn, count in mol_counts.items():
        print(f"  {rn:<12s}  {count}")

    # sysinfo.json for downstream analysis
    sysinfo = {
        "run_dir":      run_dir,
        "gro_file":     "system.gro",
        "lx": lx, "ly": ly, "lz": lz,
        "n_layers":     n_layers,
        "n_sheet":      n_sheet,
        "n_top":        n_top,
        "n_bot":        n_bot,
        "top_resname":  top_resname,
        "bot_resname":  bot_resname,
        "z_gra_top":    float(z_gra_top),
        "z_gra_bottom": float(z_gra_bottom),
        "z_top_term":   float(z_top_term),
        "z_bot_tip":    float(z_bot_tip),
        "z_solv_top":   float(z_solv_top),
        "z_solv_bot":   float(z_solv_bot),
        "s_gap":        s_gap,
        "odorant":      odorant,
        "is_droplet":   is_droplet,
        "n_odorants":   n_odorants,
        "n_water":      n_water,
        "n_afwater":    n_afwater,
        "mol_counts":   mol_counts,
        "itp_files":    itp_files,
        "bead_types":   bead_types,
        "input_hair":   args.input_hair,
    }
    sysinfo_path = os.path.join(run_dir, "sysinfo.json")
    with open(sysinfo_path, "w") as f:
        json.dump(sysinfo, f, indent=2)
    print(f"Wrote: {sysinfo_path}")

    # Write GROMACS index file
    ndx_file = os.path.join(run_dir, "index.ndx")
    write_index_ndx(ndx_file, all_beads, odorant=odorant)

    # Copy mdp files to folder
    for i in range(3):
        shutil.copy(f"./inputs/runfiles/min{i+1}.mdp", os.path.join(run_dir, f"min{i+1}.mdp"))
    for i in range(2):
        shutil.copy(f"./inputs/runfiles/nvt{i+1}.mdp", os.path.join(run_dir, f"nvt{i+1}.mdp"))

    # Copy .itp files
    shutil.copy("./inputs/forcefields/odorants.itp", os.path.join(run_dir, "odorants.itp"))
    shutil.copy("./inputs/forcefields/lipidchain_bot.itp", os.path.join(run_dir, "lipidchain_bot.itp"))
    shutil.copy("./inputs/forcefields/lipidchain_top.itp", os.path.join(run_dir, "lipidchain_top.itp"))
    shutil.copy("./inputs/forcefields/forcefield.itp", os.path.join(run_dir, "forcefield.itp"))
    
if __name__ == "__main__":
    main()