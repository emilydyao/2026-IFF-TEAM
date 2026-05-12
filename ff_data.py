import math

BOND_LJ       = 0.47    # nm, standard MARTINI C-C bond length (b0 in ITP must match placed spacing)
INTERLAYER_DZ = 0.47    # nm  graphite interlayer spacing
GRA_RESNAME   = "GRA"
MEA_RESNAME   = "MEA"

GRA_CC   = 0.2836                       # nm, graphene C-C (MARTINI)
GRA_A    = math.sqrt(84) / 4 * GRA_CC   # = 0.6498 nm, grafting column pitch
GRA_B    = math.sqrt(3) * GRA_A         # = 1.1255 nm, grafting row pitch
GRA_PHI  = math.atan(math.sqrt(3) / 2 / 4.5)  # = 10.893 degrees, commensurate tilt

LIPID_BEADS = [
    # (atomname, MARTINI-type, charge, mass)
    ("CGH", "CGH",  0.0, 24.0),   # 1 ghost / anchor
    ("CNH", "P5",   0.0, 72.0),   # 2 head (amine/carboxyl)
    ("C5B", "C5",   0.0, 72.0),   # 3 C5 base
    ("C1A", "C1",   0.0, 72.0),   # 4 chain
    ("C1B", "C1",   0.0, 72.0),   # 5 chain
    ("C1C", "C1",   0.0, 72.0),   # 6 chain
    ("C1D", "C1",   0.0, 72.0),   # 7 chain
    ("C1T", "C1",   0.0, 90.0),   # 8 terminal methyl
]

N_LIPID_BEADS = len(LIPID_BEADS)
LIPID_BONDS = [(i, i + 1) for i in range(1, N_LIPID_BEADS)]
LIPID_ANGLES = [(i, i + 1, i + 2) for i in range(1, N_LIPID_BEADS - 1)]

DMS_TEMP = [
    ("S01", "C2",  0.0, 72.0),   # dimethyl sulphide — single apolar bead
]

LIN_TEMP = [
    ("C1A", "P1",  0.0, 72.0),   # linalool head (hydroxyl)
    ("C2A", "C2",  0.0, 72.0),
    ("C3A", "C1",  0.0, 72.0),
    ("C4A", "C1",  0.0, 72.0),
    ("C5A", "C1",  0.0, 72.0),
]

LIN_BONDS  = [(1,3),(2,3),(3,4),(4,5)]         # P1,C2 both on C3A; C3A–C4A–C5A chain
LIN_ANGLES = [(1,3,4),(2,3,4),(3,4,5)]

SYR_TEMP = [
    ("O1A", "P1",  0.0, 72.0),   # syringaldehyde polar head
    ("O2A", "P2",  0.0, 72.0),
    ("R1A", "SC4", 0.0, 45.0),   # ring
    ("R2A", "SC4", 0.0, 45.0),   # ring
    ("C1A", "P2",  0.0, 72.0),   # aldehyde
]

SYR_BONDS  = [(1,3),(2,3),(3,4),(4,5)]         # P1–SC4(3)–SC4(4)–P2(5); P2(2) on SC4(3)
SYR_ANGLES = [(1,3,4),(2,3,4),(3,4,5)]

PHE_TEMP = [
    ("C1A", "P1",  0.0, 72.0),   # phenethyl alcohol head
    ("C2A", "C1",  0.0, 72.0),
    ("R1A", "SC4", 0.0, 45.0),
    ("R2A", "SC4", 0.0, 45.0),
]

PHE_BONDS  = [(1,3),(2,3),(3,4)]               # P1,C1 both on SC4(3); SC4(3)–SC4(4)
PHE_ANGLES = [(1,3,4),(2,3,4)]

ODORANT_TEMPLATES = {
    "DMS": DMS_TEMP,
    "LIN": LIN_TEMP,
    "SYR": SYR_TEMP,
    "PHE": PHE_TEMP,
}

ODORANT_BONDS = {
    "DMS": [],
    "LIN": LIN_BONDS,
    "SYR": SYR_BONDS,
    "PHE": PHE_BONDS,
}

ODORANT_ANGLES = {
    "DMS": [],
    "LIN": LIN_ANGLES,
    "SYR": SYR_ANGLES,
    "PHE": PHE_ANGLES,
}

# MARTINI v2.2 LJ interaction matrix.
# The values below are (epsilon_kJ, sigma_nm) for each level:
_MARTINI_LEVELS = {
    "O":   (5.6, 0.47),
    "I":   (5.0, 0.47),
    "II":  (4.5, 0.47),
    "III": (4.0, 0.47),
    "IV":  (3.5, 0.47),
    "V":   (3.1, 0.47),
    "VI":  (2.7, 0.47),
    "VII": (2.3, 0.47),
    "VIII":(2.0, 0.47),
    # Ghost anchor level: near-zero LJ so CGH beads don't interact
    "GHOST": (0.0001, 0.47),
}

# Only the bead types actually used in this hair model are listed:
#   C1, C5, P5, P4 (water)
# Ghost / anchor beads (CGH) use C1 parameters; they are excluded from
# non-bonded interactions via the [ exclusions ] block in their .itp.

_PAIR_LEVELS = {
    # ── apolar chain interactions ──
    ("C1",  "C1"):  "IV",
    ("C1",  "C2"):  "IV",    # C2 treated same as C1 for apolar–apolar
    ("C1",  "C5"):  "V",
    ("C2",  "C2"):  "IV",
    ("C2",  "C5"):  "V",
    ("C5",  "C5"):  "IV",
    # ── apolar vs polar (water-type) ──
    ("C1",  "P4"):  "VIII",
    ("C1",  "P5"):  "VIII",
    ("C1",  "P1"):  "VIII",
    ("C1",  "P2"):  "VIII",
    ("C2",  "P4"):  "VIII",
    ("C2",  "P5"):  "VIII",
    ("C2",  "P1"):  "VII",   # C2 slightly less repelled than C1 from polar
    ("C2",  "P2"):  "VII",
    ("C5",  "P4"):  "V",
    ("C5",  "P5"):  "V",
    ("C5",  "P1"):  "V",
    ("C5",  "P2"):  "V",
    # ── polar self and cross ──
    ("P1",  "P1"):  "II",
    ("P1",  "P2"):  "II",
    ("P1",  "P4"):  "O",
    ("P1",  "P5"):  "O",
    ("P2",  "P2"):  "II",
    ("P2",  "P4"):  "O",
    ("P2",  "P5"):  "O",
    ("P4",  "P4"):  "I",     # water self
    ("P4",  "P5"):  "O",
    ("P5",  "P5"):  "O",
    # ── ring beads (SC4, sigma=0.43 nm) vs everything ──
    ("C1",  "SC4"): "V",
    ("C2",  "SC4"): "V",
    ("C5",  "SC4"): "IV",
    ("P1",  "SC4"): "VI",
    ("P2",  "SC4"): "V",
    ("P4",  "SC4"): "VII",
    ("P5",  "SC4"): "VII",
    ("SC4", "SC4"): "IV",
    # ── CGH ghost anchor bead: near-zero LJ with everything ──
    # The ghost bead is a grafting anchor only; it must not exert
    # meaningful forces on any other bead.  The [ exclusions ] block
    # in the MEA .itp handles intra-chain exclusions; these nonbond_params
    # suppress all inter-chain and inter-molecule ghost interactions.
    ("CGH", "CGH"):  "GHOST",   # C < C — sorted ok
    ("C1",  "CGH"):  "GHOST",   # C1 < CG — sorted ok
    ("C2",  "CGH"):  "GHOST",   # C2 < CG — sorted ok
    ("C5",  "CGH"):  "GHOST",   # C5 < CG — sorted ok
    ("CGH", "P1"):   "GHOST",   # CGH < P1 — was ("P1","CGH") WRONG
    ("CGH", "P2"):   "GHOST",   # CGH < P2 — was ("P2","CGH") WRONG
    ("CGH", "P4"):   "GHOST",   # CGH < P4 — was ("P4","CGH") WRONG
    ("CGH", "P5"):   "GHOST",   # CGH < P5 — was ("P5","CGH") WRONG
    ("CGH", "SC4"):  "GHOST",   # CGH < SC4 — was ("SC4","CGH") WRONG
}

_SC4_SIGMA = 0.43  # nm


def _sorted_pair(a, b):
    return (a, b) if a <= b else (b, a)

def _get_epsilon_sigma(typeA, typeB):
    key   = _sorted_pair(typeA, typeB)
    level = _PAIR_LEVELS.get(key)
    if level is None:
        # strip S-prefix fallback (e.g. SC4 → C4, not ideal but safe)
        tA2 = typeA[1:] if typeA.startswith("S") else typeA
        tB2 = typeB[1:] if typeB.startswith("S") else typeB
        level = _PAIR_LEVELS.get(_sorted_pair(tA2, tB2), "IV")
    eps, sig = _MARTINI_LEVELS[level]
    # override sigma for SC4 ring pairs
    if "SC4" in (typeA, typeB):
        sig = _SC4_SIGMA
    return eps, sig

def write_forcefield_itp(filename, bead_types_used):
    _MARTINI_MASS = {
        "C1": 72.0, "C2": 72.0, "C5": 72.0,
        "P1": 72.0, "P2": 72.0, "P4": 72.0, "P5": 72.0,
        "SC4": 45.0,
        "CGH": 24.0,   # ghost anchor bead — same mass as in LIPID_BEADS
    }
    # All beads in this system are neutral
    _MARTINI_CHARGE = {t: 0.0 for t in _MARTINI_MASS}

    with open(filename, 'w') as f:
        f.write("; MARTINI v2.2 force-field parameters\n")
        f.write("; System: 18-MEA / graphene / non-polarisable water / odorants\n")
        f.write("; Generated by gen_sys.py — do not edit manually\n\n")

        f.write("[ defaults ]\n")
        f.write("; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n")
        f.write("  1       1          no         1.0      1.0\n\n")

        f.write("[ atomtypes ]\n")
        f.write("; name  mass    charge  ptype  sigma(nm)  epsilon(kJ/mol)\n")
        for t in bead_types_used:
            mass   = _MARTINI_MASS.get(t, 72.0)
            charge = _MARTINI_CHARGE.get(t, 0.0)
            eps, sig = _get_epsilon_sigma(t, t)
            f.write(f"  {t:<6s}  {mass:6.1f}  {charge:+.1f}  A  "
                    f"{sig:.4f}  {eps:.4f}\n")
        f.write("\n")

        f.write("[ nonbond_params ]\n")
        f.write("; i      j      func  sigma(nm)    epsilon(kJ/mol)  ; level\n")
        for i, tA in enumerate(bead_types_used):
            for tB in bead_types_used[i:]:
                eps, sig = _get_epsilon_sigma(tA, tB)
                key   = _sorted_pair(tA, tB)
                level = _PAIR_LEVELS.get(key, "IV")
                f.write(f"  {tA:<6s}  {tB:<6s}  1  "
                        f"{sig:.4f}  {eps:.4f}  ; {level}\n")
        f.write("\n")