"""Isolate AT's fringe against a FIXED Xsuite reference.

compare() normally converts the same AT element it tracks, so toggling
FringeQuad moves both codes together (the converter derives edge_entry_model
from it).  Here the Xsuite line is pinned to the FringeQuad=1 element
("full" edge model) while only AT's flag varies, which says whether AT's
multipole_fringe is what Xsuite's "full" model actually expects.
"""
from xs_compare import compare, sector_bend, add_field_errors


def build(fq, b0=None):
    e = sector_bend(fringe_quad=fq)
    if b0 is not None:
        e = add_field_errors(e, rel=b0, orders=(0,), skew=False)
    return e


ref_full = build(1)                      # Xsuite side: edge model "full"
ref_full_b0 = build(1, b0=1e-4)

print("=== Xsuite pinned to FringeQuad=1 ('full'); AT flag varies ===")
print("--- no field error ---")
compare(build(1), ref_elem=ref_full, label="AT fq=1 (matched)  ")
compare(build(0), ref_elem=ref_full, label="AT fq=0 (mismatch) ")

print("--- with B0 = 1e-4 ---")
compare(build(1, b0=1e-4), ref_elem=ref_full_b0, label="AT fq=1 (matched)  ")
compare(build(0, b0=1e-4), ref_elem=ref_full_b0, label="AT fq=0 (mismatch) ")

print()
print("=== Xsuite pinned to FringeQuad=0 ('dipole-only'); AT flag varies ===")
ref_dip = build(0)
ref_dip_b0 = build(0, b0=1e-4)
compare(build(0), ref_elem=ref_dip, label="AT fq=0 (matched)  ")
compare(build(1), ref_elem=ref_dip, label="AT fq=1 (mismatch) ")
compare(build(0, b0=1e-4), ref_elem=ref_dip_b0, label="AT fq=0 +B0 matched")
