#!/usr/bin/env python3
"""Beam geometry quoted in Sections 2.1 and 3.2 and Table 1 of the manuscript, derived from the constants of
freefall_dose_model_monte.py: beam half-angles and numerical aperture (Eq. 2), the input beam at the lens, the 1/e
beam widths at the stem entrance, stem exit and slide, and the fraction of the stem cross-section inside the 1/e
boundary (Eq. 5) at each plane.  Prints only; no files written."""
import numpy as np
from freefall_dose_model_monte import w0x, w0y, W1x, W1y, z1, f_lens, thx, thy, nx, ny, R_tube, z_top, stem_len, L, v_z, dz, dt, r_f

print(f"half-angles: theta_x = {thx*1e3:.1f} mrad, theta_y = {thy*1e3:.1f} mrad;  NA_x = {np.sin(thx):.3f}, NA_y = {np.sin(thy):.3f}")
print(f"input beam at the lens (similar triangles, z1 = {z1:.0f} mm -> f = {f_lens:.0f} mm): {2*W1x*f_lens/z1:.1f} x {2*W1y*f_lens/z1:.1f} mm")
g = np.arange(-R_tube, R_tube, 0.005); GX, GY = np.meshgrid(g, g, indexing="ij"); bore = (GX**2 + GY**2) <= R_tube**2
for label, z in (("stem entrance", z_top), ("stem exit", z_top + stem_len), ("slide", z_top + L)):
    wx, wy = w0x + thx * z, w0y + thy * z
    ellipse = ((GX / wx)**2 + (GY / wy)**2) <= 1
    contour = (np.abs(GX / wx)**nx + np.abs(GY / wy)**ny) <= 1
    print(f"{label:14s} z = {z:4.0f} mm: 1/e full width {2*wx:.2f} x {2*wy:.2f} mm; "
          f"1/e ellipse covers {(bore & ellipse).sum()/bore.sum()*100:.0f}% of the bore, super-Gaussian 1/e contour {(bore & contour).sum()/bore.sum()*100:.0f}%")
print(f"transit: v_z = {v_z:.0f} mm/s, dz = {dz:.2f} mm, dt = {dt*1e3:.1f} ms; fly radius r_f = {r_f:.2f} mm; "
      f"free fall from rest over {L:.0f} mm = {np.sqrt(2*L/1e3/9.81)*1e3:.0f} ms")
