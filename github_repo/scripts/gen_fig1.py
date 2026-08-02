#!/usr/bin/env python3
"""Generate clean STROBE Figure 1 SVG."""
lines = []
lines.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 960 580">')
lines.append('<defs>')
lines.append('  <marker id="aB" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto"><path d="M0,0 L7,2.5 L0,5 Z" fill="#5C6BC0"/></marker>')
lines.append('  <marker id="aG" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto"><path d="M0,0 L7,2.5 L0,5 Z" fill="#43A047"/></marker>')
lines.append('  <marker id="aO" markerWidth="7" markerHeight="5" refX="7" refY="2.5" orient="auto"><path d="M0,0 L7,2.5 L0,5 Z" fill="#EF6C00"/></marker>')
lines.append('</defs>')
lines.append('<rect width="960" height="580" fill="#FFFFFF"/>')

# Title
lines.append('<text x="480" y="26" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="13" font-weight="bold" fill="#212121">Figure 1. Study Flow Diagram</text>')

# ═══ LEFT: CHARLS ═══
Lx, Lw = 36, 426
lines.append(f'<rect x="{Lx}" y="38" width="{Lw}" height="26" rx="3" fill="#283593"/>')
lines.append(f'<text x="{Lx+Lw/2}" y="56" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="10" font-weight="bold" fill="#ffffff">CHARLS 2015  Derivation Cohort</text>')

ch = [
    (Lx+8, 72, Lw-16, 32, "Participants with blood sample (Wave 3)", "N = 13,420", "#E8EAF6", "#283593"),
    (Lx+8, 114, Lw-16, 24, "Excluded: age &lt; 45 (n = 945)", "", "#FFF8E1", "#E65100"),
    (Lx+8, 146, Lw-16, 24, "Excluded: baseline CVD (n = 129), cancer (n = 32), eGFR &lt; 15 (n = 24)", "", "#FFF8E1", "#E65100"),
    (Lx+8, 178, Lw-16, 24, "Excluded: missing 3+ core biomarkers (n = 38)", "", "#FFF8E1", "#E65100"),
    (Lx+8, 220, Lw-16, 36, "Final analytic sample", "N = 12,436", "#C5CAE9", "#1A237E"),
    (Lx+8, 276, Lw-16, 44, "852 all-cause deaths (6.9%)", "Median follow-up 5.4 years", "#F5F5F5", "#546E7A"),
]
for x, y, w, h, label, sub, fill, stroke in ch:
    lines.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="0.7"/>')
    if sub:
        lines.append(f'<text x="{x+w/2}" y="{y+15}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8.5" fill="#263238">{label}</text>')
        lines.append(f'<text x="{x+w/2}" y="{y+28}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8" font-weight="bold" fill="{stroke}">{sub}</text>')
    else:
        lines.append(f'<text x="{x+w/2}" y="{y+h/2+3}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8" fill="{stroke}">{label}</text>')

for ay in [104, 138, 170]:
    lines.append(f'<line x1="{Lx+Lw/2}" y1="{ay}" x2="{Lx+Lw/2}" y2="{ay+8}" stroke="#5C6BC0" stroke-width="1" marker-end="url(#aB)"/>')
lines.append(f'<line x1="{Lx+Lw/2}" y1="202" x2="{Lx+Lw/2}" y2="218" stroke="#5C6BC0" stroke-width="1" marker-end="url(#aB)"/>')
lines.append(f'<line x1="{Lx+Lw/2}" y1="256" x2="{Lx+Lw/2}" y2="274" stroke="#5C6BC0" stroke-width="1" marker-end="url(#aB)"/>')

# ═══ RIGHT: NHANES ═══
Rx, Rw = 498, 426
lines.append(f'<rect x="{Rx}" y="38" width="{Rw}" height="26" rx="3" fill="#2E7D32"/>')
lines.append(f'<text x="{Rx+Rw/2}" y="56" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="10" font-weight="bold" fill="#ffffff">NHANES 1999-2016  External Validation Cohort</text>')

nh = [
    (Rx+8, 72, Rw-16, 32, "Eligible participants, 8 continuous waves", "N = 47,810", "#E8F5E9", "#2E7D32"),
    (Rx+8, 114, Rw-16, 24, "Excluded: age &lt; 45 (n = 17,832)", "", "#FFF8E1", "#E65100"),
    (Rx+8, 146, Rw-16, 24, "Excluded: baseline CVD (n = 9,412), cancer (n = 2,704)", "", "#FFF8E1", "#E65100"),
    (Rx+8, 178, Rw-16, 24, "Excluded: implausible anthropometry (n = 258)", "", "#FFF8E1", "#E65100"),
    (Rx+8, 220, Rw-16, 36, "Final analytic sample", "N = 17,804", "#C8E6C9", "#1B5E20"),
    (Rx+8, 276, Rw-16, 44, "3,446 all-cause deaths (19.4%), CVD 1,012 (5.7%)", "NDI-linked, median follow-up 9.7 years", "#F5F5F5", "#546E7A"),
]
for x, y, w, h, label, sub, fill, stroke in nh:
    lines.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="{fill}" stroke="{stroke}" stroke-width="0.7"/>')
    if sub:
        lines.append(f'<text x="{x+w/2}" y="{y+15}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8.5" fill="#263238">{label}</text>')
        lines.append(f'<text x="{x+w/2}" y="{y+28}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8" font-weight="bold" fill="{stroke}">{sub}</text>')
    else:
        lines.append(f'<text x="{x+w/2}" y="{y+h/2+3}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8" fill="{stroke}">{label}</text>')

for ay in [104, 138, 170]:
    lines.append(f'<line x1="{Rx+Rw/2}" y1="{ay}" x2="{Rx+Rw/2}" y2="{ay+8}" stroke="#43A047" stroke-width="1" marker-end="url(#aG)"/>')
lines.append(f'<line x1="{Rx+Rw/2}" y1="202" x2="{Rx+Rw/2}" y2="218" stroke="#43A047" stroke-width="1" marker-end="url(#aG)"/>')
lines.append(f'<line x1="{Rx+Rw/2}" y1="256" x2="{Rx+Rw/2}" y2="274" stroke="#43A047" stroke-width="1" marker-end="url(#aG)"/>')

# ═══ BOTTOM: NHANES Wave I ═══
Bx, Bw = 210, 540
lines.append(f'<rect x="{Bx}" y="340" width="{Bw}" height="24" rx="3" fill="#E65100"/>')
lines.append(f'<text x="{Bx+Bw/2}" y="356" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="10" font-weight="bold" fill="#ffffff">NHANES 2015-2016 Wave I  CRP Sub-analysis (exploratory)</text>')

bot = [
    (Bx+8, 372, Bw-16, 26, "Wave I participants with CRP data  N = 2,349", "", "#FFF3E0", "#E65100"),
    (Bx+8, 406, Bw-16, 36, "107 all-cause deaths (4.6%), 18 CVD deaths (0.8%), median FU 3.9 y", "", "#F5F5F5", "#546E7A"),
]
for x, y, w, h, label, sub, fill, stroke in bot:
    lines.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="{fill}" stroke="{stroke}" stroke-width="0.7"/>')
    lines.append(f'<text x="{x+w/2}" y="{y+h/2+3}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8" fill="#333333">{label}</text>')

lines.append(f'<line x1="{Bx+Bw/2}" y1="364" x2="{Bx+Bw/2}" y2="370" stroke="#EF6C00" stroke-width="1" marker-end="url(#aO)"/>')

# Connector from NHANES final to Wave I
lines.append(f'<path d="M {Rx+Rw/2} 320 L {Rx+Rw/2} 332 L 480 332 L 480 338" fill="none" stroke="#9E9E9E" stroke-width="0.7" stroke-dasharray="4,3" marker-end="url(#aO)"/>')

# Footer
lines.append('<text x="36" y="480" font-family="Arial,Helvetica,sans-serif" font-size="7" fill="#9E9E9E">CVD = cardiovascular disease; eGFR = estimated glomerular filtration rate; CRP = C-reactive protein; NDI = National Death Index.</text>')
lines.append('<text x="36" y="494" font-family="Arial,Helvetica,sans-serif" font-size="7" fill="#9E9E9E">CHARLS = China Health and Retirement Longitudinal Study; NHANES = National Health and Nutrition Examination Survey.</text>')
lines.append('<text x="36" y="508" font-family="Arial,Helvetica,sans-serif" font-size="7" fill="#9E9E9E">Network-CMIN = Mahalanobis distance-based biomarker network dysregulation score computed from 6-10 routine clinical biomarkers.</text>')

lines.append('</svg>')

path = 'C:/Users/kkkk/.claude/projects/C--Users-kkkk/network_cmin/figures/fig1_strobe.svg'
with open(path, 'w') as f:
    f.write('\n'.join(lines))
print(f"SVG written to {path}")

# Validate
import xml.etree.ElementTree as ET
ET.parse(path)
print("XML valid")
