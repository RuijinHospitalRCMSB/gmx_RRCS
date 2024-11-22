# -*- coding: utf-8 -*-

import sys
from MDAnalysis import Universe

sys.setdefaultencoding('utf-8')

md_traj = Universe('..\examples\md0.pdb', '..\examples\md0.xtc')
# print(md_traj)
# print(md_traj.trajectory[0].positions.shape)
# print(md_traj.trajectory[0].positions[100])
# print(md_traj.trajectory[1].positions[100])
# print(md_traj.trajectory[5].positions[100])
# print(md_traj.trajectory[5].positions[200])
# help(md_traj.trajectory[0])

prefixs = ['', '1', '2', '3']
h_atoms = ""
for prefix in prefixs:
    h_atoms += f"and not name {prefix}H* "

print(h_atoms)
protein = md_traj.select_atoms("protein")
atoms = protein.select_atoms(f"record_type ATOM and resid 103 and segindex 0 {h_atoms}") # not name *HH* and not name *HD*")
print([a.name for a in atoms])