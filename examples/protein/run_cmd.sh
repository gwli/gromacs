# steps from https://jerkwin.coding.me/GMX/GMXtut-0/#%E7%AC%AC%E4%B8%80%E6%AD%A5:%E8%8E%B7%E5%8F%96%E5%B9%B6%E5%A4%84%E7%90%86pdb%E6%96%87%E4%BB%B6
gmx pdb2gmx -ignh -ff amber99sb-ildn -f 1OMB_fws.pdb -o fws.gro -p fws.top -water tip3p
gmx editconf -f fws.gro -o fws-PBC.gro -bt dodecahedron -d 1.2
gmx grompp -f em-vac-pme.mdp -c fws-PBC.gro -p fws.top -o em-vac.tpr -maxwarn 99999999
gmx mdrun -v -deffnm em-vac
gmx solvate -cp em-vac.gro -cs spc216.gro -p fws.top -o fws-b4ion.gro
gmx grompp -f em-sol-pme.mdp -c fws-b4ion.gro -p fws.top -o ion.tpr -maxwarn 99999999
gmx genion -s ion.tpr -o fws-b4em.gro -neutral -conc 0.15 -p fws.top
gmx grompp -f em-sol-pme.mdp -c fws-b4em.gro -p fws.top -o em-sol.tpr
gmx mdrun -v -deffnm em-sol
gmx grompp -f nvt-pr-md.mdp -c em-sol.gro -p fws.top -o nvt-pr.tpr
gmx mdrun -v -deffnm nvt-pr
