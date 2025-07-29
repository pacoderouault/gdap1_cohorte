#!/usr/bin/env python
from pymol import cmd
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple
import numpy as np

MODEL = "/mnt/Data/CHU_services/biochimie_genetique_moleculaire/Recherche/GDAP1/Camille_S/AF-Q8TB36-F1-model_v4.pdb"  # pdb identifier or pdb file path
# MODEL = "7AIA"  # pdb identifier or pdb file path
MUTATION_FILE = "/mnt/Data/CHU_services/biochimie_genetique_moleculaire/Recherche/GDAP1/Camille_S/mutations_GDAP1.tsv"  # mutation file, tsv format with 2 required columns (res and pos) and one optional columns for group.
VISIBLE_CHAINS = ["A"]  # chains to draw
SHOW_LABELS = False


def load_model(model, object_name="prot"):
    if Path(model).exists():
        cmd.load(MODEL, object_name)
        print(f"Chargement local réussi : {MODEL}")
    else:
        cmd.fetch(MODEL, name=object_name, async_=0)
        print(f"Fetch réussi : {MODEL}")


def hide_objects(object_name="prot"):
    # hide chains
    chains = cmd.get_chains(object_name)
    for chain in chains:
        if chain not in VISIBLE_CHAINS:
            cmd.remove(f"chain {chain}")
    # hide solvent
    cmd.remove("solvent")
    cmd.remove("hydrogens")
    cmd.remove("hetatm")


def parse_mutation_file(file: Path):
    mutations = {}
    with open(file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.rstrip().split("\t")
                if len(cols) < 2:
                    raise ValueError("Invalid mutation file format.")
                if not cols[1].isdigit():
                    raise ValueError("The columns position must contains only digits.")

                res = cols[0]
                pos = int(cols[1])
                if len(cols) > 2:
                    group = cols[2]
                else:
                    group = None
                if len(cols) > 3:
                    color = cols[3]
                else:
                    color = None
                if not group in mutations:
                    mutations[group] = []
                mutations[group].append((res, pos, color))
    return mutations


def main():
    # Réinitialisation de la scène
    cmd.reinitialize("everything")
    load_model(MODEL)

    # nettoyange des objects
    hide_objects()

    cmd.bg_color("white")

    # visualisation des chains:
    for chain in cmd.get_chains("prot"):
        if chain in VISIBLE_CHAINS:
            cmd.show("cartoon", f"chain {chain}")
            cmd.show("surface", f"chain {chain}")
            cmd.set("transparency", 0.5)
            cmd.set("cartoon_transparency", 0.0)
            cmd.set("surface_color", "white")
            cmd.color("bluewhite", f"cartoon and chain {chain}")
            # color bluewhite, cartoon and chain A

    if MUTATION_FILE:
        mutations = parse_mutation_file(MUTATION_FILE)
        # print(mutations)

        for group, muts in mutations.items():

            # créer la sélection de toutes les mutations
            resi = "resi "
            for mut in muts:
                resi += f"+{mut[1]}"
            cmd.select(group, f"chain {chain} and {resi}")

            # Colorer la surface de ces résidus
            cmd.set("surface_color", mut[2], group)

            # Montrer sphères et colorer les atomes mutés
            cmd.show("spheres", group + " and name CA")
            cmd.set("sphere_scale", 1, group + " and name CA")
            cmd.color(mut[2], group)

            if SHOW_LABELS:
                label = '"%s%s" % (resn, resi)'
                cmd.label(selection=group + " and name CA", expression=label)

        if SHOW_LABELS:
            cmd.center("prot")
            cmd.zoom("prot", buffer=30)
            cmd.clip("slab", 300)

            cmd.set("label_bg_color", "pink")
            cmd.set("label_color", "black")
            cmd.set("label_outline_color", "black")
            cmd.set("label_connector", "on")
            cmd.set("label_connector_width", 1.0)

            xyz_center = np.asarray(cmd.get_position())
            xyz_mut = cmd.get_coords("prot", 1)
            xyz_diff = xyz_mut[0] - xyz_center
            new_x = 70 * (
                xyz_diff[0] / (abs(xyz_diff[0]) + abs(xyz_diff[1]) + abs(xyz_diff[2]))
            )
            new_y = 70 * (
                xyz_diff[1] / (abs(xyz_diff[0]) + abs(xyz_diff[1]) + abs(xyz_diff[2]))
            )
            new_z = 70 * (
                xyz_diff[2] / (abs(xyz_diff[0]) + abs(xyz_diff[1]) + abs(xyz_diff[2]))
            )
            cmd.set("label_placement_offset", (new_x, new_y, new_z), "prot")


main()
