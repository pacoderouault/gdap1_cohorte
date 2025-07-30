#!/usr/bin/env python
from pymol import cmd
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple
import numpy as np

# MODEL = "AF-Q8TB36-F1-model_v4.pdb"  # pdb identifier or pdb file path
MODEL = "7AIA"  # pdb identifier or pdb file path
# MUTATION_FILE = "nguyen_mutations.tsv"  # mutation file, tsv format with 2 required columns (res and pos) and one optional columns for group.
MUTATION_FILE = "cohorte_mutations.tsv"  # mutation file, tsv format with 2 required columns (res and pos) and one optional columns for group.
MUTATION_AS = "ballAndStick"  # spheres, stick, ballAndStick
VISIBLE_CHAINS = ["AAA"]  # chains to draw
SHOW_LABELS = False


def setPretty():
    # Workspace settings
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("orthoscopic", 0)
    cmd.set("transparency", 0.5)
    cmd.set("dash_gap", 0)
    cmd.set("ray_trace_mode", 1)  # normal color + black outline
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", 0.1)
    cmd.space("cmyk")
    cmd.set("ray_trace_fog", 0)
    cmd.set("depth_cue", 0)


def BallnStick(arg1, stick_color="marine", sphere_color=None):
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85", "elem C and " + arg1)
    cmd.color("gray98", "elem H and " + arg1)
    cmd.color("slate", "elem N and " + arg1)
    cmd.set("stick_radius", 0.15, arg1)
    cmd.set("sphere_scale", 0.18, arg1)
    cmd.set("sphere_scale", 0.13, arg1 + " and elem H")
    if sphere_color:
        cmd.set("sphere_color", sphere_color)
    cmd.set("stick_color", stick_color, arg1)
    cmd.set("dash_gap", 0.01)
    cmd.set("dash_radius", 0.035)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)
    cmd.set("valence", 1, arg1)
    cmd.set("valence_size", 0.2, arg1)
    cmd.zoom(arg1)
    cmd.hide("labels")


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


def show_disulfide_bonds(selection="all", cutoff=2.2):
    """
    Détecte et affiche les ponts disulfure dans la sélection spécifiée.
    Affiche un bâtonnet pour chaque liaison détectée.
    """
    # Sélectionne tous les atomes de soufre des cystéines
    cmd.select("cys_sg", f"{selection} and resn CYS and name SG")
    # Récupère les coordonnées de chaque atome SG
    sg_atoms = cmd.get_model("cys_sg").atom

    count = 0
    for i, atom1 in enumerate(sg_atoms):
        for j, atom2 in enumerate(sg_atoms):
            if j <= i:
                continue  # évite les doublons
            dist = (
                (atom1.coord[0] - atom2.coord[0]) ** 2
                + (atom1.coord[1] - atom2.coord[1]) ** 2
                + (atom1.coord[2] - atom2.coord[2]) ** 2
            ) ** 0.5
            if dist < cutoff:
                bond_name = f"ssbond_{count}"
                cmd.distance(bond_name, f"index {atom1.index}", f"index {atom2.index}")
                cmd.set("dash_color", "yellow", bond_name)
                cmd.set("dash_width", 2.5, bond_name)
                count += 1
    print(f"{count} ponts disulfure affichés.")


def main():
    # Réinitialisation de la scène
    cmd.reinitialize("everything")
    load_model(MODEL)

    # nettoyange des objects
    hide_objects()

    setPretty()
    # cmd.bg_color("white")

    # visualisation des chains:
    for chain in cmd.get_chains("prot"):
        if chain in VISIBLE_CHAINS:
            cmd.show("cartoon", f"chain {chain}")
            cmd.show("surface", f"chain {chain}")
            cmd.set("transparency", 0.7)
            cmd.set("cartoon_transparency", 0.0)
            cmd.set("surface_color", "white")
            cmd.color("skyblue", f"cartoon and chain {chain}")
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
            # cmd.show("spheres", group + " and name CA")
            # cmd.set("sphere_scale", 1, group + " and name CA")
            if MUTATION_AS == "spheres":
                cmd.show("spheres", group)
                cmd.set("sphere_scale", 0.6, group)
                cmd.color(mut[2], group)
            elif MUTATION_AS == "sticks":
                cmd.show("sticks", group)
                cmd.color(mut[2], group)
            elif MUTATION_AS == "ballAndStick":
                BallnStick(group, stick_color=mut[2])

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

    show_disulfide_bonds()

    cmd.orient("prot")


main()
