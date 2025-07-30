from pymol import cmd


def superposer_structures():
    cmd.reinitialize()

    # Charger la structure X-ray
    cmd.fetch("7AIA", name="7AIA", async_=0)
    # cmd.fetch("7ALM", name="7ALM", async_=0)
    # cmd.fetch("7YWD", name="7YWD", async_=0)

    # Charger le modèle AlphaFold
    cmd.load("AF-Q8TB36-F1-model_v4.pdb", "af")

    # supprimer la chain B de la structure xray
    cmd.remove("7AIA and chain BBB")
    # cmd.remove("7ALM and chain B")
    # cmd.remove("7YWD and chain B")
    # cmd.remove("7YWD and chain C")
    # cmd.remove("7YWD and chain D")

    # Afficher cartoon uniquement
    cmd.hide("everything")
    cmd.show("cartoon", "7AIA")
    # cmd.show("cartoon", "7ALM")
    # cmd.show("cartoon", "7YWD")
    cmd.show("cartoon", "af")

    # Colorer les structure
    cmd.color("cyan", "7AIA")
    # cmd.color("magenta", "7ALM")
    # cmd.color("white", "7YWD")
    cmd.color("yellow", "af")

    # Superposition des structures
    cmd.align("af and chain A and resi 23-302", "7AIA and chain AAA")
    # cmd.align("7ALM and chain A", "7AIA and chain AAA")
    # cmd.align("7YWD and chain A", "7AIA and chain AAA")

    # transparence
    cmd.set("cartoon_transparency", 0.3, "af")
    cmd.set("cartoon_transparency", 0.3, "7ALM")
    # cmd.set("cartoon_transparency", 0.3, "7YWD")

    # Zoom sur l’ensemble
    cmd.zoom("7AIA or af")


superposer_structures()
