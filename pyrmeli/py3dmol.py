import py3Dmol
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import rdMolAlign

from typing import Union, List


def show_mol(
    mol,
    confIds: Union[int, List[int]] = [0],
    alignConfs=True,
    view=None,
    color="lightgrey",
):
    """
    Add conformers of a molecule to py3Dmol view.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol
        Molecule to be visualized
    confIds: Union[int, List[int]]=[0]
        Conformer IDs to be visualized.
    view:
        py3Dmol.view. If :code:`None`, a new one will be created.

    Returns
    -------
    py3Dmol.view
        Upadated view with conformers.
    """
    if isinstance(confIds, int):
        confIds = [confIds]

    # Process supported molecules
    if isinstance(mol, rdchem.Mol):
        if alignConfs:
            rdMolAlign.AlignMolConformers(mol, confIds=confIds)
        smols = [Chem.MolToMolBlock(mol, confId=id) for id in confIds]
    else:
        # TODO: Support OpenBabel molecule
        raise ArgumentError(f"Molecule of type {type(mol)} is not supported.")

    # Create py3Dmol view
    if view is None:
        # TODO: Allow to define view size
        view = py3Dmol.view()

    for idx, smol in enumerate(smols):
        view.addModel(smol, "mol")

    # TODO: Change only color of added models (if view is not None)
    view.setStyle({"stick": {"colorscheme": f"{color}Carbon"}})
    view.zoomTo()

    return view
