from typing import List, Union

from rdkit import Chem
from rdkit.Chem import rdDistGeom


def generate_conformers(mol, numConfs: int, randomSeed: int = 42):
    """
    Generate conformers using ETKDG.

    Parameters
    ----------
    mol:
        RDKit molecule
    numConfs: int
        Number of conformers to be generated
    randomSeed: int
        Seed for the random number generator

    Returns
    -------
    RDKit molecule with added hydrogen and conformers and conformers' IDs

    Notes
    -----
    The actual number of generated conformers might differ from :code:`numConfs` because
    pruning is performed in order to return dissimilar conformers.
    """
    # Define ETKDG parameter set
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = randomSeed

    # Add hydrogens
    molh = Chem.AddHs(mol)

    # Generate conformers
    cids = rdDistGeom.EmbedMultipleConfs(molh, numConfs=numConfs, params=ps)

    return molh, list(cids)


def mol_from_conformers(mol, confIds: Union[int, List[int]]):
    """
    Create molecule with a subset of conformers.

    Parameters
    ----------
    mol:
        RDKit molecule
    confIds: Union[int, List[int]]
        Conformer index or list of conformer indices

    Returns
    -------
    ROMol
        RDKit molecule with conformers specified by :code:`confIds`

    Notes
    -----
    Clean solution provided by @ptosco on _`RDKit Discussions #4520`.

    _`RDKit Discussions #4520` : https://github.com/rdkit/rdkit/discussions/4520
    """
    if isinstance(confIds, int):
        confIds = [confIds]

    # Copy original molecule without conformers
    cmol = Chem.Mol(mol, True)

    # Add selected conformers to molecule
    for confId in confIds:
        cmol.AddConformer(mol.GetConformer(confId), assignId=True)

    return cmol
