"""
RDKit tools.
"""

from typing import List, Union, Optional

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import Descriptors


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
    Union[ROMol, List[int]]
        RDKit molecule with added hydrogen and conformers and conformers' IDs

    Notes
    -----
    The actual number of generated conformers might differ from :code:`numConfs`
    because pruning is performed in order to return dissimilar conformers.
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
    Clean solution provided by @ptosco on `RDKit Discussions #4520`_.

    .. _`RDKit Discussions #4520`: https://github.com/rdkit/rdkit/discussions/4520
    """
    if isinstance(confIds, int):
        confIds = [confIds]

    # Copy original molecule without conformers
    cmol = Chem.Mol(mol, True)

    # Add selected conformers to molecule
    for confId in confIds:
        cmol.AddConformer(mol.GetConformer(confId), assignId=True)

    return cmol


def smi2props(smi: str, properties: Optional[List[str]] = None):
    """
    Compute properties from SMILES.

    Parameters
    ----------
    smi: str
        SMILES
    props: List[str]
        List of properties to be computed

    Returns
    -------
    dict
        Dictionary of properties

    Notes
    -----
    RDKit properties are case sensitive.

    Some code and ideas borrowed from `Calculate RDKit descriptors with Dask`_

    .. _`Calculate RDKit descriptors with Dask`: https://gist.github.com/PatWalters/beb437da364a2c4bf8de724a2039b903
    """
    mol = Chem.MolFromSmiles(smi)

    res = {}
    if mol:
        for name, prop in Descriptors.descList:
            if properties is None or name in properties:
                res[name] = prop(mol)

        if properties is not None and len(res) != len(properties):
            raise RuntimeError(
                f"Not all required properties can be computed: {properties}"
            )
    else:
        raise ValueError(f"Molecule can be created from SMILES: {smi}")

    return res
