"""
Notes
-----
Tests on graphical functionality only check that the function can be called correctly
and that the output is of the correct type.
"""
import pytest

import py3Dmol
from rdkit import Chem

from pyrmeli import rdkit as rmrdkit
from pyrmeli import py3dmol as rmpy3dmol


@pytest.fixture(scope="module", params=["CCC", "CCCC"])
def rdmol(request):
    """
    RDKit molecule.
    """
    mol = Chem.MolFromSmiles(request.param)
    molh, _ = rmrdkit.generate_conformers(mol, numConfs=10)

    return molh


def test_show_mol(rdmol):
    view = rmpy3dmol.show_mol(rdmol)
    assert isinstance(view, py3Dmol.view)


@pytest.mark.parametrize("confIds", [[0], [0, 1]])
def test_show_mol(rdmol, confIds):
    view = rmpy3dmol.show_mol(rdmol, confIds=confIds)
    assert isinstance(view, py3Dmol.view)
