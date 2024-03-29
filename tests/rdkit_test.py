import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors


from pyrmeli import rdkit as rmrdkit


@pytest.fixture(scope="module", params=[10, 100])
def numConfs(request):
    """
    Number of conformers to generate.
    """
    return request.param


@pytest.fixture(scope="module", params=["CCC", "CCCC"])
def mol(request):
    """
    RDKit molecule.
    """
    return Chem.MolFromSmiles(request.param)


def test_generate_conformers(mol, numConfs):
    cmol, cids = rmrdkit.generate_conformers(mol, numConfs=numConfs, randomSeed=42)

    assert cmol.GetNumHeavyAtoms() == mol.GetNumHeavyAtoms()
    assert cmol.GetNumConformers() > 0
    assert cmol.GetNumConformers() <= numConfs
    assert len(cids) == cmol.GetNumConformers()


def test_mol_from_conformers_subsample(mol, numConfs):
    cmol, cids = rmrdkit.generate_conformers(mol, numConfs=numConfs, randomSeed=42)

    # Select every other conformer
    ids = cids[::2]
    scmol = rmrdkit.mol_from_conformers(cmol, ids)

    assert scmol.GetNumHeavyAtoms() == cmol.GetNumHeavyAtoms()
    assert scmol.GetNumConformers() == len(ids)


def test_mol_from_conformers_one(mol, numConfs):
    cmol, cids = rmrdkit.generate_conformers(mol, numConfs=numConfs, randomSeed=42)

    # Select single conformer
    scmol = rmrdkit.mol_from_conformers(cmol, cids[0])

    assert scmol.GetNumHeavyAtoms() == cmol.GetNumHeavyAtoms()
    assert scmol.GetNumConformers() == 1


@pytest.mark.parametrize(
    "desc", [[d[0] for d in Descriptors.descList], ["qed", "MolWt", "HeavyAtomCount"]]
)
@pytest.mark.parametrize("smi", ["CCC", "CCCC"])
def test_smi2props(smi, desc):
    props = rmrdkit.smi2props(smi, desc)

    assert len(props) == len(desc)
    assert props["HeavyAtomCount"] == len(smi)


@pytest.mark.parametrize("smi", ["CCC", "CCCC"])
def test_smi2props_wrong_properties(smi):
    with pytest.raises(RuntimeError):
        rmrdkit.smi2props(smi, properties=["X", "Y"])


def test_smi2props_wrong_smiles():
    with pytest.raises(ValueError):
        rmrdkit.smi2props("ccc", properties=["MolWt", "HeavyAtomCount"])
