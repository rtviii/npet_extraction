from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from Bio.PDB.MMCIFParser import FastMMCIFParser

from libnpet.core.ribosome_types import ConstrictionInfo, PTCInfo, RibosomeProfile
from libnpet.core.config import SETTINGS

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{rcsb_id}.cif"


def _download_mmcif(rcsb_id: str) -> Path:
    import requests
    rcsb_id = rcsb_id.upper()
    dest = SETTINGS.npet2_root / "mmcif" / rcsb_id / f"{rcsb_id}.cif"
    if dest.exists():
        print(f"  [mmcif] using cached {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_DOWNLOAD_URL.format(rcsb_id=rcsb_id)
    print(f"  [mmcif] downloading {url} ...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    print(f"  [mmcif] saved to {dest}")
    return dest


def _flatten_assembly_map(raw_map: list) -> list[dict]:
    """
    The API returns assembly_map with instances nested under
    rcsb_polymer_entity_instance_container_identifiers. Our AssemblyInfo
    model expects them flat. This converts between the two.
    """
    result = []
    for asm in raw_map:
        flat: dict = {
            "rcsb_id": asm["rcsb_id"],
            "polymer_entity_instances": [],
            "nonpolymer_entity_instances": None,
        }
        for inst in asm.get("polymer_entity_instances", []):
            ids = inst["rcsb_polymer_entity_instance_container_identifiers"]
            flat["polymer_entity_instances"].append({
                "entity_id": ids["entity_id"],
                "auth_asym_id": ids["auth_asym_id"],
            })
        nonpolys = asm.get("nonpolymer_entity_instances")
        if nonpolys:
            flat["nonpolymer_entity_instances"] = []
            for inst in nonpolys:
                ids = inst["rcsb_nonpolymer_entity_instance_container_identifiers"]
                flat["nonpolymer_entity_instances"].append({
                    "entity_id": ids["entity_id"],
                    "auth_asym_id": ids["auth_asym_id"],
                    "auth_seq_id": ids.get("auth_seq_id", ""),
                })
        result.append(flat)
    return result


def _parse_profile(data: dict) -> RibosomeProfile:
    """
    Parse a RibosomeProfile from either a local profile JSON or a full
    RibosomeStructure response from the API. Handles the assembly_map
    nesting difference.
    """
    if "assembly_map" in data and data["assembly_map"]:
        first = data["assembly_map"][0]
        # detect API format by presence of nested container identifiers
        if (
            first.get("polymer_entity_instances")
            and "rcsb_polymer_entity_instance_container_identifiers"
            in first["polymer_entity_instances"][0]
        ):
            data = {**data, "assembly_map": _flatten_assembly_map(data["assembly_map"])}
    return RibosomeProfile.model_validate(data)


class FileStructureProvider:
    """
    Loads atoms from a local mmCIF file.
    Profile can come from a local JSON or the API.
    """

    def __init__(
        self,
        mmcif_path: str | Path,
        profile_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.mmcif_path = Path(mmcif_path)
        self.profile_path = Path(profile_path) if profile_path else None
        self.api_base = api_base or SETTINGS.riboxyz_api_base

        if not self.mmcif_path.exists():
            raise FileNotFoundError(f"mmCIF not found: {self.mmcif_path}")
        if self.profile_path and not self.profile_path.exists():
            raise FileNotFoundError(f"Profile not found: {self.profile_path}")

    def fingerprint(self, rcsb_id: str) -> str:
        return f"mmcif:{self.mmcif_path}"

    def _load_profile(self, rcsb_id: str) -> RibosomeProfile:
        if self.profile_path:
            return _parse_profile(json.loads(self.profile_path.read_text()))

        import requests
        url = f"{self.api_base}/structures/profile"
        resp = requests.get(url, params={"rcsb_id": rcsb_id.upper()}, timeout=30)
        resp.raise_for_status()
        return _parse_profile(resp.json())

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        parser = FastMMCIFParser(QUIET=True)
        structure = parser.get_structure(rcsb_id, str(self.mmcif_path))
        atoms = list(structure[0].get_atoms())

        if not atoms:
            raise ValueError(f"No atoms found in {self.mmcif_path}")

        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray([getattr(a, "element", "") or a.get_id()[0] for a in atoms])
        profile = self._load_profile(rcsb_id)

        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(self.mmcif_path),
            "profile": profile,
        }


class FileLandmarkProvider:
    """
    Loads PTC and constriction site from local JSON files or the API.
    Each is a separate file matching the API response schema exactly:
      PTC:         {"location": [x,y,z], "residues": [...]}
      Constriction: {"location": [x,y,z]}
    """

    def __init__(
        self,
        ptc_path: Optional[str | Path] = None,
        constriction_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.ptc_path = Path(ptc_path) if ptc_path else None
        self.constriction_path = Path(constriction_path) if constriction_path else None
        self.api_base = api_base or SETTINGS.riboxyz_api_base

        if self.ptc_path and not self.ptc_path.exists():
            raise FileNotFoundError(f"PTC file not found: {self.ptc_path}")
        if self.constriction_path and not self.constriction_path.exists():
            raise FileNotFoundError(f"Constriction file not found: {self.constriction_path}")

    def fingerprint(self, rcsb_id: str) -> str:
        if self.ptc_path:
            return f"landmarks_file:ptc={self.ptc_path},constr={self.constriction_path}"
        return f"landmarks_api:{self.api_base}"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        import requests
        rcsb_id = rcsb_id.upper()

        if self.ptc_path:
            ptc_info = PTCInfo.model_validate(json.loads(self.ptc_path.read_text()))
        else:
            resp = requests.get(
                f"{self.api_base}/loci/ptc",
                params={"rcsb_id": rcsb_id},
                timeout=30,
            )
            resp.raise_for_status()
            ptc_info = PTCInfo.model_validate(resp.json())

        if self.constriction_path:
            constr_info = ConstrictionInfo.model_validate(
                json.loads(self.constriction_path.read_text())
            )
        else:
            resp = requests.get(
                f"{self.api_base}/loci/constriction_site",
                params={"rcsb_id": rcsb_id},
                timeout=30,
            )
            resp.raise_for_status()
            constr_info = ConstrictionInfo.model_validate(resp.json())

        return {
            "ptc_xyz": np.array(ptc_info.location, dtype=np.float32),
            "constriction_xyz": np.array(constr_info.location, dtype=np.float32),
        }