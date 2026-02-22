from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from libnpet.core.pipeline import Stage
from libnpet.core.ribosome_types import RibosomeProfile
from libnpet.core.types import StageContext, ArtifactType

WATER_NAMES = {"HOH", "WAT", "DOD"}

class Stage80TunnelLining(Stage):
    key = "80_tunnel_lining"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "proximity_A": float(c.lining_proximity_A),
            "include_nonpolymers": bool(c.lining_include_nonpolymers),
            "include_waters": bool(c.lining_include_waters),
        }

    def run(self, ctx: StageContext) -> None:
        import gemmi

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        proximity_A = float(c.lining_proximity_A)

        # --- mesh vertices as KDTree ---
        mesh_path = ctx.inputs.get("tunnel_mesh_path")
        if not mesh_path or not Path(mesh_path).exists():
            raise ValueError(
                f"[{self.key}] tunnel_mesh_path missing from ctx; Stage70 must run first"
            )
        mesh_verts = np.asarray(pv.read(mesh_path).points, dtype=np.float32)
        tree = cKDTree(mesh_verts)

        # --- load structure and profile ---
        mmcif_path = ctx.require("mmcif_path")
        profile: RibosomeProfile = ctx.require("profile")

        from libnpet.core.structure_selection import tunnel_debris_chains
        excluded_chains: set[str] = set(tunnel_debris_chains(ctx.rcsb_id, profile))
        if excluded_chains:
            print(f"[{self.key}] excluding debris/nascent chains: {sorted(excluded_chains)}")

        st = gemmi.read_structure(mmcif_path)
        st.setup_entities()

        # build auth_asym_id -> profile metadata for polymers
        polymer_meta: dict[str, dict] = {}
        for p in profile.proteins:
            polymer_meta[p.auth_asym_id] = {
                "type": "protein",
                "nomenclature": [n.value for n in p.nomenclature],
                "description": p.rcsb_pdbx_description,
            }
        for r in profile.rnas:
            polymer_meta[r.auth_asym_id] = {
                "type": "rna",
                "nomenclature": [n.value for n in r.nomenclature],
                "description": r.rcsb_pdbx_description,
            }
        for o in profile.other_polymers:
            polymer_meta[o.auth_asym_id] = {
                "type": "other_polymer",
                "nomenclature": [n.value for n in o.nomenclature],
                "description": o.rcsb_pdbx_description,
            }

        # --- batch atom -> mesh distance query ---
        model = st[0]

        all_coords = []
        all_keys = []
        for chain in model:
            for res in chain:
                rkey = (chain.name, str(res.seqid), res.name)
                for atom in res:
                    p = atom.pos
                    all_coords.append([p.x, p.y, p.z])
                    all_keys.append(rkey)

        dists, _ = tree.query(all_coords, k=1)

        # min distance per (chain, seqid, resname)
        res_min_dist: dict[tuple, float] = {}
        for rkey, d in zip(all_keys, dists):
            cur = res_min_dist.get(rkey, float("inf"))
            if d < cur:
                res_min_dist[rkey] = float(d)

        # --- find qualifying chain IDs (any atom within threshold) ---
        qualifying_chain_ids: set[str] = set()
        for (cid, seqid_str, res_name), min_dist in res_min_dist.items():
            if cid in excluded_chains:
                continue
            if min_dist > proximity_A:
                continue
            if res_name in WATER_NAMES:
                if c.lining_include_waters:
                    qualifying_chain_ids.add(cid)
            elif cid in polymer_meta:
                qualifying_chain_ids.add(cid)
            else:
                if c.lining_include_nonpolymers:
                    qualifying_chain_ids.add(cid)

        # --- closest approach per qualifying chain ---
        chain_min_dist: dict[str, float] = {}
        for (cid, _, _), min_dist in res_min_dist.items():
            if cid in qualifying_chain_ids:
                cur = chain_min_dist.get(cid, float("inf"))
                if min_dist < cur:
                    chain_min_dist[cid] = min_dist

        # --- classify into report buckets ---
        lining_polymers: dict[str, dict] = {}
        lining_nonpolymers: dict[str, dict] = {}
        lining_waters: list[str] = []

        for cid in qualifying_chain_ids:
            entry = {
                "auth_asym_id": cid,
                "closest_atom_A": round(chain_min_dist[cid], 3),
            }
            if cid in polymer_meta:
                lining_polymers[cid] = {
                    **entry,
                    "type": polymer_meta[cid]["type"],
                    "nomenclature": polymer_meta[cid]["nomenclature"],
                    "description": polymer_meta[cid]["description"],
                }
            else:
                # check if this chain is purely waters
                chain_res_names = {
                    res_name for (c_, _, res_name) in res_min_dist if c_ == cid
                }
                if all(rn in WATER_NAMES for rn in chain_res_names):
                    lining_waters.append(cid)
                else:
                    lining_nonpolymers[cid] = entry

        # --- build report ---
        report = {
            "rcsb_id": ctx.rcsb_id,
            "mesh_path": str(mesh_path),
            "proximity_A": proximity_A,
            "options": {
                "include_nonpolymers": bool(c.lining_include_nonpolymers),
                "include_waters": bool(c.lining_include_waters),
            },
            "summary": {
                "n_polymer_chains": len(lining_polymers),
                "n_nonpolymer_chains": len(lining_nonpolymers),
                "n_water_chains": len(lining_waters),
            },
            "polymers": sorted(
                lining_polymers.values(), key=lambda x: x["auth_asym_id"]
            ),
            "nonpolymers": sorted(
                lining_nonpolymers.values(), key=lambda x: x["auth_asym_id"]
            ),
            "waters": sorted(lining_waters),
        }

        report_path = stage_dir / "tunnel_lining.json"
        report_path.write_text(json.dumps(report, indent=2))
        ctx.store.register_file(
            name="tunnel_lining_report",
            stage=self.key,
            type=ArtifactType.JSON,
            path=report_path,
        )

        # --- write mmcif with full qualifying chains ---
        st_lining = gemmi.Structure()
        st_lining.name = st.name
        st_lining.cell = st.cell
        st_lining.spacegroup_hm = st.spacegroup_hm

        model_out = gemmi.Model("1")
        for chain in model:
            if chain.name in qualifying_chain_ids:
                model_out.add_chain(chain)

        if len(model_out) > 0:
            st_lining.add_model(model_out)

        cif_path = stage_dir / "tunnel_lining.cif"
        st_lining.make_mmcif_document().write_file(str(cif_path))
        ctx.store.register_file(
            name="tunnel_lining_mmcif",
            stage=self.key,
            type=ArtifactType.MMCIF,
            path=cif_path,
        )

        # copy both to run root
        shutil.copy2(str(report_path), str(ctx.store.run_dir / "tunnel_lining.json"))
        shutil.copy2(str(cif_path), str(ctx.store.run_dir / "tunnel_lining.cif"))

        ctx.inputs["tunnel_lining_report"] = report

        print(
            f"[{self.key}] {len(lining_polymers)} polymer chains, "
            f"{len(lining_nonpolymers)} nonpolymer chains, "
            f"{len(lining_waters)} water chains within {proximity_A}A of mesh"
        )
