So i have this pipeline for processing the ribosome npet tunnel and it seems to work pretty well locally (produces desirable results), but for whatever reason it keeps failing silently at stage 2 when i try to pacakge it into docker container. 

Can we investigate why? This is the last step for me to ship this shit so i would really like to crack it...

Let me show you the code. First of all im even having trouble with mounting this shit correctly.. 

```
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  npet2 run 5afi                                                                                                                    [master]
npet2: processing 1 structure(s), workers=1
  [mmcif] using cached /Users/rtviii/dev/riboxyz/NPET2/mmcif/5AFI/5AFI.cif

============================================================
  npet2 | 5AFI | run 20260211_20260221_192150_55d3ac7a92f87a78
  stages: 8 | config: RunConfig
============================================================

--- [1/8] 00_inputs ---------------------------------------
  [00_inputs] done in 2.17s

--- [2/8] 10_landmarks ------------------------------------
  [10_landmarks] PTC-Constriction distance=26.6A, cylinder z=[-20.0, 120.0]A
  [10_landmarks] done in 0.01s

--- [3/8] 20_exterior_shell -------------------------------
  [20_exterior_shell] d3d_alpha=200, d3d_tol=10, d3d_offset=3, kdtree_radius=40
  [20_exterior_shell] max_nn=60, tangent_k=20, poisson_depth=6, poisson_ptweight=4, fill_holes=2000
  [quick_surface_points] subsampled 147,762 -> 80,000 points
  [quick_surface_points] surface: 296 pts from alpha shape
Open3D Poisson Reconstruction: /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78/stage/20_exterior_shell/alpha_normals.ply -> /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78/stage/20_exterior_shell/alpha_shell.ply
>>Wrote /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78/stage/20_exterior_shell/alpha_shell.ply
  [20_exterior_shell] done in 8.82s

--- [4/8] 30_region_atoms ---------------------------------
  [30_region_atoms] radius_A=35, height_A=120
[30_region_atoms] seed_atoms=15,851 occ_atoms=15,851 occ_chains=56
[30_region_atoms] excluded: v(tRNA (auto-detected)), w(tRNA (auto-detected)), y(tRNA (auto-detected))
  [30_region_atoms] done in 7.91s

--- [5/8] 40_empty_space ----------------------------------
  [40_empty_space] radius_A=35, height_A=120
  [40_empty_space] grid_levels[0]: name=level_0, voxel_size_A=1, backend=legacy_kdtree, atom_radius_mode=uniform, uniform_atom_radius_A=2
/Users/rtviii/dev/npet_extraction/libnpet/stages/legacy_minimal.py:582: PyVistaDeprecationWarning: This filter is deprecated. Use `select_interior_points` instead.
  sel = pts_poly.select_enclosed_points(shell, check_surface=watertight)
  [40_empty_space] done in 3.10s

--- [6/8] 50_clustering -----------------------------------
  [50_clustering] coarse_eps_A=5.5, coarse_min_samples=600, refine_eps_A=3.5
  [50_clustering] refine_min_samples=175, mesh_enable=True, cluster_selection=axial_proximity
[50_clustering] empty_points n=204,924
Running DBSCAN on 204924 points. eps=5.5, min_samples=600, distance_metric=euclidean
[50_clustering] coarse DBSCAN: 4.82s, 10 clusters
  [cluster_select] picked cluster 0 (n=85,399, dist_to_constriction=0.4A)
Running DBSCAN on 85399 points. eps=3.5, min_samples=175, distance_metric=euclidean
[50_clustering] refine DBSCAN: 1.19s, 4 clusters
  [cluster_select] picked cluster 1 (n=69,754, dist_to_constriction=0.4A)
[50_clustering] winner: coarse=85,399 -> refine=69,754
[50_clustering] generating mesh for level_0...
[50_clustering]   MC mesh: 0.37s, 28,856 pts, 57,712 faces, watertight=True
[50_clustering] mesh saved: /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78/stage/50_clustering/mesh_level_0.ply
  [50_clustering] done in 6.69s

--- [7/8] 55_grid_refine ----------------------------------
  [55_grid_refine] voxel_size_A=0.5, roi_pad_A=10, atom_radius_A=2, keep_within_A=6, occ_close_iters=0, void_open_iters=1, forbid_roi_boundary=True, coarse_eps_A=3
  [55_grid_refine] coarse_min_samples=30, refine_eps_A=3, refine_min_samples=20, dbscan_max_points=0, dbscan_seed=0, mesh_enable=True, mesh_poisson_depth=8, mesh_poisson_ptweight=0.5
[55_grid_refine] voxel=0.5Å, ROI pad=10.0Å, atom_r=2.0Å
[55_grid_refine] selected 14,383 atoms near ROI (ALL atoms, prevents interference)
/Users/rtviii/dev/npet_extraction/libnpet/stages/grid_refine.py:207: PyVistaDeprecationWarning: This filter is deprecated. Use `select_interior_points` instead.
  sel = pv.PolyData(empty_pts_c0).select_enclosed_points(shell_c0, check_surface=watertight)
[55_grid_refine] boundary points: 137,228
[55_grid_refine] coarse DBSCAN: 0.44s, 1 clusters
[55_grid_refine] refine DBSCAN: 0.43s, 1 clusters
[55_grid_refine] creating selected void mask for voxel fallback...
[55_grid_refine]   selected mask: 784,921 voxels
[55_grid_refine] DBSCAN refined-grid: coarse=0 → refine=0, final_pts=137,178
[55_grid_refine] generating mesh for level_1...
[55_grid_refine]   MC mesh: 2.52s, 206,987 pts, 414,110 faces, ...
[55_grid_refine] mesh saved: /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78/stage/55_grid_refine/mesh_level_1.ply
  [55_grid_refine] done in 16.08s

--- [8/8] 70_mesh_validate --------------------------------
[70_mesh_validate] using level_1 mesh (0.5A grid)
[70_mesh_validate] final mesh: {'n_points': 206987, 'n_faces': 414110, 'open_edges': 0, 'is_manifold': True, 'bounds': [151.07286071777344, 224.272705078125, 165.71090698242188, 262.77581787109375, 68.97809600830078, 160.3324737548828]}, watertight=True
[70_mesh_validate]   copied level_0 mesh (1.0A grid) + pre-smooth
[70_mesh_validate]   copied level_1 mesh (0.5A grid) + pre-smooth
  [70_mesh_validate] done in 2.08s

============================================================
  npet2 | 5AFI | completed in 46.89s
  run_dir: /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78
============================================================

  5AFI: success -> /Users/rtviii/dev/riboxyz/NPET2/runs/5AFI/20260211_20260221_192150_55d3ac7a92f87a78

npet2: 1 succeeded, 0 failed out of 1
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  docker run \                                                                                                                      [master]
  -v $(pwd)/npet2_data:/data \
  -v /path/to/your/structures:/structures:ro \
  npet2 run 5NWY \
    --mmcif ./structures/5NWY.cif \
    --profile ./structures/5NWY.json \
    --ptc ./structures/5NWY_PTC.json \
    --constriction ./structures/5NWY_CONSTRICTION_SITE.json
docker: Error response from daemon: Mounts denied:
The path /path/to/your/structures is not shared from the host and is not known to Docker.
You can configure shared paths from Docker -> Preferences... -> Resources -> File Sharing.
See https://docs.docker.com/desktop/mac for more info.
ERRO[0000] error waiting for container:
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  docker run \                                                                                                                      [master]
  -v $(pwd)/npet2_data:/data \
  -v ./structures:/structures:ro \
  npet2 run 5NWY \
    --mmcif ./structures/5NWY.cif \
    --profile ./structures/5NWY.json \
    --ptc ./structures/5NWY_PTC.json \
    --constriction ./structures/5NWY_CONSTRICTION_SITE.json
npet2: processing 1 structure(s), workers=1
Error: --mmcif structures/5NWY.cif does not exist
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]                                                                                                                                    [master]
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  l                                                                                                                                 [master]
total 56
drwxr-xr-x  16 rtviii  staff   512B Feb 21 13:17 .
drwxr-xr-x@ 28 rtviii  staff   896B Feb 19 14:36 ..
-rw-r--r--   1 rtviii  staff    44B Feb 19 15:34 .dockerignore
drwxr-xr-x  14 rtviii  staff   448B Feb 20 14:21 .git
-rw-r--r--   1 rtviii  staff    66B Feb 20 11:28 .gitignore
drwxr-xr-x   8 rtviii  staff   256B Feb 19 14:44 .venv
-rw-r--r--   1 rtviii  staff   295B Feb 19 15:32 docker-compose.yml
-rw-r--r--   1 rtviii  staff   704B Feb 21 13:17 Dockerfile
drwxr-xr-x  10 rtviii  staff   320B Feb 20 11:38 libnpet
drwxr-xr-x   8 rtviii  staff   256B Feb 21 13:18 npet_extraction.egg-info
drwxr-xr-x   5 rtviii  staff   160B Feb 19 16:32 npet2_data
-rw-r--r--   1 rtviii  staff   458B Feb 19 16:42 pyproject.toml
-rw-r--r--   1 rtviii  staff   382B Feb 21 19:21 q_eliminate_poissonrecon.md
-rw-r--r--   1 rtviii  staff     0B Feb 21 17:51 q_package_docker.md
-rw-r--r--   1 rtviii  staff   705B Feb 21 13:25 README.md
drwxr-xr-x  21 rtviii  staff   672B Feb 19 16:41 structures
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  cd structures                                                                                                                     [master]
(.venv) ᢹ saeta.rtviii[ npet_extraction/structures ]  l                                                                                                                          [master]
total 75056
drwxr-xr-x  21 rtviii  staff   672B Feb 19 16:41 .
drwxr-xr-x  16 rtviii  staff   512B Feb 21 13:17 ..
-rw-r--r--   1 rtviii  staff   319K Feb 19 16:41 5NWY_ALPHA_SHAPE_ascii.ply
-rw-r--r--   1 rtviii  staff   177K Feb 19 16:41 5NWY_ALPHA_SHAPE.ply
-rw-r--r--   1 rtviii  staff    72B Feb 19 16:41 5NWY_CONSTRICTION_SITE.json
-rw-r--r--   1 rtviii  staff   109K Feb 19 16:41 5NWY_convex_hull.npy
-rw-r--r--   1 rtviii  staff   218K Feb 19 16:41 5NWY_normal_estimated_surf.ply
-rw-r--r--   1 rtviii  staff   1.1M Feb 19 16:41 5NWY_NPET_MESH_ascii.ply
-rw-r--r--   1 rtviii  staff   544K Feb 19 16:41 5NWY_NPET_MESH.ply
-rw-r--r--   1 rtviii  staff   606K Feb 19 16:41 5NWY_poisson_recon_ascii.ply
-rw-r--r--   1 rtviii  staff   301K Feb 19 16:41 5NWY_poisson_recon.ply
-rw-r--r--   1 rtviii  staff   1.7K Feb 19 16:41 5NWY_PTC.json
-rw-r--r--   1 rtviii  staff    12M Feb 19 16:41 5NWY_spheres_expanded_pointset.npy
-rw-r--r--   1 rtviii  staff   5.2M Feb 19 16:41 5NWY_tunnel_atoms_bbox.json
-rw-r-xr-x   1 rtviii  staff    15M Feb 19 16:41 5NWY.cif
-rw-r-xr-x   1 rtviii  staff    66K Feb 19 16:41 5NWY.json
-rw-r--r--   1 rtviii  staff   535K Feb 19 16:41 5NWY.png
drwxr-xr-x  12 rtviii  staff   384B Feb 19 16:41 artifacts
-rw-r--r--   1 rtviii  staff   219K Feb 19 16:41 classification_report_5NWY.json
-rwx------   1 rtviii  staff    42K Feb 19 16:41 tunnel_5NWY.csv
drwxr-xr-x   9 rtviii  staff   288B Feb 19 16:41 TUNNELS
(.venv) ᢹ saeta.rtviii[ npet_extraction/structures ]                                                                                                                             [master]
(.venv) ᢹ saeta.rtviii[ npet_extraction/structures ]  l                                                                                                                          [master]
total 75056
drwxr-xr-x  21 rtviii  staff   672B Feb 19 16:41 .
drwxr-xr-x  16 rtviii  staff   512B Feb 21 13:17 ..
-rw-r--r--   1 rtviii  staff   319K Feb 19 16:41 5NWY_ALPHA_SHAPE_ascii.ply
-rw-r--r--   1 rtviii  staff   177K Feb 19 16:41 5NWY_ALPHA_SHAPE.ply
-rw-r--r--   1 rtviii  staff    72B Feb 19 16:41 5NWY_CONSTRICTION_SITE.json
-rw-r--r--   1 rtviii  staff   109K Feb 19 16:41 5NWY_convex_hull.npy
-rw-r--r--   1 rtviii  staff   218K Feb 19 16:41 5NWY_normal_estimated_surf.ply
-rw-r--r--   1 rtviii  staff   1.1M Feb 19 16:41 5NWY_NPET_MESH_ascii.ply
-rw-r--r--   1 rtviii  staff   544K Feb 19 16:41 5NWY_NPET_MESH.ply
-rw-r--r--   1 rtviii  staff   606K Feb 19 16:41 5NWY_poisson_recon_ascii.ply
-rw-r--r--   1 rtviii  staff   301K Feb 19 16:41 5NWY_poisson_recon.ply
-rw-r--r--   1 rtviii  staff   1.7K Feb 19 16:41 5NWY_PTC.json
-rw-r--r--   1 rtviii  staff    12M Feb 19 16:41 5NWY_spheres_expanded_pointset.npy
-rw-r--r--   1 rtviii  staff   5.2M Feb 19 16:41 5NWY_tunnel_atoms_bbox.json
-rw-r-xr-x   1 rtviii  staff    15M Feb 19 16:41 5NWY.cif
-rw-r-xr-x   1 rtviii  staff    66K Feb 19 16:41 5NWY.json
-rw-r--r--   1 rtviii  staff   535K Feb 19 16:41 5NWY.png
drwxr-xr-x  12 rtviii  staff   384B Feb 19 16:41 artifacts
-rw-r--r--   1 rtviii  staff   219K Feb 19 16:41 classification_report_5NWY.json
-rwx------   1 rtviii  staff    42K Feb 19 16:41 tunnel_5NWY.csv
drwxr-xr-x   9 rtviii  staff   288B Feb 19 16:41 TUNNELS
(.venv) ᢹ saeta.rtviii[ npet_extraction/structures ]  cd ..                                                                                                                      [master]
(.venv) ᢹ saeta.rtviii[ dev/npet_extraction ]  docker run \                                                                                                                      [master]
  -v $(pwd)/npet2_data:/data \
  -v ./structures:/structures:ro \
  npet2 run 5NWY \
    --mmcif ./structures/5NWY.cif \
    --profile ./structures/5NWY.json \
    --ptc ./structures/5NWY_PTC.json \
    --constriction ./structures/5NWY_CONSTRICTION_SITE.json
npet2: processing 1 structure(s), workers=1
Error: --mmcif structures/5NWY.cif does not exist

```


Let me know if you want to see any other files...

```
Dockerfile
```
FROM ubuntu:22.04

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3-pip \
    libgomp1 \
    libgl1 \
    libglib2.0-0 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /npet2

COPY pyproject.toml .
COPY libnpet/ libnpet/

RUN pip3 install --no-cache-dir --upgrade pip setuptools && \
    pip3 install --no-cache-dir .

ENV NPET2_ROOT=/data
ENV NPET2_RUNS_ROOT=/data/runs
ENV NPET2_CACHE_ROOT=/data/cache
ENV NPET2_RIBOXYZ_API_URL=https://api.ribosome.xyz

VOLUME ["/data"]

ENTRYPOINT ["npet2"]
CMD ["--help"]
```

docker-compose.yml
```yml
services:
  npet2:
    build: .
    image: npet2:latest
    volumes:
      - ${NPET2_DATA:-./npet2_data}:/data
    environment:
      - NPET2_RIBOXYZ_API_URL=${NPET2_RIBOXYZ_API_URL:-http://host.docker.internal:8000}
    # pass arbitrary CLI args after --
    # docker compose run npet2 run 5NWY
```

libnpet/run.py
```py
# npet2/run.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

from libnpet.adapters.standalone_providers import FileStructureProvider, FileLandmarkProvider
from libnpet.core.config import RunConfig, SETTINGS
from libnpet.core.manifest import RunManifest
from libnpet.core.run_id import compute_run_id
from libnpet.core.store import LocalRunStore
from libnpet.core.types import StageContext
from libnpet.core.pipeline import Pipeline
from libnpet.core.cache import LocalStageCache

from libnpet.stages.bootstrap import Stage00Inputs, Stage10Landmarks
from libnpet.stages.grid_refine import Stage55GridRefine
from libnpet.stages.legacy_minimal import (
    Stage20ExteriorShell,
    Stage30RegionAtoms,
    Stage40EmptySpace,
    Stage50Clustering,
    Stage70MeshValidate,
)


def _pipeline_version() -> str:
    return "npet2-dev"


def run_npet2(
    rcsb_id: str,
    config: Optional[RunConfig] = None,
    *,
    structure_provider=None,
    landmark_provider=None,
) -> StageContext:
    rcsb_id = rcsb_id.upper()
    config = config or RunConfig()

    if structure_provider is None or landmark_provider is None:
        raise ValueError(
            "Both structure_provider and landmark_provider are required. "
            "Use FileStructureProvider / FileLandmarkProvider for standalone mode."
        )

    config_resolved = asdict(config)
    inputs_fp = {
        "structure": structure_provider.fingerprint(rcsb_id),
        "landmarks": landmark_provider.fingerprint(rcsb_id),
    }

    struct_runs_dir = SETTINGS.runs_root / rcsb_id
    struct_runs_dir.mkdir(parents=True, exist_ok=True)

    run_id = compute_run_id(
        rcsb_id=rcsb_id,
        pipeline_version=_pipeline_version(),
        inputs_fp=inputs_fp,
        config_resolved=config_resolved,
        runs_dir=struct_runs_dir,
    )

    run_dir = struct_runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    manifest = RunManifest(
        rcsb_id=rcsb_id,
        run_id=run_id,
        pipeline_version=_pipeline_version(),
        inputs={"fingerprints": inputs_fp},
        config_resolved=config_resolved,
    )
    store = LocalRunStore(run_dir=run_dir, manifest=manifest)

    ctx = StageContext(
        run_id=run_id,
        rcsb_id=rcsb_id,
        config=config,
        store=store,
        inputs={
            "structure_provider": structure_provider,
            "landmark_provider": landmark_provider,
            "stage_cache": LocalStageCache(SETTINGS.cache_root),
            "inputs_fp": inputs_fp,
        },
    )

    pipeline = Pipeline([
        Stage00Inputs(),
        Stage10Landmarks(),
        Stage20ExteriorShell(),
        Stage30RegionAtoms(),
        Stage40EmptySpace(),
        Stage50Clustering(),
        Stage55GridRefine(),
        Stage70MeshValidate(),
    ])

    pipeline.run(ctx)
    return ctx

```

libnpet/__main__.py
```py
# libnpet/__main__.py
"""
npet2 - Ribosome nascent peptide exit tunnel geometry pipeline

Two modes:

  API mode   -- pass an RCSB ID, everything is fetched automatically.
                mmCIF is downloaded from RCSB, profile/landmarks from ribosome.xyz.

    npet2 run 5NWY --api-url http://localhost:8000

  File mode  -- point at local files explicitly.
                Any omitted file falls back to the API.

    npet2 run 5NWY \\
        --mmcif   /data/5NWY.cif \\
        --profile /data/5NWY.json \\
        --ptc     /data/5NWY_PTC.json \\
        --constriction /data/5NWY_CONSTRICTION_SITE.json

  Mixed      -- local mmCIF, API for the rest (useful when you have the
                structure locally but not the landmarks yet)

    npet2 run 5NWY --mmcif /data/5NWY.cif --api-url http://localhost:8000

Advanced config overrides are available via --config-json or individual flags.
Run `npet2 show-config` to see all defaults.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Optional


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="npet2",
        description="Ribosome exit tunnel geometry pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = p.add_subparsers(dest="command")

    # -------------------------------------------------------------------------
    # setup
    # -------------------------------------------------------------------------
    setup_p = sub.add_parser("setup", help="Build required external binaries (PoissonRecon)")
    setup_p.add_argument("--data-dir", metavar="DIR",
                         help="Root directory for npet2 data (default: ~/.npet2 or $NPET2_ROOT)")

    # -------------------------------------------------------------------------
    # run
    # -------------------------------------------------------------------------
    run_p = sub.add_parser(
        "run",
        help="Run the pipeline on one or more structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Target structures
    id_grp = run_p.add_argument_group("target structures (provide at least one)")
    id_grp.add_argument("rcsb_ids", nargs="*", metavar="RCSB_ID",
                        help="One or more RCSB IDs (e.g. 5NWY 7K00)")
    id_grp.add_argument("--rcsb_id", dest="rcsb_id_flag", metavar="RCSB_ID",
                        help="Single RCSB ID (alternative to positional)")
    id_grp.add_argument("--from-file", metavar="FILE",
                        help="Text file with one RCSB ID per line")

    # Data sources
    src_grp = run_p.add_argument_group(
        "data sources",
        "All inputs default to the ribosome.xyz API when not provided as files.\n"
        "The mmCIF is downloaded from RCSB if --mmcif is omitted.",
    )
    src_grp.add_argument("--api-url", metavar="URL", default=None,
                         help="ribosome.xyz API base URL (default: $NPET2_RIBOXYZ_API_URL "
                              "or http://localhost:8000)")
    src_grp.add_argument("--mmcif", metavar="PATH",
                         help="mmCIF file path. If omitted, downloaded from RCSB and cached.")
    src_grp.add_argument("--profile", metavar="PATH",
                         help="RibosomeStructure JSON. If omitted, fetched from API.")
    src_grp.add_argument("--ptc", metavar="PATH",
                         help="PTC landmark JSON {\"location\":[x,y,z],...}. "
                              "If omitted, fetched from API.")
    src_grp.add_argument("--constriction", metavar="PATH",
                         help="Constriction site JSON {\"location\":[x,y,z]}. "
                              "If omitted, fetched from API.")
    src_grp.add_argument("--data-dir", metavar="DIR",
                         help="Root directory for all npet2 data: mmcif cache, runs, build "
                              "(default: ~/.npet2 or $NPET2_ROOT)")

    # Output
    run_p.add_argument("--output-dir", metavar="DIR",
                       help="Root directory for run outputs (default: $NPET2_RUNS_ROOT)")

    # Parallelism
    run_p.add_argument("--workers", "-j", type=int, default=1, metavar="N",
                       help="Parallel workers for batch runs (default: 1)")

    # Advanced config overrides
    cfg = run_p.add_argument_group(
        "advanced config overrides",
        "Override individual pipeline parameters. "
        "Use --config-json to supply a full RunConfig JSON instead."
    )
    cfg.add_argument("--config-json", metavar="PATH",
                     help="Path to a RunConfig JSON (see `npet2 show-config`)")
    cfg.add_argument("--cylinder-radius", type=float, metavar="A")
    cfg.add_argument("--cylinder-height", type=float, metavar="A")
    cfg.add_argument("--ptc-extension", type=float, metavar="A")
    cfg.add_argument("--voxel-size", type=float, metavar="A",
                     help="Level-0 grid voxel size in Angstroms")
    cfg.add_argument("--refine-voxel-size", type=float, metavar="A")
    cfg.add_argument("--no-mesh", action="store_true", help="Disable mesh generation")
    cfg.add_argument("--no-refine", action="store_true", help="Skip Stage55 grid refinement")
    cfg.add_argument("--dbscan-coarse-eps", type=float, metavar="A")
    cfg.add_argument("--dbscan-coarse-min-samples", type=int)
    cfg.add_argument("--dbscan-refine-eps", type=float, metavar="A")
    cfg.add_argument("--dbscan-refine-min-samples", type=int)

    # -------------------------------------------------------------------------
    # show-config
    # -------------------------------------------------------------------------
    sub.add_parser("show-config", help="Print the default RunConfig as JSON")

    return p


def _build_config(args) -> "RunConfig":
    from libnpet.core.config import RunConfig, GridLevelConfig

    if args.config_json:
        data = json.loads(Path(args.config_json).read_text())
        if "grid_levels" in data:
            data["grid_levels"] = [GridLevelConfig(**gl) for gl in data["grid_levels"]]
        return RunConfig(**data)

    kwargs = {}
    if args.cylinder_radius is not None:
        kwargs["cylinder_radius_A"] = args.cylinder_radius
    if args.cylinder_height is not None:
        kwargs["cylinder_height_A"] = args.cylinder_height
    if args.ptc_extension is not None:
        kwargs["cylinder_ptc_extension_A"] = args.ptc_extension
    if args.refine_voxel_size is not None:
        kwargs["refine_voxel_size_A"] = args.refine_voxel_size
    if args.no_mesh:
        kwargs["mesh_level0_enable"] = False
        kwargs["mesh_level1_enable"] = False
    if args.dbscan_coarse_eps is not None:
        kwargs["dbscan_level0_coarse_eps_A"] = args.dbscan_coarse_eps
    if args.dbscan_coarse_min_samples is not None:
        kwargs["dbscan_level0_coarse_min_samples"] = args.dbscan_coarse_min_samples
    if args.dbscan_refine_eps is not None:
        kwargs["dbscan_level0_refine_eps_A"] = args.dbscan_refine_eps
    if args.dbscan_refine_min_samples is not None:
        kwargs["dbscan_level0_refine_min_samples"] = args.dbscan_refine_min_samples
    if args.voxel_size is not None:
        kwargs["grid_levels"] = [
            GridLevelConfig(name="level_0", voxel_size_A=args.voxel_size,
                            occupancy_backend="legacy_kdtree"),
        ]
    return RunConfig(**kwargs)


def _collect_rcsb_ids(args) -> list[str]:
    ids = list(args.rcsb_ids) if args.rcsb_ids else []
    if getattr(args, "rcsb_id_flag", None):
        ids.append(args.rcsb_id_flag)
    if args.from_file:
        p = Path(args.from_file)
        if not p.exists():
            print(f"Error: --from-file {p} does not exist", file=sys.stderr)
            sys.exit(1)
        for line in p.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line)
    if not ids:
        print("Error: no RCSB IDs specified.", file=sys.stderr)
        print("  Provide as positional args, --rcsb_id, or --from-file.", file=sys.stderr)
        sys.exit(1)
    return [x.upper() for x in ids]


def _make_providers(args, rcsb_id: str):
    from libnpet.adapters.standalone_providers import (
        FileStructureProvider,
        FileLandmarkProvider,
        _download_mmcif,
    )
    from libnpet.core.config import SETTINGS

    api_base = args.api_url or SETTINGS.riboxyz_api_base

    # mmCIF: explicit path, or auto-download from RCSB
    if args.mmcif:
        mmcif_path = Path(args.mmcif)
        if not mmcif_path.exists():
            print(f"Error: --mmcif {mmcif_path} does not exist", file=sys.stderr)
            sys.exit(1)
    else:
        mmcif_path = _download_mmcif(rcsb_id)

    profile_path = Path(args.profile) if args.profile else None
    ptc_path     = Path(args.ptc)     if args.ptc     else None
    constr_path  = Path(args.constriction) if args.constriction else None

    sp = FileStructureProvider(
        mmcif_path=mmcif_path,
        profile_path=profile_path,
        api_base=api_base,
    )
    lp = FileLandmarkProvider(
        ptc_path=ptc_path,
        constriction_path=constr_path,
        api_base=api_base,
    )
    return sp, lp


def _apply_data_dir(args) -> None:
    """Override SETTINGS.npet2_root and derived paths if --data-dir was given."""
    data_dir = getattr(args, "data_dir", None)
    if not data_dir:
        return
    import libnpet.core.config as cfg_mod
    root = Path(data_dir)
    cfg_mod.SETTINGS = cfg_mod.Settings(
        npet2_root        = root,
        runs_root         = root / "runs",
        cache_root        = root / "cache",
        poisson_recon_bin = str(root / "bin" / "PoissonRecon"),
        riboxyz_api_base  = cfg_mod.SETTINGS.riboxyz_api_base,
    )
def _run_single(rcsb_id: str, args, config, output_root: Optional[Path]) -> dict:
    from libnpet.run import run_npet2

    _apply_data_dir(args) 
    try:
        sp, lp = _make_providers(args, rcsb_id)

        if output_root:
            import libnpet.core.config as cfg_mod
            cfg_mod.SETTINGS = cfg_mod.Settings(
                runs_root=output_root,
                npet2_root=cfg_mod.SETTINGS.npet2_root,
                cache_root=cfg_mod.SETTINGS.cache_root,
                poisson_recon_bin=cfg_mod.SETTINGS.poisson_recon_bin,
                riboxyz_api_base=cfg_mod.SETTINGS.riboxyz_api_base,
            )

        ctx = run_npet2(rcsb_id, config, structure_provider=sp, landmark_provider=lp)
        return {"rcsb_id": rcsb_id, "status": "success", "run_dir": str(ctx.store.run_dir)}
    except Exception as e:
        return {
            "rcsb_id": rcsb_id,
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


def _run_worker(packed_args: tuple) -> dict:
    rcsb_id, args_ns, config_dict, output_root_str = packed_args
    from libnpet.core.config import RunConfig, GridLevelConfig

    if "grid_levels" in config_dict:
        config_dict["grid_levels"] = [GridLevelConfig(**gl) for gl in config_dict["grid_levels"]]
    config = RunConfig(**config_dict)
    output_root = Path(output_root_str) if output_root_str else None
    return _run_single(rcsb_id, args_ns, config, output_root)

def _cmd_setup(args=None):
    import platform
    import shutil
    import subprocess

    if args:
        _apply_data_dir(args)

    from libnpet.core.config import SETTINGS

    dest = Path(SETTINGS.poisson_recon_bin)

    print("npet2 setup: checking for PoissonRecon...")
    if dest.exists():
        print(f"  already built: {dest}")
        return

    system = platform.system()

    # prerequisites
    missing = []
    for tool in ("git", "make"):
        if not shutil.which(tool):
            missing.append(tool)
    if not (shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")):
        missing.append("c++ compiler (g++ or clang++)")

    if missing:
        print("  missing required tools:", file=sys.stderr)
        for m in missing:
            print(f"    - {m}", file=sys.stderr)
        print(
            "\n  Install them and retry.\n"
            "\n  macOS:\n"
            "    xcode-select --install\n"
            "\n  Ubuntu/Debian:\n"
            "    apt-get install -y git build-essential",
            file=sys.stderr,
        )
        sys.exit(1)

    # libjpeg: required by PoissonRecon's image I/O
    if system == "Darwin":
        brew = shutil.which("brew")
        # check if jpeg headers are findable
        import subprocess as sp
        jpeg_ok = sp.run(
            ["clang++", "-x", "c++", "-include", "jpeglib.h", "-E", "-"],
            input="", capture_output=True, text=True
        ).returncode == 0
        if not jpeg_ok:
            if brew:
                print("  libjpeg not found, installing via brew ...")
                subprocess.run(["brew", "install", "jpeg"], check=True)
            else:
                print(
                    "  libjpeg not found and brew is not available.\n"
                    "  Install brew (https://brew.sh) then run: brew install jpeg",
                    file=sys.stderr,
                )
                sys.exit(1)
    else:
        # on Linux check for the header directly
        jpeg_header = Path("/usr/include/jpeglib.h")
        if not jpeg_header.exists():
            print(
                "  libjpeg-dev not found.\n"
                "  Install it with: apt-get install -y libjpeg-dev",
                file=sys.stderr,
            )
            sys.exit(1)

    build_root = Path(SETTINGS.npet2_root) / "build" / "PoissonRecon"
    build_root.mkdir(parents=True, exist_ok=True)
    dest.parent.mkdir(parents=True, exist_ok=True)

    repo_dir = build_root / "src"
    if not repo_dir.exists():
        print("  cloning mkazhdan/PoissonRecon ...")
        subprocess.run(
            ["git", "clone", "--depth", "1",
             "https://github.com/mkazhdan/PoissonRecon.git", str(repo_dir)],
            check=True,
        )
    else:
        print("  source already cloned, skipping git clone")

    if system == "Darwin":
        compiler = "clang"
        print("  macOS detected: building with COMPILER=clang (no OpenMP)")
        # brew jpeg lives under a non-standard prefix on Apple Silicon / Intel
        brew_prefix = subprocess.run(
            ["brew", "--prefix", "jpeg"],
            capture_output=True, text=True
        ).stdout.strip()
        extra_env = {
            **__import__("os").environ,
            "CPATH": f"{brew_prefix}/include",
            "LIBRARY_PATH": f"{brew_prefix}/lib",
        }
    else:
        compiler = "gcc"
        print("  Linux detected: building with COMPILER=gcc (OpenMP enabled)")
        extra_env = None

    print("  building PoissonRecon (this takes a minute) ...")
    subprocess.run(
        ["make", "-j4", "poissonrecon", f"COMPILER={compiler}"],
        cwd=str(repo_dir),
        env=extra_env,
        check=True,
    )

    candidates = [c for c in (repo_dir / "Bin").rglob("PoissonRecon") if c.is_file()]
    if not candidates:
        print("  build finished but PoissonRecon binary not found", file=sys.stderr)
        print(f"  look in {repo_dir / 'Bin'} and set NPET2_POISSON_RECON_BIN manually",
              file=sys.stderr)
        sys.exit(1)

    shutil.copy2(str(candidates[0]), str(dest))
    dest.chmod(0o755)
    print(f"  PoissonRecon installed to {dest}")
    print(f"  all npet2 data will live under {SETTINGS.npet2_root}")
    print("  setup complete. You can now run: npet2 run <RCSB_ID>")


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "setup":
        _cmd_setup(args)
        return

    if args.command == "show-config":
        from libnpet.core.config import RunConfig
        print(json.dumps(asdict(RunConfig()), indent=2))
        return

    if args.command == "run":
        _apply_data_dir(args)
        rcsb_ids    = _collect_rcsb_ids(args)
        config      = _build_config(args)
        output_root = Path(args.output_dir) if args.output_dir else None
        n_workers   = min(args.workers, len(rcsb_ids))

        print(f"npet2: processing {len(rcsb_ids)} structure(s), workers={n_workers}")

        if n_workers <= 1:
            results = []
            for rid in rcsb_ids:
                r = _run_single(rid, args, config, output_root)
                results.append(r)
                status = r["status"]
                print(f"  {rid}: {status}" + (
                    f" -> {r.get('run_dir', '')}" if status == "success"
                    else f"\n{r.get('traceback', r.get('error', ''))}"
                ))
        else:
            config_dict     = asdict(config)
            output_root_str = str(output_root) if output_root else None
            packed = [(rid, args, config_dict, output_root_str) for rid in rcsb_ids]
            results = []
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_worker, p): p[0] for p in packed}
                for fut in as_completed(futures):
                    rid = futures[fut]
                    try:
                        r = fut.result()
                    except Exception as e:
                        r = {"rcsb_id": rid, "status": "failed", "error": str(e)}
                    results.append(r)
                    status = r["status"]
                    print(f"  {rid}: {status}" + (
                        f" -> {r.get('run_dir', '')}" if status == "success"
                        else f" ({r.get('error', '')})"
                    ))

        ok   = sum(1 for r in results if r["status"] == "success")
        fail = len(results) - ok
        print(f"\nnpet2: {ok} succeeded, {fail} failed out of {len(results)}")
        if fail > 0:
            sys.exit(1)



if __name__ == "__main__":
    main()
```

libnpet/__init__.py
```py

```

libnpet/stages/legacy_minimal.py
```py
# ribctl/lib/npet2/stages/legacy_minimal.py
# npet2/stages/legacy_minimal.py

from __future__ import annotations
import json
from pathlib import Path
import time
from typing import Any, Dict, List, Tuple
import numpy as np
import pyvista as pv
import open3d as o3d

from libnpet.backends.grid_occupancy import (
    connected_components_3d,
    occupancy_via_edt,
)
from libnpet.backends.meshing import save_mesh_with_ascii
from libnpet.core.cache import StageCacheKey
from libnpet.core.pipeline import Stage
from libnpet.core.ribosome_types import RibosomeProfile
from libnpet.core.structure_selection import (
    intersect_with_first_assembly,
    ribosome_wall_auth_asym_ids,
    tunnel_debris_chains,
    atom_inclusion_policy,
)
from libnpet.core.types import StageContext, ArtifactType

from scipy import ndimage

from libnpet.backends.geometry import (
    cif_to_point_cloud,
    fast_normal_estimation,
    quick_surface_points,
    validate_mesh_pyvista,
    apply_poisson_reconstruction,
    filter_residues_parallel,
    transform_points_to_C0,
    transform_points_from_C0,
    create_point_cloud_mask,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    estimate_normals,
)
from libnpet.stages.grid_refine import (
    _make_bbox_grid,
    _points_to_ijk,
    _valid_ijk,
    _voxel_centers_from_indices,
)

def _residues_from_chain_ids(structure, chain_ids: set[str]):
    model = structure[0]
    residues = []
    for cid in chain_ids:
        if cid not in model:
            continue
        chain = model[cid]
        for r in chain.get_residues():
            if len(getattr(r, "child_list", [])) == 0:
                continue
            residues.append(r)
    return residues


def _pick_tunnel_cluster(
    clusters: dict[int, list],
    constr: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Pick the cluster whose points are closest to the constriction site."""
    constr = np.asarray(constr, dtype=np.float32).reshape(1, 3)
    best_id = -1
    best_dist = float("inf")
    for cid, pts_list in clusters.items():
        if cid == -1:
            continue
        pts = np.asarray(pts_list, dtype=np.float32)
        if pts.shape[0] == 0:
            continue
        dists = np.linalg.norm(pts - constr, axis=1)
        min_dist = float(dists.min())
        if min_dist < best_dist:
            best_dist = min_dist
            best_id = cid
    if best_id == -1:
        raise ValueError("No valid clusters found")
    print(f"  [cluster_select] picked cluster {best_id} "
          f"(n={len(clusters[best_id]):,}, dist_to_constriction={best_dist:.1f}A)")
    return np.asarray(clusters[best_id], dtype=np.float32), best_id

def _get_biopython_structure(ctx: StageContext):
    """Get biopython structure from ctx, or parse mmcif on demand."""
    bs = ctx.inputs.get("biopython_structure")
    if bs is not None:
        return bs
    from Bio.PDB.MMCIFParser import FastMMCIFParser
    mmcif_path = ctx.require("mmcif_path")
    bs = FastMMCIFParser(QUIET=True).get_structure(ctx.rcsb_id, mmcif_path)
    ctx.inputs["biopython_structure"] = bs
    return bs

class Stage20ExteriorShell(Stage):
    key = "20_exterior_shell"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "d3d_alpha": c.alpha_d3d_alpha,
            "d3d_tol": c.alpha_d3d_tol,
            "d3d_offset": c.alpha_d3d_offset,
            "kdtree_radius": c.alpha_kdtree_radius,
            "max_nn": c.alpha_max_nn,
            "tangent_k": c.alpha_tangent_planes_k,
            "poisson_depth": c.alpha_poisson_depth,
            "poisson_ptweight": c.alpha_poisson_ptweight,
            "fill_holes": c.alpha_fill_holes,
        }

    def run(self, ctx: StageContext) -> None:
        stage_cache = ctx.require("stage_cache")
        inputs_fp = ctx.require("inputs_fp")
        params = self.params(ctx)

        key = StageCacheKey(
            stage=self.key,
            inputs_fp={"structure": inputs_fp["structure"]},
            params=params,
            impl_version="v1",
        )

        stage_dir = ctx.store.stage_dir(self.key)
        cached_files = [
            "alpha_shell.ply",
            "alpha_shell_quality.json",
            "alpha_normals.ply",
            "alpha_surface_points.npy",
            "ribosome_ptcloud.npy",
        ]

        if stage_cache.has(
            key, required=["alpha_shell.ply", "alpha_shell_quality.json"]
        ):
            stage_cache.copy_into(key, stage_dir, cached_files)

            quality = json.loads((stage_dir / "alpha_shell_quality.json").read_text())
            ctx.inputs["alpha_shell_path"] = str(stage_dir / "alpha_shell.ply")
            ctx.inputs["alpha_shell_watertight"] = bool(
                quality.get("watertight", False)
            )

            ctx.store.register_file(
                name="alpha_shell_mesh",
                stage=self.key,
                type=ArtifactType.PLY_MESH,
                path=stage_dir / "alpha_shell.ply",
            )
            ctx.store.register_file(
                name="alpha_shell_quality",
                stage=self.key,
                type=ArtifactType.JSON,
                path=stage_dir / "alpha_shell_quality.json",
            )
            return

        c = ctx.config
        profile: RibosomeProfile = ctx.require("profile")
        cifpath = Path(ctx.require("mmcif_path"))

        ptcloud_path = stage_dir / "ribosome_ptcloud.npy"
        surface_pts_path = stage_dir / "alpha_surface_points.npy"
        normals_pcd_path = stage_dir / "alpha_normals.ply"
        mesh_path = stage_dir / "alpha_shell.ply"
        quality_path = stage_dir / "alpha_shell_quality.json"

        wall = ribosome_wall_auth_asym_ids(
            profile,
            exclude_trna=bool(getattr(ctx.config, "occupancy_exclude_trna", True)),
            extra_exclude=tunnel_debris_chains(ctx.rcsb_id, profile),
        )
        wall = intersect_with_first_assembly(profile, wall)

        ptcloud = cif_to_point_cloud(str(cifpath), sorted(wall), do_atoms=True)

        np.save(ptcloud_path, ptcloud)
        ctx.store.register_file(
            name="ribosome_ptcloud",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=ptcloud_path,
        )

        surface_pts = quick_surface_points(
            ptcloud, c.alpha_d3d_alpha, c.alpha_d3d_tol, c.alpha_d3d_offset
        ).astype(np.float32)
        np.save(surface_pts_path, surface_pts)
        ctx.store.register_file(
            name="alpha_surface_points",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=surface_pts_path,
        )

        normal_estimated_pcd = fast_normal_estimation(
            surface_pts, c.alpha_kdtree_radius, c.alpha_max_nn, c.alpha_tangent_planes_k
        )

        center = normal_estimated_pcd.get_center()
        normal_estimated_pcd.orient_normals_towards_camera_location(
            camera_location=center
        )
        normal_estimated_pcd.normals = o3d.utility.Vector3dVector(
            -np.asarray(normal_estimated_pcd.normals)
        )

        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        ctx.store.register_file(
            name="alpha_normals_pcd",
            stage=self.key,
            type=ArtifactType.PLY_PCD,
            path=normals_pcd_path,
        )

        apply_poisson_reconstruction(
            str(normals_pcd_path),
            mesh_path,
            recon_depth=c.alpha_poisson_depth,
            recon_pt_weight=c.alpha_poisson_ptweight,
        )

        mesh = pv.read(mesh_path)
        mesh = mesh.fill_holes(c.alpha_fill_holes)
        # mesh = mesh.connectivity(largest=True).triangulate()
        mesh = mesh.connectivity(extraction_mode='largest').triangulate()
        mesh.save(mesh_path)

        watertight = validate_mesh_pyvista(mesh)

        quality = {
            "watertight": bool(watertight),
            "n_points": int(mesh.n_points),
            "n_faces": int(mesh.n_faces_strict),
            "open_edges": int(mesh.n_open_edges),
            "is_manifold": bool(mesh.is_manifold),
            "bounds": list(mesh.bounds),
        }
        quality_path.write_text(json.dumps(quality, indent=2))
        ctx.store.register_file(
            name="alpha_shell_quality",
            stage=self.key,
            type=ArtifactType.JSON,
            path=quality_path,
        )

        ctx.store.register_file(
            name="alpha_shell_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
        )
        ctx.inputs["alpha_shell_path"] = str(mesh_path)
        ctx.inputs["alpha_shell_watertight"] = bool(watertight)
        if watertight:
            stage_cache.put_from(key, stage_dir, cached_files)

class Stage30RegionAtoms(Stage):
    key = "30_region_atoms"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {"radius_A": c.cylinder_radius_A, "height_A": c.cylinder_height_A}

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        profile: RibosomeProfile = ctx.require("profile")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        policy = atom_inclusion_policy(profile, c, ctx.rcsb_id)

        occ_chain_ids = policy["wall_chain_ids"]
        seed_chain_ids = set(occ_chain_ids)

        structure = _get_biopython_structure(ctx)

        residues_seed = _residues_from_chain_ids(structure, seed_chain_ids)
        residues_occ = _residues_from_chain_ids(structure, occ_chain_ids)

        residues_seed = filter_residues_parallel(
            residues=residues_seed,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
        )
        residues_occ = filter_residues_parallel(
            residues=residues_occ,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
        )

        seed_points = np.asarray(
            [atom.get_coord() for r in residues_seed for atom in r.child_list],
            dtype=np.float32,
        )
        occ_points = np.asarray(
            [atom.get_coord() for r in residues_occ for atom in r.child_list],
            dtype=np.float32,
        )

        stage_dir = ctx.store.stage_dir(self.key)

        out_seed = stage_dir / "region_atom_xyz.npy"
        np.save(out_seed, seed_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_seed,
            meta={
                "n": int(seed_points.shape[0]),
                "note": "seed atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz"] = seed_points

        out_occ = stage_dir / "region_atom_xyz_occ.npy"
        np.save(out_occ, occ_points)
        ctx.store.register_file(
            name="region_atom_xyz_occ",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_occ,
            meta={
                "n": int(occ_points.shape[0]),
                "note": "occupancy atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz_occ"] = occ_points

        policy_record = {
            "occupancy_chain_mode": getattr(c, "occupancy_chain_mode", "walls_only"),
            "exclude_trna": bool(getattr(c, "occupancy_exclude_trna", True)),
            "wall_chain_ids": sorted(occ_chain_ids),
            "excluded_chains": {k: v for k, v in policy["reasons"].items()},
            "policy_summary": (
                "INCLUDED: ribosomal proteins + rRNAs (with modified residues). "
                "EXCLUDED: waters, ions, nonpolymer ligands, tRNAs, debris chains."
            ),
        }
        (stage_dir / "atom_selection_policy.json").write_text(
            json.dumps(policy_record, indent=2)
        )

        print(
            f"[{self.key}] seed_atoms={seed_points.shape[0]:,} "
            f"occ_atoms={occ_points.shape[0]:,} occ_chains={len(occ_chain_ids)}"
        )
        if policy["reasons"]:
            excluded_summary = ", ".join(
                f"{k}({v})" for k, v in sorted(policy["reasons"].items())
            )
            print(f"[{self.key}] excluded: {excluded_summary}")

class Stage40EmptySpace(Stage):
    key = "40_empty_space"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "grid_levels": [
                {
                    "name": gl.name,
                    "voxel_size_A": gl.voxel_size_A,
                    "backend": gl.occupancy_backend,
                    "atom_radius_mode": getattr(gl, "atom_radius_mode", "uniform"),
                    "uniform_atom_radius_A": getattr(gl, "uniform_atom_radius_A", None),
                }
                for gl in c.grid_levels
            ],
            "radius_A": c.cylinder_radius_A,
            "height_A": c.cylinder_height_A,
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import numpy as np
        import pyvista as pv


        from libnpet.backends.grid_occupancy import (
            make_cylinder_grid,
            cylinder_mask,
            occupancy_via_edt,
            empty_points_from_mask,
            save_grid_npy,
            get_occupied_voxel_centers,
        )

        c = ctx.config

        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        # Use filtered atoms for clustering seed reference
        region_xyz_filtered = np.asarray(
            ctx.require("region_atom_xyz"), dtype=np.float32
        )

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr)

        shell = pv.read(alpha_shell_path)
        if not isinstance(shell, pv.PolyData):
            shell = shell.extract_surface()
        if not shell.is_all_triangles:
            shell = shell.triangulate()

        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        stage_dir = ctx.store.stage_dir(self.key)

        clip_note = {
            "alpha_shell_path": str(alpha_shell_path),
            "alpha_shell_watertight": watertight,
            "clipping_mode": "select_enclosed_points(check_surface=True)"
            if watertight
            else "select_enclosed_points(check_surface=False) [fallback; shell not watertight]",
            "occupancy_atoms": "region_atom_xyz_all (includes ALL chains)",
        }
        p_clip_note = stage_dir / "clipping_note.json"
        p_clip_note.write_text(json.dumps(clip_note, indent=2))
        ctx.store.register_file(
            name="clipping_note",
            stage=self.key,
            type=ArtifactType.JSON,
            path=p_clip_note,
        )

        last_empty = None

        for gl in c.grid_levels:
            backend = gl.occupancy_backend

            if backend == "legacy_kdtree":
                mask, (x, y, z) = create_point_cloud_mask(
                    region_c0,
                    radius=c.cylinder_radius_A,
                    height=c.cylinder_height_A,
                    voxel_size=gl.voxel_size_A,
                    radius_around_point=gl.uniform_atom_radius_A,
                    z_min=z_min,
                )

                idx = np.where(~mask)
                empty_c0 = np.column_stack((x[idx[0]], y[idx[1]], z[idx[2]])).astype(
                    np.float32
                )

            elif backend == "edt":
                grid = make_cylinder_grid(
                    radius_A=float(c.cylinder_radius_A),
                    height_A=float(c.cylinder_height_A),
                    voxel_A=float(gl.voxel_size_A),
                    z_min=z_min,
                )

                occ = occupancy_via_edt(
                    region_c0,
                    grid,
                    atom_radius_A=float(gl.uniform_atom_radius_A),
                )

                cyl2d = cylinder_mask(grid, radius_A=float(c.cylinder_radius_A))
                cyl = np.broadcast_to(cyl2d, grid.shape)

                occ = occ | (~cyl)

                empty_mask = ~occ
                empty_c0 = empty_points_from_mask(grid, empty_mask & cyl)

                save_grid_npy(
                    grid, occ, stage_dir / f"occupancy_grid_{gl.name}", compress=False
                )
                save_grid_npy(
                    grid,
                    (empty_mask & cyl),
                    stage_dir / f"empty_mask_{gl.name}",
                    compress=False,
                )

                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"occupancy_grid_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"occupancy_grid_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"empty_mask_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"empty_mask_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )

                try:
                    occ_centers_c0 = get_occupied_voxel_centers(grid, occ).astype(
                        np.float32
                    )
                    occ_centers_world = transform_points_from_C0(
                        occ_centers_c0, ptc, constr
                    ).astype(np.float32)
                    p_occ = stage_dir / f"occupied_voxels_{gl.name}.npy"
                    np.save(p_occ, occ_centers_world)
                    ctx.store.register_file(
                        name=f"occupied_voxels_{gl.name}",
                        stage=self.key,
                        type=ArtifactType.NUMPY,
                        path=p_occ,
                        meta={
                            "voxel_size_A": gl.voxel_size_A,
                            "n": int(occ_centers_world.shape[0]),
                        },
                    )
                except Exception:
                    pass

            else:
                raise ValueError(
                    f"Grid level {gl.name}: unsupported backend {backend} "
                    f"(supported: legacy_kdtree, edt)"
                )

            empty_world = transform_points_from_C0(empty_c0, ptc, constr).astype(
                np.float32
            )

            p_pre = stage_dir / f"empty_points_{gl.name}_preclip.npy"
            np.save(p_pre, empty_world)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}_preclip",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=p_pre,
                meta={"voxel_size_A": gl.voxel_size_A, "n": int(empty_world.shape[0])},
            )

            if empty_world.shape[0] == 0:
                inside = empty_world
            else:
                pts_poly = pv.PolyData(empty_world)
                sel = pts_poly.select_enclosed_points(shell, check_surface=watertight)
                inside = empty_world[sel["SelectedPoints"] == 1].astype(np.float32)

            out = stage_dir / f"empty_points_{gl.name}.npy"
            np.save(out, inside)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "voxel_size_A": gl.voxel_size_A,
                    "n": int(inside.shape[0]),
                    "backend": backend,
                    "alpha_shell_watertight": watertight,
                },
            )

            ctx.inputs[f"empty_points_{gl.name}"] = inside
            last_empty = inside

        ctx.inputs["empty_points"] = last_empty

class Stage50Clustering(Stage):
    """
    DBSCAN clustering on level_0 (coarse grid, typically 1.0A).
    
    Two-pass strategy:
      1. Coarse DBSCAN: merge regions, bridge gaps
      2. Refine DBSCAN: tighten on largest cluster from pass 1
    
    Cluster selection uses axial proximity to the PTC-Constriction axis
    rather than raw point count, which prevents the inter-subunit space
    from being picked over the actual tunnel.
    
    Optionally generates a mesh from the refined cluster.
    
    Outputs:
      stage/50_clustering/
        coarse/{points.npy, labels.npy, cluster_*.npy, index.json}
        refine/{points.npy, labels.npy, cluster_*.npy, index.json}
        largest_cluster.npy
        refined_cluster.npy
        mesh_level_0.ply (if enabled)
    """
    key = "50_clustering"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "coarse_eps_A": c.dbscan_level0_coarse_eps_A,
            "coarse_min_samples": c.dbscan_level0_coarse_min_samples,
            "refine_eps_A": c.dbscan_level0_refine_eps_A,
            "refine_min_samples": c.dbscan_level0_refine_min_samples,
            "mesh_enable": bool(getattr(c, "mesh_level0_enable", True)),
            "cluster_selection": "axial_proximity",
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import time

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        # Input: empty points from Stage40
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)
        if empty_pts.ndim != 2 or empty_pts.shape[1] != 3 or empty_pts.shape[0] == 0:
            raise ValueError(f"[{self.key}] empty_points must be (N,3) and non-empty, got {empty_pts.shape}")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        print(f"[{self.key}] empty_points n={empty_pts.shape[0]:,}")

        # -----------------------
        # PASS 1: Coarse DBSCAN
        # -----------------------
        t0 = time.perf_counter()
        db_coarse, clusters_coarse = DBSCAN_capture(
            empty_pts, 
            c.dbscan_level0_coarse_eps_A, 
            c.dbscan_level0_coarse_min_samples
        )
        labels_coarse = np.asarray(db_coarse.labels_, dtype=np.int32)
        dt0 = time.perf_counter() - t0
        
        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        self._save_dbscan_pass(
            stage_dir / "coarse",
            empty_pts,
            labels_coarse,
            clusters_coarse,
            eps=c.dbscan_level0_coarse_eps_A,
            min_samples=c.dbscan_level0_coarse_min_samples,
        )

        # Pick tunnel cluster from coarse (axial proximity, not largest)


        largest, largest_id = _pick_tunnel_cluster(clusters_coarse, constr)
        if largest.shape[0] == 0:
            raise ValueError(f"[{self.key}] tunnel cluster is empty")

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest)
        ctx.store.register_file(
            name="largest_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_largest,
            meta={"cluster_id": int(largest_id), "n": int(largest.shape[0]),
                  "selection": "axial_proximity"},
        )

        # -----------------------
        # PASS 2: Refine DBSCAN
        # -----------------------
        t1 = time.perf_counter()
        db_refine, clusters_refine = DBSCAN_capture(
            largest,
            c.dbscan_level0_refine_eps_A,
            c.dbscan_level0_refine_min_samples,
        )
        labels_refine = np.asarray(db_refine.labels_, dtype=np.int32)
        dt1 = time.perf_counter() - t1
        
        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        self._save_dbscan_pass(
            stage_dir / "refine",
            largest,
            labels_refine,
            clusters_refine,
            eps=c.dbscan_level0_refine_eps_A,
            min_samples=c.dbscan_level0_refine_min_samples,
        )

        # Pick tunnel cluster from refine


        refined, refined_id = _pick_tunnel_cluster(clusters_refine, constr)
        if refined.shape[0] == 0:
            raise ValueError(f"[{self.key}] refined cluster is empty")

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0]),
                  "selection": "axial_proximity"},
        )

        # Output for downstream
        ctx.inputs["refined_cluster"] = refined
        ctx.inputs["largest_cluster"] = largest

        print(f"[{self.key}] winner: coarse={largest.shape[0]:,} -> refine={refined.shape[0]:,}")

        # -----------------------
        # Optional: Mesh level_0
        # -----------------------
        if c.mesh_level0_enable:
            self._generate_mesh(ctx, refined, level_name="level_0")

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        clusters_dict: dict[int, list],
        eps: float,
        min_samples: int,
    ) -> None:
        """Save DBSCAN pass artifacts."""
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        for lab, plist in clusters_dict.items():
            lab = int(lab)
            if lab == -1:
                continue
            arr = np.asarray(plist, dtype=np.float32)
            if arr.size > 0:
                np.save(pass_dir / f"cluster_id{lab}.npy", arr)

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        import time
        from libnpet.backends.meshing import (
            mesh_from_binary_volume,
            voxelize_points,
            clip_mesh_to_atom_clearance,
            save_mesh_with_ascii,
        )

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        print(f"[{self.key}] generating mesh for {level_name}...")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        pts_c0 = transform_points_to_C0(points, ptc, constr).astype(np.float32)
        voxel = 1.0

        t0 = time.perf_counter()
        mask, origin = voxelize_points(pts_c0, voxel_size=voxel, pad_voxels=2)

        try:
            surf_c0, pre_smooth_c0 = mesh_from_binary_volume(
                mask, origin, voxel,
                gaussian_sigma_voxels=c.mesh_level0_gaussian_sigma,
                smooth_method=c.mesh_smooth_method,
                smooth_iters=c.mesh_level0_smooth_iters,
                taubin_pass_band=c.mesh_taubin_pass_band,
                fill_holes_size=c.mesh_fill_holes_A,
            )
        except ValueError as e:
            print(f"[{self.key}] MC mesh failed for {level_name}: {e}")
            return

        # Transform both meshes to world coordinates
        def _to_world(mesh_c0: pv.PolyData) -> pv.PolyData:
            pts_w = transform_points_from_C0(
                np.asarray(mesh_c0.points, dtype=np.float32), ptc, constr
            ).astype(np.float32)
            m = mesh_c0.copy(deep=True)
            m.points = pts_w
            return m

        pre_smooth_w = _to_world(pre_smooth_c0)
        pre_smooth_path = stage_dir / f"mesh_{level_name}_pre_smooth.ply"
        save_mesh_with_ascii(pre_smooth_w, pre_smooth_path, tag=f"{level_name}-pre-smooth")

        surf_w = _to_world(surf_c0)

        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)
        surf_w = clip_mesh_to_atom_clearance(surf_w, region_xyz, min_clearance_A=c.mesh_atom_clearance_A)

        dt = time.perf_counter() - t0
        # print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
        #       f"{surf_w.n_faces:,} faces, watertight={surf_w.is_manifold and surf_w.n_open_edges == 0}")

        print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
            f"{surf_w.n_faces_strict:,} faces, watertight={surf_w.is_manifold and surf_w.n_open_edges == 0}")

        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        save_mesh_with_ascii(surf_w, mesh_path, tag=level_name)

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name, "method": "marching_cubes_taubin"},
        )
        print(f"[{self.key}] mesh saved: {mesh_path}")

class Stage70MeshValidate(Stage):
    key = "70_mesh_validate"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        import json
        import shutil

        stage_dir = ctx.store.stage_dir(self.key)
        mesh_path = stage_dir / "npet2_tunnel_mesh.ply"

        def _mesh_stats(m: pv.PolyData) -> dict:
            return {
                "n_points": int(m.n_points),

                "n_faces": int(m.n_faces_strict),
                "open_edges": int(m.n_open_edges),
                "is_manifold": bool(m.is_manifold),
                "bounds": [float(x) for x in m.bounds],
            }

        # Try level_1 first (higher detail), fall back to level_0
        chosen_src = None
        chosen_label = None

        l1_path = ctx.inputs.get("level_1_mesh_path")
        if l1_path and Path(l1_path).exists():
            m = pv.read(l1_path)
            if m.is_manifold and m.n_open_edges == 0 and m.n_points > 0:
                chosen_src = l1_path
                chosen_label = "level_1"
                print(f"[{self.key}] using level_1 mesh (0.5A grid)")

        if chosen_src is None:
            # Look for level_0
            l0_path = ctx.store.run_dir / "stage" / "50_clustering" / "mesh_level_0.ply"
            if l0_path.exists():
                m = pv.read(str(l0_path))
                if m.n_points > 0:
                    chosen_src = str(l0_path)
                    chosen_label = "level_0"
                    print(f"[{self.key}] falling back to level_0 mesh (1.0A grid)")

        if chosen_src is None:
            raise ValueError(f"[{self.key}] no valid mesh found from any stage")

        shutil.copy2(chosen_src, mesh_path)
        final = pv.read(str(mesh_path))

        # Save final mesh as both binary and ASCII
        save_mesh_with_ascii(final, mesh_path, tag="final")

        st = _mesh_stats(final)
        watertight = final.is_manifold and final.n_open_edges == 0
        print(f"[{self.key}] final mesh: {st}, watertight={watertight}")

        if not watertight:
            raise ValueError(f"[{self.key}] final mesh is not watertight")

        ctx.store.register_file(
            name="tunnel_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"watertight": True, "source": chosen_label},
        )
        ctx.inputs["tunnel_mesh_path"] = str(mesh_path)

        # Also copy final meshes to run root for convenience
        root_mesh = ctx.store.run_dir / "tunnel_mesh.ply"
        root_mesh_ascii = ctx.store.run_dir / "tunnel_mesh_ascii.ply"
        shutil.copy2(str(mesh_path), str(root_mesh))
        ascii_src = mesh_path.parent / f"{mesh_path.stem}_ascii.ply"
        if ascii_src.exists():
            shutil.copy2(str(ascii_src), str(root_mesh_ascii))

        self._copy_comparison_meshes(ctx, stage_dir)

    def _copy_comparison_meshes(self, ctx, stage_dir):
        import shutil

        run_root = ctx.store.run_dir

        for stage_name, level, voxel in [
            ("50_clustering", "level_0", 1.0),
            ("55_grid_refine", "level_1", 0.5),
        ]:
            src_dir = ctx.store.run_dir / "stage" / stage_name

            # Post-smooth mesh
            for suffix in [f"mesh_{level}.ply", f"mesh_{level}_ascii.ply"]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    # Also put in run root
                    shutil.copy2(src, run_root / suffix)

            # Pre-smooth mesh
            for suffix in [
                f"mesh_{level}_pre_smooth.ply",
                f"mesh_{level}_pre_smooth_ascii.ply",
            ]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    shutil.copy2(src, run_root / suffix)

            if (src_dir / f"mesh_{level}.ply").exists():
                print(
                    f"[{self.key}]   copied {level} mesh ({voxel}A grid) + pre-smooth"
                )

```

libnpet/stages/grid_refine.py
```py
# grid_refine.py

from __future__ import annotations
from libnpet.backends.clustering_io import clusters_from_labels

# npet2/stages/grid_refine.py
from dataclasses import asdict
import json
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pyvista as pv
from scipy import ndimage
import time
import open3d as o3d

from libnpet.core.pipeline import Stage
from libnpet.core.types import StageContext, ArtifactType

from libnpet.backends.geometry import (
    transform_points_to_C0,
    transform_points_from_C0,
    estimate_normals,
)

from libnpet.backends.grid_occupancy import (
    GridSpec,
    occupancy_via_edt,
)

# ... rest of the file is identical, except:
# - In _save_dbscan_pass, replace the dynamic import with:
#     from libnpet.backends.clustering_io import clusters_from_labels
# - In _generate_mesh, replace:
#     from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
#   with nothing (already imported at top)
#   and replace:
#     from libnpet.backends.meshing import ...
#   with:
#     from libnpet.backends.meshing import ...


def _make_bbox_grid(lo: np.ndarray, hi: np.ndarray, voxel: float) -> GridSpec:
    lo = np.asarray(lo, dtype=np.float32)
    hi = np.asarray(hi, dtype=np.float32)
    voxel = float(voxel)

    span = hi - lo
    shape = tuple((np.ceil(span / voxel).astype(np.int32) + 1).tolist())
    return GridSpec(origin=lo, voxel_size=voxel, shape=shape)


def _voxel_centers_from_indices(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin[None, :] + ijk * float(grid.voxel_size)


def _points_to_ijk(grid: GridSpec, pts_c0: np.ndarray) -> np.ndarray:
    """Nearest-voxel mapping for points in C0 -> ijk indices."""
    v = float(grid.voxel_size)
    ijk = np.floor((pts_c0 - grid.origin[None, :]) / v + 0.5).astype(np.int32)
    return ijk


def _valid_ijk(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    nx, ny, nz = grid.shape
    m = (
        (ijk[:, 0] >= 0)
        & (ijk[:, 0] < nx)
        & (ijk[:, 1] >= 0)
        & (ijk[:, 1] < ny)
        & (ijk[:, 2] >= 0)
        & (ijk[:, 2] < nz)
    )
    return m




class Stage55GridRefine(Stage):
    key = "55_grid_refine"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A"         : float(c.refine_voxel_size_A),
            "roi_pad_A"            : float(c.refine_roi_pad_A),
            "atom_radius_A"        : float(c.refine_atom_radius_A),
            "keep_within_A"        : float(c.refine_keep_within_A),
            "occ_close_iters"      : int(c.refine_occ_close_iters),
            "void_open_iters"      : int(c.refine_void_open_iters),
            "forbid_roi_boundary"  : bool(c.refine_forbid_roi_boundary),
            "coarse_eps_A"         : float(c.dbscan_level1_coarse_eps_A),
            "coarse_min_samples"   : int(c.dbscan_level1_coarse_min_samples),
            "refine_eps_A"         : float(c.dbscan_level1_refine_eps_A),
            "refine_min_samples"   : int(c.dbscan_level1_refine_min_samples),
            "dbscan_max_points"    : int(c.refine_dbscan_max_points),
            "dbscan_seed"          : int(c.refine_dbscan_seed),
            "mesh_enable"          : bool(getattr(c, "mesh_level1_enable", True)),
            "mesh_poisson_depth"   : int(getattr(c, "mesh_level1_poisson_depth", 8)),
            "mesh_poisson_ptweight": float(getattr(c, "mesh_level1_poisson_ptweight", 0.5)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        from pathlib import Path

        import numpy as np
        import pyvista as pv
        from scipy import ndimage
        from scipy.spatial import cKDTree
        from sklearn.cluster import DBSCAN

        c = ctx.config
        stage_dir = Path(ctx.store.stage_dir(self.key))
        stage_dir.mkdir(parents=True, exist_ok=True)

        refined_world = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)
        
        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = str(ctx.require("alpha_shell_path"))
        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        voxel = float(c.refine_voxel_size_A)
        pad = float(c.refine_roi_pad_A)
        atom_r = float(c.refine_atom_radius_A)
        keep_within_A = float(c.refine_keep_within_A)
        occ_close_iters = int(c.refine_occ_close_iters)
        void_open_iters = int(c.refine_void_open_iters)
        forbid_roi_boundary = bool(c.refine_forbid_roi_boundary)

        eps_c = float(c.dbscan_level1_coarse_eps_A)
        ms_c = int(c.dbscan_level1_coarse_min_samples)
        eps_r = float(c.dbscan_level1_refine_eps_A)
        ms_r = int(c.dbscan_level1_refine_min_samples)

        dbscan_max_points = int(c.refine_dbscan_max_points)
        dbscan_seed = int(c.refine_dbscan_seed)

        print(f"[{self.key}] voxel={voxel}Å, ROI pad={pad}Å, atom_r={atom_r}Å")

        refined_c0 = transform_points_to_C0(refined_world, ptc, constr).astype(np.float32)
        lo = refined_c0.min(axis=0) - pad
        hi = refined_c0.max(axis=0) + pad

        roi_obj = {
            "roi_id": "bbox_pad_stage50",
            "frame": "C0",
            "pad_A": float(pad),
            "lo": [float(x) for x in lo.tolist()],
            "hi": [float(x) for x in hi.tolist()],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
            "source": {"stage": "50_clustering", "artifact": "refined_cluster"},
        }
        (stage_dir / "roi_bbox_c0.json").write_text(json.dumps(roi_obj, indent=2))

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr).astype(np.float32)
        lo_sel = lo - atom_r
        hi_sel = hi + atom_r
        m_atoms = np.all((region_c0 >= lo_sel[None, :]) & (region_c0 <= hi_sel[None, :]), axis=1)
        atoms_roi_c0 = region_c0[m_atoms]
        
        print(f"[{self.key}] selected {atoms_roi_c0.shape[0]:,} atoms near ROI (ALL atoms, prevents interference)")

        grid = _make_bbox_grid(lo, hi, voxel)
        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))
        z_max = z_min + float(c.cylinder_height_A)

        cyl = self._cylinder_mask_bbox_grid(
            grid,
            radius_A=float(c.cylinder_radius_A),
            zmin_A=z_min,
            zmax_A=z_max,
        )

        occupied = occupancy_via_edt(atoms_roi_c0, grid, atom_radius_A=atom_r)

        if occ_close_iters > 0:
            occupied = ndimage.binary_closing(occupied, iterations=occ_close_iters)

        occupied = occupied | (~cyl)
        np.save(stage_dir / "occupied_mask_level_1.npy", occupied.astype(np.uint8))

        empty_mask = (~occupied) & cyl
        if empty_mask.sum() == 0:
            raise ValueError(f"[{self.key}] no empty voxels in ROI")

        empty_idx = np.argwhere(empty_mask)
        empty_pts_c0 = _voxel_centers_from_indices(grid, empty_idx).astype(np.float32)

        shell_world = pv.read(alpha_shell_path).triangulate()
        shell_c0 = shell_world.copy(deep=True)
        shell_c0.points = transform_points_to_C0(
            np.asarray(shell_world.points, dtype=np.float32), ptc, constr
        )

        sel = pv.PolyData(empty_pts_c0).select_enclosed_points(shell_c0, check_surface=watertight)
        inside_flags = (np.asarray(sel["SelectedPoints"], dtype=np.int8) == 1)

        void_mask = np.zeros_like(empty_mask, dtype=np.bool_)
        inside_idx = empty_idx[inside_flags]
        if inside_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] no empty voxels inside shell")
        
        void_mask[inside_idx[:, 0], inside_idx[:, 1], inside_idx[:, 2]] = True

        if forbid_roi_boundary:
            void_mask[0, :, :]  = False
            void_mask[-1, :, :] = False
            void_mask[:, 0, :]  = False
            void_mask[:, -1, :] = False
            void_mask[:, :, 0]  = False
            void_mask[:, :, -1] = False

        if keep_within_A > 0.0:
            coarse_ijk = _points_to_ijk(grid, refined_c0)
            m_valid = _valid_ijk(grid, coarse_ijk)
            coarse_ijk = coarse_ijk[m_valid]
            if coarse_ijk.shape[0] > 0:
                seed = np.zeros(grid.shape, dtype=np.bool_)
                seed[coarse_ijk[:, 0], coarse_ijk[:, 1], coarse_ijk[:, 2]] = True
                dist_vox = ndimage.distance_transform_edt(~seed)
                r_vox = float(keep_within_A) / float(voxel)
                void_mask = void_mask & (dist_vox <= r_vox)

        if void_open_iters > 0:
            st = ndimage.generate_binary_structure(3, 1)
            void_mask = ndimage.binary_opening(void_mask, structure=st, iterations=void_open_iters)

        np.save(stage_dir / "void_mask_level_1.npy", void_mask.astype(np.uint8))

        st_er = ndimage.generate_binary_structure(3, 1)
        er = ndimage.binary_erosion(void_mask, structure=st_er, iterations=1)
        boundary_mask = void_mask & (~er)
        boundary_idx = np.argwhere(boundary_mask)
        
        if boundary_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] boundary extraction produced 0 voxels")

        boundary_pts_c0 = _voxel_centers_from_indices(grid, boundary_idx).astype(np.float32)
        boundary_pts_w = transform_points_from_C0(boundary_pts_c0, ptc, constr).astype(np.float32)

        print(f"[{self.key}] boundary points: {boundary_pts_w.shape[0]:,}")

        boundary_pts_c0_cap, boundary_pts_w_cap, idx_cap = self._maybe_cap(
            boundary_pts_c0, boundary_pts_w, dbscan_max_points, dbscan_seed
        )

        tree = cKDTree(refined_c0)

        t0 = time.perf_counter()
        db_coarse = DBSCAN(eps=eps_c, min_samples=ms_c, metric="euclidean", n_jobs=-1)
        labels_coarse = db_coarse.fit_predict(boundary_pts_c0_cap).astype(np.int32)
        dt0 = time.perf_counter() - t0

        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        d_coarse = stage_dir / "coarse"
        self._save_dbscan_pass(
            d_coarse,
            boundary_pts_w_cap,
            labels_coarse,
            eps=eps_c,
            min_samples=ms_c,
        )

        stats_coarse = self._cluster_stats(boundary_pts_c0_cap, labels_coarse, tree)
        best_coarse = self._choose_best_label(stats_coarse)
        
        if best_coarse == -1:
            raise ValueError(
                f"[{self.key}] coarse DBSCAN produced no clusters. "
                f"Try eps={eps_c} → {eps_c*1.5:.1f} or min_samples={ms_c} → {ms_c//2}"
            )

        m_best = labels_coarse == best_coarse
        pts_refine_c0 = boundary_pts_c0_cap[m_best]
        pts_refine_w = boundary_pts_w_cap[m_best]

        t1 = time.perf_counter()
        db_refine = DBSCAN(eps=eps_r, min_samples=ms_r, metric="euclidean", n_jobs=-1)
        labels_refine = db_refine.fit_predict(pts_refine_c0).astype(np.int32)
        dt1 = time.perf_counter() - t1

        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        d_refine = stage_dir / "refine"
        self._save_dbscan_pass(
            d_refine,
            pts_refine_w,
            labels_refine,
            eps=eps_r,
            min_samples=ms_r,
        )

        stats_refine = self._cluster_stats(pts_refine_c0, labels_refine, tree)
        best_refine = self._choose_best_label(stats_refine)
        
        if best_refine == -1:
            raise ValueError(
                f"[{self.key}] refine DBSCAN produced no clusters. "
                f"Try eps={eps_r} → {eps_r*1.5:.1f} or min_samples={ms_r} → {ms_r//2}"
            )

        final_mask = labels_refine == best_refine
        final_surface_w = pts_refine_w[final_mask].astype(np.float32)
        
        if final_surface_w.shape[0] == 0:
            raise ValueError(f"[{self.key}] final cluster is empty")

        np.save(stage_dir / "refined_surface_points_level_1.npy", final_surface_w)

        print(f"[{self.key}] creating selected void mask for voxel fallback...")
        
        final_surface_c0 = transform_points_to_C0(final_surface_w, ptc, constr).astype(np.float32)
        final_ijk = _points_to_ijk(grid, final_surface_c0)
        m_valid = _valid_ijk(grid, final_ijk)
        final_ijk = final_ijk[m_valid]
        
        boundary_selected = np.zeros(grid.shape, dtype=np.bool_)
        if final_ijk.shape[0] > 0:
            boundary_selected[final_ijk[:, 0], final_ijk[:, 1], final_ijk[:, 2]] = True
        
        from scipy.ndimage import binary_fill_holes, binary_dilation
        
        boundary_dilated = binary_dilation(boundary_selected, iterations=2)
        selected_mask = binary_fill_holes(boundary_dilated)
        selected_mask = selected_mask & void_mask
        
        selected_mask_path = stage_dir / "selected_void_component_mask_level_1.npy"
        np.save(selected_mask_path, selected_mask.astype(np.uint8))
        
        print(f"[{self.key}]   selected mask: {selected_mask.sum():,} voxels")
        
        ctx.inputs["selected_void_component_mask_level_1_path"] = str(selected_mask_path)
        ctx.inputs["grid_spec_level_1_path"] = str(stage_dir / "grid_spec_level_1.json")

        spec_obj = {
            "frame": "C0",
            "origin": [float(x) for x in grid.origin.tolist()],
            "voxel_size_A": float(grid.voxel_size),
            "shape": [int(x) for x in grid.shape],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
        }
        (stage_dir / "grid_spec_level_1.json").write_text(json.dumps(spec_obj, indent=2))

        diag = {
            "voxel_size_A": voxel,
            "roi_pad_A": pad,
            "atom_radius_A": atom_r,
            "keep_within_A": keep_within_A,
            "boundary_points_total": int(boundary_pts_c0.shape[0]),
            "boundary_points_used_for_dbscan": int(boundary_pts_c0_cap.shape[0]),
            "dbscan_coarse": {
                "eps_A": eps_c,
                "min_samples": ms_c,
                "best_label": int(best_coarse),
                "n_clusters": int(len(stats_coarse)),
                "clusters": stats_coarse[:25],
            },
            "dbscan_refine": {
                "eps_A": eps_r,
                "min_samples": ms_r,
                "best_label": int(best_refine),
                "n_clusters": int(len(stats_refine)),
                "clusters": stats_refine[:25],
            },
            "final_surface_points": int(final_surface_w.shape[0]),
        }
        (stage_dir / "dbscan_diagnostics.json").write_text(json.dumps(diag, indent=2))

        print(
            f"[{self.key}] DBSCAN refined-grid: "
            f"coarse={best_coarse} → refine={best_refine}, final_pts={final_surface_w.shape[0]:,}"
        )

        ctx.inputs["refined_cluster_surface"] = True
        ctx.inputs["refined_cluster"] = final_surface_w
        ctx.artifacts["refined_surface_points_level_1"] = str(stage_dir / "refined_surface_points_level_1.npy")

        if c.mesh_level1_enable:
            self._generate_mesh(ctx, final_surface_w, level_name="level_1")

    def _cylinder_mask_bbox_grid(
        self, grid: GridSpec, radius_A: float, zmin_A: float, zmax_A: float
    ) -> np.ndarray:
        nx, ny, nz = grid.shape
        v = float(grid.voxel_size)
        ox, oy, oz = grid.origin

        x = ox + np.arange(nx, dtype=np.float32) * v
        y = oy + np.arange(ny, dtype=np.float32) * v
        z = oz + np.arange(nz, dtype=np.float32) * v

        X, Y = np.meshgrid(x, y, indexing="ij")
        inside_r = (X * X + Y * Y) <= (radius_A * radius_A)
        inside_z = (z >= zmin_A) & (z <= zmax_A)
        return inside_r[:, :, None] & inside_z[None, None, :]

    def _maybe_cap(
        self, 
        points_c0: np.ndarray, 
        points_w: np.ndarray, 
        cap: int, 
        seed: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = points_c0.shape[0]
        if cap and n > cap:
            rng = np.random.default_rng(seed)
            idx = rng.choice(n, size=int(cap), replace=False)
            idx.sort()
            return points_c0[idx], points_w[idx], idx
        idx = np.arange(n, dtype=np.int64)
        return points_c0, points_w, idx

    def _cluster_stats(
        self, points_c0: np.ndarray, labels: np.ndarray, tree: cKDTree
    ) -> list[dict]:
        stats = []
        for lab in np.unique(labels):
            if lab == -1:
                continue
            m = labels == lab
            pts = points_c0[m]
            if pts.shape[0] == 0:
                continue
            d, _ = tree.query(pts, k=1)
            stats.append({
                "label": int(lab),
                "size": int(pts.shape[0]),
                "median_dist_to_stage50_A": float(np.median(d)),
                "p05_dist_to_stage50_A": float(np.percentile(d, 5)),
                "p95_dist_to_stage50_A": float(np.percentile(d, 95)),
            })
        stats.sort(key=lambda x: (x["median_dist_to_stage50_A"], -x["size"]))
        return stats

    def _choose_best_label(self, stats: list[dict]) -> int:
        return int(stats[0]["label"]) if stats else -1

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        eps: float,
        min_samples: int,
    ) -> None:
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        clusters = clusters_from_labels(pts, labels)
        for cid, cpts in clusters.items():
            if cid == -1:
                continue
            if cpts.shape[0] > 0:
                np.save(pass_dir / f"cluster_id{cid}.npy", cpts.astype(np.float32))

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        import time
        import json
        from libnpet.backends.meshing import (
            mesh_from_binary_volume,
            clip_mesh_to_atom_clearance,
            save_mesh_with_ascii,
        )

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        print(f"[{self.key}] generating mesh for {level_name}...")

        mask_path = stage_dir / "selected_void_component_mask_level_1.npy"
        if not mask_path.exists():
            mask_path = stage_dir / "void_mask_level_1.npy"
        if not mask_path.exists():
            print(f"[{self.key}] no void mask found for {level_name}, skipping mesh")
            return

        mask = np.load(mask_path).astype(bool)
        spec = json.loads((stage_dir / "grid_spec_level_1.json").read_text())
        origin = np.asarray(spec["origin"], dtype=np.float32)
        voxel = float(spec["voxel_size_A"])
        ptc = np.asarray(spec["transform"]["ptc"], dtype=np.float32)
        constr = np.asarray(spec["transform"]["constriction"], dtype=np.float32)

        t0 = time.perf_counter()
        try:
            surf_c0, pre_smooth_c0 = mesh_from_binary_volume(
                mask, origin, voxel,
                gaussian_sigma_voxels=c.mesh_level1_gaussian_sigma,
                smooth_method=c.mesh_smooth_method,
                smooth_iters=c.mesh_level1_smooth_iters,
                taubin_pass_band=c.mesh_taubin_pass_band,
                fill_holes_size=c.mesh_fill_holes_A,
            )
        except ValueError as e:
            print(f"[{self.key}] MC mesh failed for {level_name}: {e}")
            return

        def _to_world(mesh_c0: pv.PolyData) -> pv.PolyData:
            pts_w = transform_points_from_C0(
                np.asarray(mesh_c0.points, dtype=np.float32), ptc, constr
            ).astype(np.float32)
            m = mesh_c0.copy(deep=True)
            m.points = pts_w
            return m

        pre_smooth_w = _to_world(pre_smooth_c0)
        pre_smooth_path = stage_dir / f"mesh_{level_name}_pre_smooth.ply"
        save_mesh_with_ascii(pre_smooth_w, pre_smooth_path, tag=f"{level_name}-pre-smooth")

        surf_w = _to_world(surf_c0)

        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)
        surf_w = clip_mesh_to_atom_clearance(surf_w, region_xyz, min_clearance_A=c.mesh_atom_clearance_A)

        dt = time.perf_counter() - t0
        is_watertight = surf_w.is_manifold and surf_w.n_open_edges == 0
        print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
            f"{surf_w.n_faces_strict:,} faces, ...")

        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        save_mesh_with_ascii(surf_w, mesh_path, tag=level_name)

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name, "method": "marching_cubes_taubin",
                  "watertight": is_watertight, "voxel_size_A": voxel},
        )
        print(f"[{self.key}] mesh saved: {mesh_path}")

        ctx.inputs["level_1_mesh_path"] = str(mesh_path)
        ctx.inputs["level_1_mesh_watertight"] = is_watertight

    def _clip_mesh_to_atoms(
        self,
        mesh: pv.PolyData,
        atom_xyz: np.ndarray,
        min_clearance_A: float = 1.5,
    ) -> pv.PolyData:
        """
        Push any mesh vertices that are closer than min_clearance_A to any atom
        back along the atom->vertex direction until they are at min_clearance_A.
        
        This prevents Poisson/smoothing overshoot from eating into atom-occupied space.
        """
        from scipy.spatial import cKDTree

        tree = cKDTree(atom_xyz)
        pts = np.asarray(mesh.points, dtype=np.float64)

        dist, idx = tree.query(pts, k=1)
        violating = dist < min_clearance_A
        n_violations = int(violating.sum())

        if n_violations == 0:
            print(f"[{self.key}]   atom clearance check: all vertices OK (min_clearance={min_clearance_A}A)")
            return mesh

        nearest_atom = atom_xyz[idx[violating]]
        direction = pts[violating] - nearest_atom
        norms = np.linalg.norm(direction, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-8)
        direction = direction / norms

        pts[violating] = nearest_atom + direction * min_clearance_A

        result = mesh.copy()
        result.points = pts.astype(np.float32)

        print(f"[{self.key}]   atom clearance check: pushed {n_violations:,} vertices "
            f"(of {pts.shape[0]:,}) to min {min_clearance_A}A clearance")

        return result
```

libnpet/stages/bootstrap.py
```py
# ribctl/lib/npet2/stages/bootstrap.py (updated)
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from libnpet.core.pipeline import Stage
from libnpet.core.ribosome_types import RibosomeProfile
from libnpet.core.types import StageContext


class Stage00Inputs(Stage):
    key = "00_inputs"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        structure_provider = ctx.require("structure_provider")
        data = structure_provider.load_atoms(ctx.rcsb_id)

        atom_xyz = np.asarray(data["atom_xyz"], dtype=np.float32)
        ctx.inputs["atom_xyz"] = atom_xyz
        ctx.inputs["atom_element"] = data.get("atom_element", None)
        ctx.inputs["mmcif_path"] = data["mmcif_path"]

        # Profile: validate it's the right type
        profile = data["profile"]
        if not isinstance(profile, RibosomeProfile):
            profile = RibosomeProfile.model_validate(profile)
        ctx.inputs["profile"] = profile

        # If the provider gave us a biopython structure or RibosomeOps, stash it.
        # Standalone providers won't -- Stage30 will parse mmcif on demand.
        if "ro" in data:
            ro = data["ro"]
            ctx.inputs["ro"] = ro
            ctx.inputs["biopython_structure"] = ro.assets.biopython_structure()
        elif "biopython_structure" in data:
            ctx.inputs["biopython_structure"] = data["biopython_structure"]

        ctx.artifacts["atom_xyz"] = ctx.store.put_numpy(
            name="atom_xyz",
            stage=self.key,
            arr=atom_xyz,
            meta={"shape": list(atom_xyz.shape), "dtype": str(atom_xyz.dtype)},
        )

        mins = atom_xyz.min(axis=0)
        maxs = atom_xyz.max(axis=0)
        ctx.stats["atom_bounds"] = {"min": mins.tolist(), "max": maxs.tolist()}
        ctx.stats["n_atoms"] = int(atom_xyz.shape[0])


class Stage10Landmarks(Stage):
    key = "10_landmarks"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        landmark_provider = ctx.require("landmark_provider")
        lm = landmark_provider.get_landmarks(ctx.rcsb_id)

        ptc = np.asarray(lm["ptc_xyz"], dtype=np.float32)
        constr = np.asarray(lm["constriction_xyz"], dtype=np.float32)

        if ptc.shape != (3,):
            raise ValueError(f"PTC must be shape (3,), got {ptc.shape}")
        if constr.shape != (3,):
            raise ValueError(f"Constriction must be shape (3,), got {constr.shape}")

        ctx.inputs["ptc_xyz"] = ptc
        ctx.inputs["constriction_xyz"] = constr

        D = float(np.linalg.norm(constr - ptc))
        if D < 1.0:
            raise ValueError(
                f"PTC and constriction are too close ({D:.1f}A) -- check coordinates"
            )

        z_min = -float(ctx.config.cylinder_ptc_extension_A)
        z_max = float(ctx.config.cylinder_height_A)

        ctx.inputs["cylinder_z_min"] = z_min
        ctx.inputs["cylinder_z_max"] = z_max
        ctx.inputs["landmark_distance"] = D

        print(
            f"  [10_landmarks] PTC-Constriction distance={D:.1f}A, "
            f"cylinder z=[{z_min:.1f}, {z_max:.1f}]A"
        )

        ctx.artifacts["ptc"] = ctx.store.put_json(
            name="ptc",
            stage=self.key,
            obj={"location": ptc.tolist()},
            meta={"units": "A"},
        )
        ctx.artifacts["constriction_site"] = ctx.store.put_json(
            name="constriction_site",
            stage=self.key,
            obj={"location": constr.tolist()},
            meta={
                "units": "A",
                "landmark_distance_A": D,
                "cylinder_z_min": z_min,
                "cylinder_z_max": z_max,
            },
        )

```

libnpet/core/interfaces.py
```py
# ribctl/lib/npet2/core/interfaces.py
"""
Provider protocols for npet2.

These define the boundary between the pipeline and any data source.
Implement these to plug in riboxyz, a local file system, or an API.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Protocol

import numpy as np

from .ribosome_types import RibosomeProfile, PTCInfo, ConstrictionInfo
from .types import ArtifactRef, ArtifactType


class StructureProvider(Protocol):
    """Provides atom coordinates and the ribosome profile for a structure."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        """
        Must return:
          - atom_xyz: (N, 3) float32
          - atom_element: (N,) str array (optional)
          - mmcif_path: str
          - profile: RibosomeProfile
        """
        ...


class LandmarkProvider(Protocol):
    """Provides PTC and constriction site coordinates."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        """
        Must return:
          - ptc_xyz: (3,) float32
          - constriction_xyz: (3,) float32
        """
        ...


class ArtifactStore(Protocol):
    @property
    def run_dir(self) -> Path: ...

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType,
                  data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_json(self, *, name: str, stage: str, obj: Any,
                 meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray,
                  meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def add_ref(self, ref: ArtifactRef) -> None: ...

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None: ...

```

libnpet/core/pipeline.py
```py
# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
import sys
from typing import Any, Dict, List
import traceback as tb
import time

from .types import StageContext


class Stage(ABC):
    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]: ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None: ...


def _fmt_value(v: Any) -> str:
    """Compact display for a parameter value."""
    if isinstance(v, float):
        # drop trailing zeros but keep one decimal
        return f"{v:g}"
    if isinstance(v, list) and len(v) > 3:
        return f"[{len(v)} items]"
    return str(v)


def _log_params(params: Dict[str, Any], prefix: str) -> None:
    if not params:
        return
    items = [f"{k}={_fmt_value(v)}" for k, v in params.items()
             if not isinstance(v, (dict, list))]
    # nested dicts/lists get their own lines
    nested = {k: v for k, v in params.items() if isinstance(v, (dict, list))}

    if items:
        line = ", ".join(items)
        # wrap at ~100 chars
        if len(line) > 100:
            mid = len(items) // 2
            print(f"  [{prefix}] {', '.join(items[:mid])}")
            print(f"  [{prefix}] {', '.join(items[mid:])}")
        else:
            print(f"  [{prefix}] {line}")
    for k, v in nested.items():
        if isinstance(v, list) and all(isinstance(x, dict) for x in v):
            for i, entry in enumerate(v):
                sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in entry.items())
                print(f"  [{prefix}] {k}[{i}]: {sub}")
        elif isinstance(v, dict):
            sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in v.items())
            print(f"  [{prefix}] {k}: {sub}")


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        n = len(self.stages)
        wall = 60

        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | run {ctx.run_id}")
        print(f"  stages: {n} | config: {type(ctx.config).__name__}")
        print("=" * wall)

        t_total = time.perf_counter()

        for i, stage in enumerate(self.stages, 1):
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)

            print()
            print(f"--- [{i}/{n}] {stage.key} " + "-" * max(0, wall - len(stage.key) - 12))
            _log_params(params, stage.key)

            t0 = time.perf_counter()
            try:
                stage.run(ctx)
                dt = time.perf_counter() - t0
                print(f"  [{stage.key}] done in {dt:,.2f}s")
                ctx.store.end_stage(stage.key, success=True, note=f"elapsed_s={dt:.3f}")
            except Exception as e:
                dt = time.perf_counter() - t0
                sys.stdout.flush()
                print(f"  [{stage.key}] FAILED after {dt:,.2f}s: {e}", flush=True)
                tb.print_exc()
                sys.stdout.flush()
                ctx.store.end_stage(
                    stage.key, success=False, note=f"elapsed_s={dt:.3f} err={e}"
                )
                ctx.store.finalize(success=False, error=str(e))
                raise

        dt_total = time.perf_counter() - t_total
        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | completed in {dt_total:,.2f}s")
        print(f"  run_dir: {ctx.store.run_dir}")
        print("=" * wall)
        print()

        ctx.store.finalize(success=True)
        return ctx
```

libnpet/core/manifest.py
```py
# ribctl/lib/npet2/core/manifest.py
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import time

from .types import ArtifactRef, ArtifactType


@dataclass
class StageRecord:
    name: str
    status: str = "pending"  # pending|running|success|failure|skipped
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    params: Dict[str, Any] = field(default_factory=dict)
    note: Optional[str] = None


@dataclass
class RunManifest:
    rcsb_id: str
    run_id: str
    pipeline_version: str
    created_at: float = field(default_factory=lambda: time.time())

    inputs: Dict[str, Any] = field(default_factory=dict)
    config_resolved: Dict[str, Any] = field(default_factory=dict)

    stages: Dict[str, StageRecord] = field(default_factory=dict)
    artifacts: List[Dict[str, Any]] = field(default_factory=list)

    success: Optional[bool] = None
    error: Optional[str] = None

    def add_artifact(self, ref: ArtifactRef) -> None:
        self.artifacts.append({
            "name": ref.name,
            "type": ref.type.value,
            "path": str(ref.path),
            "stage": ref.stage,
            "meta": ref.meta,
            "depends_on": list(ref.depends_on),
        })

    def to_json(self) -> str:
        # dataclasses → dict
        d = asdict(self)
        # StageRecord needs manual flatten
        d["stages"] = {k: asdict(v) for k, v in self.stages.items()}
        return json.dumps(d, indent=2)

    @staticmethod
    def from_path(path: Path) -> "RunManifest":
        data = json.loads(path.read_text())
        m = RunManifest(
            rcsb_id=data["rcsb_id"],
            run_id=data["run_id"],
            pipeline_version=data["pipeline_version"],
            created_at=data.get("created_at", time.time()),
            inputs=data.get("inputs", {}),
            config_resolved=data.get("config_resolved", {}),
        )
        m.success = data.get("success")
        m.error = data.get("error")
        # stages
        for k, v in data.get("stages", {}).items():
            m.stages[k] = StageRecord(**v)
        m.artifacts = data.get("artifacts", [])
        return m

```

libnpet/adapters/standalone_providers.py
```py
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
```

libnpet/backends/clustering_io.py
```py
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

def clusters_from_labels(points: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray]:
    clusters: dict[int, list[int]] = {}
    for i, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(i)
    out = {}
    for lab, idxs in clusters.items():
        out[lab] = points[np.asarray(idxs, dtype=np.int32)]
    return out

def write_dbscan_pass(out_dir: Path, *, prefix: str, points: np.ndarray, labels: np.ndarray) -> dict:
    out_dir = out_dir / prefix
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(out_dir / "points.npy", points.astype(np.float32))
    np.save(out_dir / "labels.npy", labels.astype(np.int32))

    clusters = clusters_from_labels(points, labels)
    index = {"prefix": prefix, "n_points": int(points.shape[0]), "clusters": []}

    for cid, pts in sorted(clusters.items(), key=lambda kv: kv[0]):
        p = out_dir / f"cluster_id{cid}.npy"
        np.save(p, pts.astype(np.float32))
        index["clusters"].append({"id": int(cid), "n": int(pts.shape[0]), "path": p.name})

    (out_dir / "index.json").write_text(json.dumps(index, indent=2))
    return index

```

libnpet/backends/geometry.py
```py
# npet2/backends/geometry.py
"""
Geometric primitives ported from riboxyz's kdtree_approach.py and alphalib.py.

Only the functions that npet2 actually uses are included here.
"""

from __future__ import annotations

import multiprocessing
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable, List, Optional, Tuple, TypeVar

import numpy as np
import open3d as o3d
import plyfile
import pyvista as pv
from Bio.PDB.MMCIFParser import MMCIFParser
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN as skDBSCAN

from libnpet.core.config import SETTINGS

T = TypeVar("T")


# ---------------------------------------------------------------------------
# Coordinate transforms (PTC/constriction axis -> canonical C0 frame)
# ---------------------------------------------------------------------------


def get_transformation_to_C0(
    base_point: np.ndarray, axis_point: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    z_axis = np.array([0, 0, 1])

    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])
    else:
        v = np.cross(axis_unit, z_axis)
        s = np.linalg.norm(v)
        c = np.dot(axis_unit, z_axis)
        v_skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + v_skew + (v_skew @ v_skew) * (1 - c) / (s * s)

    return -base_point, R


def transform_points_to_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    return (points + translation) @ rotation.T


def transform_points_from_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    return (points @ rotation) - translation


# ---------------------------------------------------------------------------
# DBSCAN wrappers
# ---------------------------------------------------------------------------


def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps: float,
    min_samples: int,
    metric: str = "euclidean",
) -> Tuple[skDBSCAN, dict]:
    print(
        f"Running DBSCAN on {len(ptcloud)} points. "
        f"eps={eps}, min_samples={min_samples}, distance_metric={metric}"
    )
    db = skDBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(ptcloud)

    clusters: dict[int, list] = {}
    for point, label in zip(ptcloud, db.labels_):
        clusters.setdefault(int(label), []).append(point)

    return db, dict(sorted(clusters.items()))


def DBSCAN_pick_largest_cluster(
    clusters_container: dict[int, list],
) -> Tuple[np.ndarray, int]:
    best_id = 0
    for k, v in clusters_container.items():
        if int(k) == -1:
            continue
        if len(v) > len(clusters_container.get(best_id, [])):
            best_id = int(k)
    return np.array(clusters_container[best_id]), best_id


# ---------------------------------------------------------------------------
# Poisson reconstruction (calls external binary)
# ---------------------------------------------------------------------------


def apply_poisson_reconstruction(
    surf_estimated_ptcloud_path: str,
    output_path: Path,
    recon_depth: int = 6,
    recon_pt_weight: float = 3.0,
) -> None:
    print(f"Open3D Poisson Reconstruction: {surf_estimated_ptcloud_path} -> {output_path}")

    pcd = o3d.io.read_point_cloud(str(surf_estimated_ptcloud_path))
    if not pcd.has_normals():
        raise ValueError(f"Point cloud at {surf_estimated_ptcloud_path} has no normals")

    mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
        pcd, depth=recon_depth
    )

    # Trim low-density vertices; these tend to be spurious extrapolations
    # at the boundary where point coverage is thin. 5% threshold is generous
    # enough for a bounding shell.
    densities = np.asarray(densities)
    verts_to_remove = densities < np.quantile(densities, 0.05)
    mesh.remove_vertices_by_mask(verts_to_remove)

    o3d.io.write_triangle_mesh(str(output_path), mesh)
    print(f">>Wrote {output_path}")


# ---------------------------------------------------------------------------
# Cylinder filtering
# ---------------------------------------------------------------------------


def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    z_min: float = 0.0,
) -> bool:
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    point_vector = point - base_point
    projection = np.dot(point_vector, axis_unit)

    projection_point = base_point + projection * axis_unit
    radial_distance = np.linalg.norm(point - projection_point)

    return (radial_distance <= radius) and (z_min <= projection <= height)


def _worker_process_chunk(chunk_data):
    positions, base_point, axis_point, radius, height, indices, z_min = chunk_data
    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(
            pos, base_point, axis_point, radius, height, z_min=z_min
        ):
            results.append(indices[i])
    return results


def get_residue_position(residue):
    return residue.center_of_mass()


def filter_residues_parallel(
    residues: list,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None,
    z_min: float = 0.0,
) -> list:
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    if chunk_size is None:
        chunk_size = max(1, len(residues) // (max_workers * 4))

    positions = np.array([get_residue_position(r) for r in residues])
    indices = list(range(len(residues)))
    index_chunks = [
        indices[i : i + chunk_size] for i in range(0, len(indices), chunk_size)
    ]

    chunks_data = [
        (positions[idx], base_point, axis_point, radius, height, idx, z_min)
        for idx in index_chunks
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))

    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]


# ---------------------------------------------------------------------------
# Grid voxel mask (legacy kdtree backend)
# ---------------------------------------------------------------------------


def generate_voxel_centers(
    radius: float, height: float, voxel_size: float, z_min: float = 0.0
) -> Tuple[np.ndarray, tuple]:
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(z_min, z_min + height, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)


def create_point_cloud_mask(
    points: np.ndarray,
    radius: float,
    height: float,
    voxel_size: float = 1.0,
    radius_around_point: float = 2.0,
    z_min: float = 0.0,
):
    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(
        radius, height, voxel_size, z_min=z_min
    )
    tree = cKDTree(points)
    indices = tree.query_ball_point(voxel_centers, radius_around_point)

    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True
    point_cloud_mask = point_cloud_mask.reshape(grid_shape)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    cylinder_mask = np.sqrt(X**2 + Y**2) <= radius
    hollow_cylinder = ~cylinder_mask

    final_mask = hollow_cylinder | point_cloud_mask
    return final_mask, (x, y, z)


# ---------------------------------------------------------------------------
# Normal estimation
# ---------------------------------------------------------------------------


def estimate_normals(
    convex_hull_surface_pts: np.ndarray,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
) -> o3d.geometry.PointCloud:
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    return pcd


# ---------------------------------------------------------------------------
# Alpha shape / surface extraction (from alphalib)
# ---------------------------------------------------------------------------


def cif_to_point_cloud(
    cif_path: str,
    chains: list[str] | None = None,
    do_atoms: bool = False,
    exclude_chains: list[str] = [],
) -> np.ndarray:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    coordinates = []

    first_model = structure[0]
    if do_atoms:
        for chain in first_model:
            if chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())
    else:
        for chain in first_model:
            if chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                coordinates.append(residue.center_of_mass())

    if not coordinates:
        raise ValueError(f"No coordinates found in {cif_path}")

    return np.array(coordinates)


def quick_surface_points(
    pointcloud: np.ndarray,
    alpha: float,
    tolerance: float,
    offset: float,
    max_points: int = 80_000,
) -> np.ndarray:
    pts = pointcloud
    if len(pts) > max_points:
        idx = np.random.choice(len(pts), max_points, replace=False)
        pts = pts[idx]
        print(
            f"  [quick_surface_points] subsampled {len(pointcloud):,} -> {len(pts):,} points"
        )

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pts.astype(np.float64))

    # open3d alpha shape is stable on large inputs unlike VTK's Delaunay3D
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(
        pcd, alpha=alpha
    )
    surface_pts = np.asarray(mesh.vertices, dtype=np.float32)

    if surface_pts.shape[0] == 0:
        print(
            f"  [quick_surface_points] alpha={alpha} gave empty surface, falling back to input points"
        )
        return pts.astype(np.float32)

    print(
        f"  [quick_surface_points] surface: {surface_pts.shape[0]:,} pts from alpha shape"
    )
    return surface_pts


def fast_normal_estimation(
    surface_pts: np.ndarray,
    kdtree_radius: float,
    max_nn: int,
    tangent_planes_k: int,
) -> o3d.geometry.PointCloud:
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(surface_pts)
    search_param = o3d.geometry.KDTreeSearchParamHybrid(
        radius=kdtree_radius, max_nn=max_nn
    )
    pcd.estimate_normals(search_param=search_param)
    pcd.orient_normals_consistent_tangent_plane(k=tangent_planes_k)
    return pcd


def validate_mesh_pyvista(mesh_or_path, stage="unknown") -> bool:
    if isinstance(mesh_or_path, (str, Path)):
        import os

        if not os.path.exists(mesh_or_path):
            print(f"WARNING: Mesh file does not exist: {mesh_or_path}")
            return False
        try:
            mesh = pv.read(str(mesh_or_path))
        except Exception as e:
            print(f"WARNING: Failed to load mesh at {mesh_or_path}: {e}")
            return False
    else:
        mesh = mesh_or_path

    if mesh is None:
        print(f"WARNING: Null mesh at stage {stage}")
        return False

    try:
        edges = mesh.extract_feature_edges(
            boundary_edges=True,
            feature_edges=False,
            manifold_edges=False,
            non_manifold_edges=False,
        )
        return edges.n_cells == 0
    except Exception as e:
        print(f"WARNING: Error checking watertightness: {e}")
        return False

```

libnpet/backends/grid_occupancy.py
```py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class GridSpec:
    # origin and voxel size define world coords: x = origin + i*voxel
    origin: np.ndarray          # (3,)
    voxel_size: float
    shape: Tuple[int, int, int] # (nx, ny, nz)


# ribctl/lib/npet2/backends/grid_occupancy.py

def make_cylinder_grid(radius_A: float, height_A: float, voxel_A: float,
                       z_min: float = 0.0) -> GridSpec:
    """
    Canonical cylinder in C0:
      x in [-R, R], y in [-R, R], z in [z_min, z_min + H]
    """
    nx = int(np.floor((2 * radius_A) / voxel_A)) + 1
    ny = int(np.floor((2 * radius_A) / voxel_A)) + 1
    nz = int(np.floor(height_A / voxel_A)) + 1
    origin = np.array([-radius_A, -radius_A, z_min], dtype=np.float32)
    return GridSpec(origin=origin, voxel_size=float(voxel_A), shape=(nx, ny, nz))


def grid_world_coords(grid: GridSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns axis coordinate arrays (x, y, z) for voxel centers along each axis.
    """
    ox, oy, oz = grid.origin
    nx, ny, nz = grid.shape
    v = grid.voxel_size
    x = ox + np.arange(nx, dtype=np.float32) * v
    y = oy + np.arange(ny, dtype=np.float32) * v
    z = oz + np.arange(nz, dtype=np.float32) * v
    return x, y, z


def cylinder_mask(grid: GridSpec, radius_A: float) -> np.ndarray:
    """
    Boolean mask of voxels inside cylinder radius (in C0).
    """
    x, y, _ = grid_world_coords(grid)
    X, Y = np.meshgrid(x, y, indexing="ij")
    inside = (X * X + Y * Y) <= (radius_A * radius_A)
    # broadcast across z
    return inside[:, :, None]


def points_to_occupied_seeds(points_c0: np.ndarray, grid: GridSpec) -> np.ndarray:
    """
    Convert points (C0 coords) to a sparse occupied grid of seeds at nearest voxels.
    """
    pts = np.asarray(points_c0, dtype=np.float32)
    v = grid.voxel_size
    origin = grid.origin

    ijk = np.floor((pts - origin[None, :]) / v + 0.5).astype(np.int32)
    nx, ny, nz = grid.shape

    valid = (
        (ijk[:, 0] >= 0) & (ijk[:, 0] < nx) &
        (ijk[:, 1] >= 0) & (ijk[:, 1] < ny) &
        (ijk[:, 2] >= 0) & (ijk[:, 2] < nz)
    )
    ijk = ijk[valid]
    occ = np.zeros(grid.shape, dtype=np.bool_)
    if ijk.shape[0] > 0:
        occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return occ


def occupancy_via_edt(points_c0: np.ndarray, grid: GridSpec, atom_radius_A: float) -> np.ndarray:
    """
    Occupancy grid: voxel is occupied if within atom_radius_A of any atom center.

    Steps:
      - seed occupied at nearest voxels
      - edt on ~occupied gives distance (in voxels) to nearest seed
      - threshold <= r_vox
    """
    seeds = points_to_occupied_seeds(points_c0, grid)

    # Distance (in voxels) from each voxel to nearest True in 'seeds'
    # distance_transform_edt computes distance to nearest zero;
    # so we compute on ~seeds, where zeros correspond to seeds.
    dist_vox = ndimage.distance_transform_edt(~seeds)

    r_vox = float(atom_radius_A) / float(grid.voxel_size)
    occupied = dist_vox <= r_vox
    return occupied


def empty_points_from_mask(grid: GridSpec, empty_mask: np.ndarray) -> np.ndarray:
    """
    Return coordinates (C0) of voxel centers where empty_mask is True.
    """
    empty_idx = np.where(empty_mask)
    x, y, z = grid_world_coords(grid)
    pts = np.column_stack((x[empty_idx[0]], y[empty_idx[1]], z[empty_idx[2]])).astype(np.float32)
    return pts


def save_grid_npy(grid: GridSpec, data: np.ndarray, path: Path, *, compress: bool = False) -> None:
    """
    Save a 3D grid along with its GridSpec metadata.
    File format: {path}_data.npy + {path}_spec.json
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if compress:
        np.savez_compressed(data_path.with_suffix('.npz'), data=data)
    else:
        np.save(data_path, data)
    
    spec_dict = {
        "origin": grid.origin.tolist(),
        "voxel_size": float(grid.voxel_size),
        "shape": list(grid.shape),
    }
    spec_path.write_text(__import__('json').dumps(spec_dict, indent=2))


def load_grid_npy(path: Path) -> tuple[GridSpec, np.ndarray]:
    """
    Load a grid saved by save_grid_npy.
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if data_path.with_suffix('.npz').exists():
        data = np.load(data_path.with_suffix('.npz'))['data']
    else:
        data = np.load(data_path)
    
    spec_dict = __import__('json').loads(spec_path.read_text())
    grid = GridSpec(
        origin=np.array(spec_dict["origin"], dtype=np.float32),
        voxel_size=float(spec_dict["voxel_size"]),
        shape=tuple(spec_dict["shape"]),
    )
    return grid, data


def voxel_to_world(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    """
    Convert voxel indices (i,j,k) to world coordinates.
    ijk: (N, 3) or (3,) array of voxel indices
    Returns: (N, 3) or (3,) world coordinates
    """
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin + ijk * grid.voxel_size


def world_to_voxel(grid: GridSpec, xyz: np.ndarray) -> np.ndarray:
    """
    Convert world coordinates to voxel indices.
    xyz: (N, 3) or (3,) world coordinates
    Returns: (N, 3) or (3,) voxel indices (floats; use floor/round as needed)
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    return (xyz - grid.origin) / grid.voxel_size


def get_occupied_voxel_centers(grid: GridSpec, occupancy: np.ndarray) -> np.ndarray:
    """
    Get world coordinates of occupied voxel centers.
    """
    occupied_idx = np.argwhere(occupancy)
    return voxel_to_world(grid, occupied_idx)
# ribctl/lib/npet2/backends/grid_occupancy.py
# ADD these functions at the end:

from scipy import ndimage

def connected_components_3d(
    binary_mask: np.ndarray, 
    connectivity: int = 26
) -> tuple[np.ndarray, int]:
    """
    Find connected components in a 3D binary mask.
    
    Args:
        binary_mask: 3D boolean array
        connectivity: 6 (face), 18 (face+edge), or 26 (face+edge+corner)
    
    Returns:
        labeled: Array same shape as input with component labels (0=background)
        n_components: Number of components found
    """
    if connectivity == 6:
        structure = ndimage.generate_binary_structure(3, 1)
    elif connectivity == 18:
        structure = ndimage.generate_binary_structure(3, 2)
    elif connectivity == 26:
        structure = ndimage.generate_binary_structure(3, 3)
    else:
        raise ValueError(f"connectivity must be 6, 18, or 26, got {connectivity}")
    
    labeled, n_components = ndimage.label(binary_mask, structure=structure)
    return labeled, n_components


def get_largest_component(labeled: np.ndarray, n_components: int) -> np.ndarray:
    """
    Extract mask of the largest connected component.
    
    Args:
        labeled: Output from connected_components_3d
        n_components: Number of components
    
    Returns:
        Binary mask of largest component only
    """
    if n_components == 0:
        return np.zeros_like(labeled, dtype=bool)
    
    # Count voxels in each component (excluding background=0)
    component_sizes = np.bincount(labeled.ravel())
    component_sizes[0] = 0  # Ignore background
    
    largest_label = np.argmax(component_sizes)
    return labeled == largest_label


def get_component_stats(labeled: np.ndarray, n_components: int) -> list[dict]:
    """
    Get statistics for all connected components.
    
    Returns:
        List of dicts with {label, size, bbox_min, bbox_max}
    """
    stats = []
    
    for label in range(1, n_components + 1):
        mask = labeled == label
        size = int(mask.sum())
        
        if size == 0:
            continue
        
        indices = np.argwhere(mask)
        bbox_min = indices.min(axis=0)
        bbox_max = indices.max(axis=0)
        
        stats.append({
            "label": int(label),
            "size": size,
            "bbox_min": bbox_min.tolist(),
            "bbox_max": bbox_max.tolist(),
        })
    
    # Sort by size descending
    stats.sort(key=lambda x: x["size"], reverse=True)
    return stats


def morphological_clean(
    binary_mask: np.ndarray, 
    operation: str = "opening",
    iterations: int = 1
) -> np.ndarray:
    """
    Apply morphological operations to clean up a binary mask.
    
    Args:
        binary_mask: 3D boolean array
        operation: "opening" (remove small bits), "closing" (fill small holes), 
                  "erosion", "dilation"
        iterations: Number of times to apply operation
    
    Returns:
        Cleaned binary mask
    """
    if operation == "opening":
        return ndimage.binary_opening(binary_mask, iterations=iterations)
    elif operation == "closing":
        return ndimage.binary_closing(binary_mask, iterations=iterations)
    elif operation == "erosion":
        return ndimage.binary_erosion(binary_mask, iterations=iterations)
    elif operation == "dilation":
        return ndimage.binary_dilation(binary_mask, iterations=iterations)
    else:
        raise ValueError(f"Unknown operation: {operation}")
```

libnpet/backends/meshing.py
```py
# ribctl/lib/npet2/backends/meshing.py
from __future__ import annotations
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.ndimage import binary_fill_holes, gaussian_filter

def save_mesh_with_ascii(mesh: pv.PolyData, path: Path, tag: str = "") -> None:
    """Save mesh as binary PLY + ASCII PLY side by side."""
    mesh.save(str(path))
    ascii_path = path.parent / f"{path.stem}_ascii.ply"
    try:
        mesh.save(str(ascii_path), binary=False)
    except Exception:
        try:
            import plyfile
            data = plyfile.PlyData.read(str(path))
            data.text = True
            data.write(str(ascii_path))
        except Exception as e:
            print(f"[meshing] ASCII PLY write failed{' (' + tag + ')' if tag else ''}: {e}")

def mesh_from_binary_volume(
    mask: np.ndarray,
    origin: np.ndarray,
    voxel_size: float,
    *,
    gaussian_sigma_voxels: float = 1.5,
    smooth_method: str = "taubin",
    smooth_iters: int = 20,
    taubin_pass_band: float = 0.1,
    fill_holes_size: float = 100.0,
) -> tuple[pv.PolyData, pv.PolyData]:
    """
    Marching cubes on a Gaussian-blurred binary volume, followed by mesh smoothing.

    Returns (smoothed_mesh, pre_smooth_mesh).
    Both are in the same coordinate frame as the input volume.
    Caller is responsible for coordinate transforms and saving.
    """
    vol = np.pad(mask.astype(np.float32), 2, constant_values=0.0)
    origin_pad = np.asarray(origin, dtype=np.float32) - 2 * voxel_size

    if gaussian_sigma_voxels > 0:
        vol = gaussian_filter(vol, sigma=gaussian_sigma_voxels)

    img = pv.ImageData(
        dimensions=vol.shape,
        spacing=(voxel_size, voxel_size, voxel_size),
        origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
    )
    img.point_data["values"] = vol.ravel(order="F")

    surf = img.contour(isosurfaces=[0.5], scalars="values").triangulate()
    if surf.n_points == 0:
        raise ValueError("Marching cubes produced empty surface")

    surf = surf.clean(tolerance=0.0)

    if fill_holes_size > 0:
        surf = surf.fill_holes(fill_holes_size)

    surf = surf.connectivity(extraction_mode='largest')

    pre_smooth = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)

    if smooth_iters > 0:
        if smooth_method == "taubin":
            surf = surf.smooth_taubin(n_iter=smooth_iters, pass_band=taubin_pass_band)
        else:
            surf = surf.smooth(n_iter=smooth_iters)

    surf = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)
    return surf, pre_smooth


def voxelize_points(
    points: np.ndarray,
    voxel_size: float,
    pad_voxels: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    lo = points.min(axis=0) - pad_voxels * voxel_size
    hi = points.max(axis=0) + pad_voxels * voxel_size

    shape = tuple(np.ceil((hi - lo) / voxel_size).astype(int) + 1)

    ijk = np.floor((points - lo) / voxel_size + 0.5).astype(np.int32)
    for d in range(3):
        ijk[:, d] = np.clip(ijk[:, d], 0, shape[d] - 1)

    mask = np.zeros(shape, dtype=bool)
    mask[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True

    mask = binary_fill_holes(mask)
    return mask, lo.astype(np.float32)


def clip_mesh_to_atom_clearance(
    mesh: pv.PolyData,
    atom_xyz: np.ndarray,
    min_clearance_A: float = 1.5,
) -> pv.PolyData:
    from scipy.spatial import cKDTree

    tree = cKDTree(atom_xyz)
    pts = np.asarray(mesh.points, dtype=np.float64)

    dist, idx = tree.query(pts, k=1)
    violating = dist < min_clearance_A

    if violating.sum() == 0:
        return mesh

    nearest = atom_xyz[idx[violating]]
    direction = pts[violating] - nearest
    norms = np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-8)
    pts[violating] = nearest + (direction / norms) * min_clearance_A

    result = mesh.copy()
    result.points = pts.astype(np.float32)
    return result
```

libnpet/backends/__init__.py
```py

```


