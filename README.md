# build once
docker build -t npet2 .

# pure API mode — data lands in ./npet2_data on the host
docker run -v $(pwd)/npet2_data:/data npet2 run 5NWY

# file mode — mount a local data directory read-only for inputs
docker run \
  -v $(pwd)/npet2_data:/data \
  -v /path/to/your/structures:/structures:ro \
  npet2 run 5NWY \
    --mmcif /structures/5NWY.cif \
    --profile /structures/5NWY.json \
    --ptc /structures/5NWY_PTC.json \
    --constriction /structures/5NWY_CONSTRICTION_SITE.json

# pointing at your local riboxyz backend
docker run \
  -v $(pwd)/npet2_data:/data \
  -e NPET2_RIBOXYZ_API_URL=http://host.docker.internal:8000 \
  npet2 run 5NWY

# show config
docker run npet2 show-config