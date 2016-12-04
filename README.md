# mitsuba_plugins
Plugins for mitsuba renderer 0.5.0.

* shapes/subd           - Alembic procedural with OpenSubdiv support
* bsdfs/wardaniso       - Ward with anisotropy rotation input
* bsdfs/zmask           - Modified orginal mask bsdf to render leaves with opacity
* textures/adjust       - Adjust texture color
* textures/blend        - Blend two textures
* textures/checker      - Modified original checker, accepts texture inputs
* textures/combine      - Multiply/Divide/Add/Subtract two textures
* textures/colorramp    - Color ramp, currently supports only v direction
* textures/oiio_texture - OpenimageIO texture with udim support
* textures/osl          - Basic OSL texture. Just for procedural textures, no closures, no tracing
* textures/variation    - Randomize color seeded by shape name or instance id

Dependencies:
Alembic,
boost,
HDF5,
zlib,
OpenEXR,
OpenSubdiv-3.0.3,
OpenimageIO,
mitsuba-0.5.0
