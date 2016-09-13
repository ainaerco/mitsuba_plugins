# mitsuba_plugins
Plugins for mitsuba renderer 0.5.0.

* shapes/subd           - Alembic procedural with OpenSubdiv support
* bsdfs/velvetcoating   - Inverted gaussian bsdf for sheen, coating material
* bsdfs/wardaniso       - Ward with anisotropy rotation input
* bsdfs/zmask           - Modified orginal mask bsdf to render leaves with opacity
* textures/adjust       - Adjust texture color
* textures/blend        - Blend two textures
* textures/checker      - Modified original checker, accepts texture inputs
* textures/combine      - Multiply/Divide/Add/Subtract two textures
* textures/oiio_texture - OpenimageIO texture

Dependencies:
Alembic,
boost,
HDF5,
zlib,
OpenEXR,
OpenSubdiv-3.0.3,
OpenimageIO,
mitsuba-0.5.0
