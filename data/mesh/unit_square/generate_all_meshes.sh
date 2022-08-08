#/bin/bash -eu
# this script creates mesh files which can be used in ParMooN. You need the
# the software 'gmsh' installed and the GEO_FILE below.

GEO_FILE=unit_square.geo
# -2 is the option to generate a 2d-grid
GMSH=gmsh
GMSH_OPTIONS="-2 -format mesh"
TMP_FILE=temporary_geo_file.geo
# first, increment, last. Passed to 'seq'. Change to, eg., '1 1 10' to get many 
# grids
FINENESS='1 5 6'

# check if gmsh is installed, see https://stackoverflow.com/a/677212
command -v $GMSH >/dev/null 2>&1 || \
  { echo >&2 "I require $GMSH but it's not installed.  Aborting."; exit 1; }

# switch to triangles, if not already set  (Mesh -> //Mesh)
sed -i 's#^Mesh#//Mesh#' $GEO_FILE

for i in $(seq -s' ' $FINENESS)
do
  # set fineness of grid
  sed "s#lc = h#lc = h/($i)#" $GEO_FILE > $TMP_FILE
  # name of resulting mesh file
  MESH_FILE="unit_square_tria$i.mesh"
  # call gmsh, create a 2d grid in mesh format
  $GMSH $GMSH_OPTIONS -o $MESH_FILE $TMP_FILE > /dev/null
done

# switch to quads  (//Mesh -> Mesh)
sed -i 's#^//Mesh#Mesh#' $GEO_FILE
for i in $(seq -s' ' $FINENESS)
do
  # set fineness of grid
  sed "s#lc = h#lc = h/($i)#" $GEO_FILE > $TMP_FILE
  # name of resulting mesh file
  MESH_FILE="unit_square_quad$i.mesh"
  # call gmsh, create a 2d grid in mesh format
  $GMSH $GMSH_OPTIONS -o $MESH_FILE $TMP_FILE > /dev/null
done

rm $TMP_FILE
