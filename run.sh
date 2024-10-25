#!/bin/bash

# Example Run
# ./run.sh --ignitions_file ./data/yosemite_sample_ignitions.csv --rx_burn_units_file ./data/rx_burn.tar.gz --toa_file ./data/toa.tar.gz --flamelen_file ./data/flamelen.tar.gz --ms_buildings_tiles_file ./data/MS_bldgs_all_tiles.json --critical_habitats_polygon_file ./data/critical_habitat_polygons.tar.gz --footprints_polygon_file ./footprints_polygons_2323 --burned_area_file ./burned_area_23423 --building_damage_file ./building_damage_35454 --habitat_damage_file ./habitat_damage_554353

# inputs
IGNITIONS_FILE=""
RX_BURN_UNITS_FILE=""
TOA_FILE=""
INTENSITY_FILE=""
MS_BUILDINGS_TILES_FILE=""
CRITICAL_HABITATS_POLYGON_FILE=""

# outputs
BURNED_AREA_FILE=""
BUILDING_DAMAGE_FILE=""
HABITAT_DAMAGE_FILE=""
FOOTPRINTS_POLYGONS_FILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -i | --ignitions_file)
        echo "Processing 'ignitions_file' option. Input argument is '$2'"
        IGNITIONS_FILE="$2"
        shift 2
        ;;
    -r | --rx_burn_units_file)
        echo "Processing 'rx_burn_units_file' option. Input argument is '$2'"
        RX_BURN_UNITS_FILE="$2"
        shift 2
        ;;
    -t | --toa_file)
        echo "Processing 'toa_file' option. Input argument is '$2'"
        TOA_FILE=$2
        shift 2
        ;;        
    -f | --intensity_file)
        echo "Processing 'intensity_file' option. Input argument is '$2'"
        INTENSITY_FILE=$2
        shift 2
        ;;
    -m | --ms_buildings_tiles_file)
        echo "Processing 'ms_buildings_tiles_file' option. Input argument is '$2'"
        MS_BUILDINGS_TILES_FILE=$2
        shift 2
        ;;
    -c | --critical_habitats_polygon_file)
        echo "Processing 'critical_habitats_polygon_file' option. Input argument is '$2'"
        CRITICAL_HABITATS_POLYGON_FILE=$2
        shift 2
        ;;
    -p | --footprints_polygon_file)
        echo "Processing 'footprints_polygon_file' option. Input argument is '$2'"
        FOOTPRINTS_POLYGONS_FILE=$2
        shift 2
        ;;
    -b | --burned_area_file)
        echo "Processing 'burned_area_file' option. Input argument is '$2'"
        BURNED_AREA_FILE=$2
        shift 2
        ;;
    -d | --building_damage_file)
        echo "Processing 'building_damage_file' option. Input argument is '$2'"
        BUILDING_DAMAGE_FILE=$2
        shift 2
        ;;
    -h | --habitat_damage_file)
        echo "Processing 'habitat_damage_file' option. Input argument is '$2'"
        HABITAT_DAMAGE_FILE=$2
        shift 2
        ;;
    -*|--*)
      echo "Unknown option $1"
      shift
      ;;
  esac
done

SCRIPTDIR=`dirname "$0"`
SCRIPTDIR=`realpath ${SCRIPTDIR}`

CURDIR=`pwd`
CURDIR=`realpath ${CURDIR}`

SCRATCH=${CURDIR}/scratch
INPUTS=${CURDIR}/inputs
OUTPUTS=${CURDIR}/outputs

# inputs
IGNITIONS_CSV=${INPUTS}/ignitions.csv
MS_BUILDINGS_TILES_JSON=${INPUTS}/ms_buildings_tiles.json
RX_BURN_UNITS_DIR=${INPUTS}/rx_burn_units
INTENSITY_DIR=${INPUTS}/intensity
TOA_DIR=${INPUTS}/toa
CRITICAL_HABITATS_POLYGON_DIR=${INPUTS}/critical_habitats_polygon

# temporary files
TMPDIR=${SCRATCH}/tmp
INTERSECTION_OUTPUT=${SCRATCH}/intersection_output.csv

# outputs
DAMAGE_RESPONSE_DIR=${OUTPUTS}/damage_response
BURNED_AREA_DIR=${OUTPUTS}/burned_area
BUILDING_DAMAGE_DIR=${OUTPUTS}/building_damage
HABITAT_DAMAGE_DIR=${OUTPUTS}/habitat_damage
FOOTPRINTS_POLYGONS_JSON=${OUTPUTS}/footprint_polygons.json

rm -f -r $SCRATCH $INPUTS $OUTPUTS

mkdir $SCRATCH $INPUTS $OUTPUTS $TMPDIR
mkdir $RX_BURN_UNITS_DIR $CRITICAL_HABITATS_POLYGON_DIR $INTENSITY_DIR $TOA_DIR

# Copy/Extract inputs
cp $IGNITIONS_FILE $IGNITIONS_CSV
cp $MS_BUILDINGS_TILES_FILE $MS_BUILDINGS_TILES_JSON
tar -xzf $RX_BURN_UNITS_FILE -C $RX_BURN_UNITS_DIR >/dev/null 2>&1
tar -xzf $CRITICAL_HABITATS_POLYGON_FILE -C $CRITICAL_HABITATS_POLYGON_DIR >/dev/null 2>&1
tar -xzf $TOA_FILE -C $TOA_DIR >/dev/null 2>&1
tar -xzf $INTENSITY_FILE -C $INTENSITY_DIR >/dev/null 2>&1


# python ./intersect.py ./data/rxburn/RX_Burn_Units.shp data/yosemite_sample_ignitions.csv data/intersection_output.csv
# python ./create_burned_area_rasters.py ./data/intersection_output.csv ./data/toa/ ./data/burned_area/ ./data/footprints_polygons.json
# python ./create_building_damage_rasters.py ./data/intersection_output.csv ./data/flamelen ./data/MS_bldgs_all_tiles.json ./data/footprints_polygons.json ./data/damage_response ./data/building_damage ./data/tmp
# python ./create_habitat_damage_rasters.py ./data/toa ./data/critical_habitat_polygons/CRITHAB_POLY.shp ./data/footprints_polygons.json ./data/habitat_damage ./data/tmp

# === Execution ===
cd $SCRATCH
# Run intersection
echo "Running Intersection .."
python3 ${SCRIPTDIR}/intersect.py ${RX_BURN_UNITS_DIR}/*.shp ${IGNITIONS_CSV} ${INTERSECTION_OUTPUT}
# Create Burned Area rasters
echo "Creating Burned Area rasters .."
python3 ${SCRIPTDIR}/create_burned_area_rasters.py ${INTERSECTION_OUTPUT} ${TOA_DIR} ${BURNED_AREA_DIR} ${FOOTPRINTS_POLYGONS_JSON}
# Create Building Damage rasters
echo "Creating Building Damage rasters .."
python3 ${SCRIPTDIR}/create_building_damage_rasters.py ${INTERSECTION_OUTPUT} ${INTENSITY_DIR} ${MS_BUILDINGS_TILES_JSON} ${FOOTPRINTS_POLYGONS_JSON} ${DAMAGE_RESPONSE_DIR} ${BUILDING_DAMAGE_DIR} ${TMPDIR}
# Create Habitat Damage rasters
echo "Creating Habitat Damage rasters .."
python3 ${SCRIPTDIR}/create_habitat_damage_rasters.py ${TOA_DIR} ${CRITICAL_HABITATS_POLYGON_DIR}/*.shp ${FOOTPRINTS_POLYGONS_JSON} ${HABITAT_DAMAGE_DIR} ${TMPDIR}

# Compress output directories
DAMAGE_RESPONSE_DIR=${OUTPUTS}/damage_response
BURNED_AREA_DIR=${OUTPUTS}/burned_area
BUILDING_DAMAGE_DIR=${OUTPUTS}/building_damage
HABITAT_DAMAGE_DIR=${OUTPUTS}/habitat_damage
FOOTPRINTS_POLYGONS_JSON=${OUTPUTS}/footprint_polygons.json

cd ${BURNED_AREA_DIR} && tar -czf ${BURNED_AREA_DIR}.tar.gz *.tif
cd ${BUILDING_DAMAGE_DIR} && tar -czf ${BUILDING_DAMAGE_DIR}.tar.gz *.tif
cd ${HABITAT_DAMAGE_DIR} && tar -czf ${HABITAT_DAMAGE_DIR}.tar.gz *.tif

cd $CURDIR

# Clean up and exit:
#rm -f -r $SCRATCH

exit 0