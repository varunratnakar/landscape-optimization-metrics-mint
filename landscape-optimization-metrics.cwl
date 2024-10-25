arguments:
- --
baseCommand: /landoptmet/run.sh
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerImageId: kcapd/landscape-optimization-metrics-mint:latest
inputs:
  ignitions_file:
    inputBinding:
      prefix: --ignitions_file
    type: File
  rx_burn_units_file:
    inputBinding:
      prefix: --rx_burn_units_file
    type: File
  intensity_file:
    inputBinding:
      prefix: --intensity_file
    type: File
  toa_file:
    inputBinding:
      prefix: --toa_file
    type: File
  ms_buildings_tiles_file:
    inputBinding:
      prefix: --ms_buildings_tiles_file
    type: File
  critical_habitats_polygon_file:
    inputBinding:
      prefix: --critical_habitats_polygon_file
    type: File
outputs:
  burned_area_file:
    outputBinding:
      glob: ./outputs/burned_area.tar.gz
    type: File
  building_damage_file:
    outputBinding:
      glob: ./outputs/building_damage.tar.gz
    type: File
  habitat_damage_file:
    outputBinding:
      glob: ./outputs/habitat_damage.tar.gz
    type: File
  footprints_polygon_file:
    outputBinding:
      glob: ./outputs/footprint_polygons.json
    type: File
requirements:
  NetworkAccess:
    networkAccess: true
