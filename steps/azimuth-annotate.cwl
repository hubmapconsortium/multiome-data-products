#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: hubmap/multi-maps:latest

inputs:
  muon_processed:
    type: File
    inputBinding:
      position: 0
  tissue:
    type: string?

outputs:
  annotated_mudata:
    type: File
    outputBinding:
      glob: '*_processed.h5ad'
  calculated_metadata_file:
    type: File?
    outputBinding:
      glob: calculated_metadata.json
  cell_type_manifest:
    type: File?
    outputBinding:
      glob: cell_type_manifest.json

baseCommand: ['python3', '/opt/pan_organ_azimuth.py']