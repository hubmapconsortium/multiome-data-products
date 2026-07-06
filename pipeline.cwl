#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: Pipeline for concatenating sc-ATAC-seq datasets into MuData object

requirements:
  ScatterFeatureRequirement: {}

inputs: 
    data_directory:
        label: "Path to directory containing cell by gene and cell by bin files"
        type: Directory
    
    uuids_file:
        label: "Path to a file containing a list of uuids and other metadata for the dataset to be indexed"
        type: File
    
    tissue:
        label: "Two letter tissue type code"
        type: string?

outputs:
    mudata_raw:
        type: File
        outputSource: concatenate/mudata_raw
    
    muon_processed: 
        type: File
        outputSource: downstream/muon_processed
    
    joint_embedding:
        type: File
        outputSource: downstream/joint_embedding
    
    final_metadata:
        type: File
        outputSource: downstream/metadata_with_cell_types

steps:

    - id: concatenate
      in: 
        - id: data_directory
          source: data_directory
        - id: uuids_file
          source: uuids_file
        - id: tissue
          source: tissue
    
      out:
        - mudata_raw
        - metadata_json
      run: steps/concatenate.cwl
      label: "Concatenates h5ad files in directory"

    - id: downstream
      in:
        - id: mudata_raw
          source: concatenate/mudata_raw
        - id: tissue
          source: tissue
        - id: metadata_json
          source: concatenate/metadata_json
      
      out:
        - muon_processed
        - mofa_out
        - joint_embedding
        - rna_embedding
        - atac_embedding
        - final_metadata_json
      run: steps/downstream.cwl


    - id: pan_organ_azimuth
      in:
        - id: muon_processed
          source: downstream/muon_processed
        - id: tissue:
          source: tissue
        - id: metadata_json
          source: downstream/final_metadata_json
      out:
        - annotated_mudata
        - metadata_with_cell_types
      run: steps/azimuth-annotate.cwl
      label: "Adds azimuth annotations"
      