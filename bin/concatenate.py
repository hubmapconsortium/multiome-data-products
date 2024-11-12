#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import muon as mu
import numpy as np
import os
import pandas as pd
import requests
import uuid
import yaml


def get_tissue_type(dataset: str) -> str:
    organ_dict = yaml.load(open("/opt/organ_types.yaml"), Loader=yaml.BaseLoader)
    organ_code = requests.get(
        f"https://entity.api.hubmapconsortium.org/dataset/{dataset}/organs/"
    )
    organ_name = organ_dict[organ_code]
    return organ_name.replace(" (Left)", "").replace(" (Right)", "")


def convert_tissue_code(tissue_code):
    with open("/opt/organ_types.yaml", 'r') as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    tissue_name = data.get(tissue_code)['description']
    return tissue_name


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath
                

def find_file_pairs(directory):
    raw_mdata_patterns = ["mudata_raw.h5mu"]
    processed_mdata_patterns = ["secondary_analysis.h5mu"]
    raw_mdata_file = find_files(directory, raw_mdata_patterns)
    processed_mdata_file = find_files(directory, processed_mdata_patterns)
    return (raw_mdata_file, processed_mdata_file)


def annotate_h5ads(
    adata_file, tissue_type: str, uuids_df: pd.DataFrame):
    # Get the directory
    data_set_dir = fspath(adata_file.parent.stem)
    # And the tissue type
    tissue_type = tissue_type if tissue_type else get_tissue_type(data_set_dir)
    dense_adata = anndata.read_h5ad(adata_file)
    adata = make_new_anndata_object(dense_adata)
    del dense_adata
    adata_copy = adata.copy()
    adata_copy.obs["barcode"] = adata.obs.index
    adata_copy.obs["barcode"] = adata_copy.obs["barcode"].str.replace("BAM_data#", "", regex=False)
    adata_copy.obs["dataset"] = data_set_dir
    
    cell_ids_list = [
        "-".join([data_set_dir, barcode]) for barcode in adata_copy.obs["barcode"]
    ]
    adata_copy.obs["cell_id"] = pd.Series(
        cell_ids_list, index=adata_copy.obs.index, dtype=str
    )
    adata_copy.obs.set_index("cell_id", drop=True, inplace=True)
    return adata_copy


def create_json(tissue, data_product_uuid, creation_time, uuids, hbmids, cell_count, file_size):
    bucket_url = f"https://hubmap-data-products.s3.amazonaws.com/{data_product_uuid}/"
    metadata = {
        "Data Product UUID": data_product_uuid,
        "Tissue": convert_tissue_code(tissue),
        "Assay": "10X Multiome",
        "URL": bucket_url + f"{tissue}.h5mu",
        "Creation Time": creation_time,
        "Dataset UUIDs": uuids,
        "Dataset HBMIDs": hbmids,
        "Total Cell Count": cell_count,
        "Raw File Size": file_size
    }
    print("Writing metadata json")
    with open(f"{data_product_uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)


def make_mudata(raw_mdata, processed_mdata):
    mdata = mu.MuData({"atac_raw_mdata": raw_mdata, "atac_processed_mdata": processed_mdata})
    mu.pp.intersect_obs(mdata)
    return mdata


def annotate_mudata(mdata, uuids_df):
    merged = uuids_df.merge(mdata.obs, left_on="uuid", right_on="dataset", how="inner")
    merged = merged.set_index(mdata.obs.index)
    merged = merged.drop(columns=["Unnamed: 0"])
    merged = merged.fillna(np.nan)
    merged["age"] = pd.to_numeric(merged["age"])
    return merged


def make_new_anndata_object(adata):
    new_adata = anndata.AnnData(X=adata.X, obs=pd.DataFrame(index=adata.obs.index), var=adata.var)
    return new_adata


def load_mudata(mdata_file):
    mdata = mu.read_h5mu(mdata_file)
    return mdata


def main(data_directory: Path, uuids_file: Path, tissue: str = None):
    output_file_name = f"{tissue}" if tissue else "atac"
    uuids_df = pd.read_csv(uuids_file, sep="\t", dtype=str)
    uuids_list = uuids_df["uuid"].to_list()
    hbmids_list = uuids_df["hubmap_id"].to_list()
    directories = [data_directory / Path(uuid) for uuid in uuids_df["uuid"]]
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories if len(listdir(directory))>1]
    print("Annotating objects")
    raw_mdatas = [
        load_mudata(file_pair[0])
        for file_pair in file_pairs
    ]

    # processed_mdatas = [     
    #     load_mudata(file_pair[1])
    #     for file_pair in file_pairs
    # ]

    print("Concatenating objects")
    raw_mdata_concat = mu.concat(raw_mdatas)

    creation_time = str(datetime.now())
    data_product_uuid = str(uuid.uuid4())
    total_cell_count = raw_mdatas.obs.shape[0]
    raw_mdata_concat.obs = annotate_mudata(raw_mdata_concat, uuids_df)
    raw_mdata_concat.uns["creation_data_time"] = creation_time
    raw_mdata_concat.uns["datasets"] = hbmids_list
    raw_mdata_concat.uns["uuid"] = data_product_uuid
    raw_mdata_concat.write(f"{output_file_name}.h5mu")
    file_size = os.path.getsize(f"{output_file_name}.h5mu")
    create_json(tissue, data_product_uuid, creation_time, uuids_list, hbmids_list, total_cell_count, file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("data_directory", type=Path)
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory, args.uuids_file, args.tissue)