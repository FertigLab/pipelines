#!/usr/bin/env python3

#from https://github.com/break-through-cancer/btc-spatial-pipelines/blob/main/.cirro/preprocess.py

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np


def set_params_as_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    samplesheet = pd.DataFrame([ds.params])

    # Save to a file
    samplesheet.to_csv("samplesheet.csv", index=None)

    ds.add_param("input", "samplesheet.csv")

    # Log the samplesheet
    ds.logger.info(samplesheet)


def main():
    ds = PreprocessDataset.from_running()

    set_params_as_samplesheet(ds)

    # log
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()