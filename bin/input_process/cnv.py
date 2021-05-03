import pandas as pd
import numpy as np
from itertools import chain
import requests


def read_cnv(inputfile):
    """Function to read CNV input file.
    input: _somatic_cnv.tsv
    output: dataframe"""

    def convert_to_int(row):
        if row['chr'].lower() in ["x", "y"]:
            return row["chr"]
        elif row['chr'] == "null":
            return row["chr"]
        else:
            return int(row["chr"])

    dataframe = pd.read_csv(inputfile, sep="\t")
    dataframe["var_type"] = "cnv"
    dataframe.fillna("null", inplace=True)
    dataframe["chr"] = dataframe["chr"].str.replace("chr", "")
    # convert chromosome number to int.
    dataframe["chr"] = dataframe.apply(convert_to_int, axis=1)
    return dataframe


def chainer(df_series):
    """Function to return flattened list of splitted values. It is used in parse_cnv(), to repeate the rows for splitted values."""
    return list(chain.from_iterable(df_series.str.split(';')))


def split_effect(row):
    """Function to return gene's effect via parsing the effect column of the input dataframe. It is used in parse_cnv() to retrieve the corresponding effect of the gene given in the gene column after ; split."""
    effect = row["effect"].split(";")
    if row["effect"] == "null":
        val = "null"
    elif row["gene"] in row["effect"]:
        for e in effect:
            if row['gene'] in e:
                val = e
            else:
                continue
    else:
        val = "null"
    return val


# room for improvement: make request calls parallel, it takes long time
def get_hgnc_id(dataframe):
    """Function to retrieve HGNC IDs via HGNC REST API"""
    def hgnc_to_str(row):
        if isinstance(row["HGNC_ID"], float):
            return str(int(row["HGNC_ID"]))
        elif isinstance(row["HGNC_ID"], int):
            return str(row["HGNC_ID"])
        else:
            return row["HGNC_ID"]

    if "HGNC_ID" in dataframe.columns:
        dataframe["HGNC_ID"] = dataframe.apply(hgnc_to_str, axis=1)
    else:
        for index, row in dataframe.iterrows():
            if row["gene"] == "null":
                dataframe.at[index, 'HGNC_ID'] = "null"
                continue
            url = "http://rest.genenames.org//search/symbol/{}".format(
                row['gene'])
            response = requests.get(
                url, headers={'Accept': 'application/json'})
            if response.status_code == 200 and response.json()["response"]["numFound"] > 0:
                value = response.json()[
                    "response"]["docs"][0]["hgnc_id"].strip("HGNC: ")
                print(value)
                dataframe.at[index, 'HGNC_ID'] = value
            else:
                dataframe.at[index, 'HGNC_ID'] = "null"
    return dataframe


# room for improvement: check for the column names before processing, it might change.
def parse_cnv(dataframe):
    """Function to process input cnv file. split gene values seperated by ;, repeat rows and reshape the dataframe, get the effect of the splitted gene via parding the effect column, retrieving hgnc ids via hgnc rest api. """
    # get repetations based on split
    lengths = dataframe['gene'].str.split(';').map(len)
    # reshape dataframe
    reshaped_data = {'size': np.repeat(dataframe['size'], lengths), 'type': np.repeat(dataframe['type'], lengths), 'copy_number': np.repeat(
        dataframe['copy_number'], lengths), 'gene': chainer(dataframe['gene']), 'exons': np.repeat(dataframe['exons'], lengths),
        'transcript': np.repeat(dataframe['transcript'], lengths), 'chr': np.repeat(dataframe['chr'], lengths),
        'start': np.repeat(dataframe['start'], lengths), 'end': np.repeat(dataframe['end'], lengths),
        'effect': np.repeat(dataframe['effect'], lengths), 'var_type': np.repeat(dataframe['var_type'], lengths)}
    if "HGNC_ID" in dataframe.columns:
        reshaped_data['HGNC_ID'] = np.repeat(dataframe['HGNC_ID'], lengths)
    reshaped_dataframe = pd.DataFrame(reshaped_data).reset_index(drop=True)
    # seperate effect based on gene column
    reshaped_dataframe["splitted_effect"] = reshaped_dataframe.apply(
        split_effect, axis=1)
    reshaped_dataframe.drop(columns=["effect"], inplace=True)
    reshaped_dataframe.rename(
        columns={"splitted_effect": "effect"}, inplace=True)
    # get hgnc ids
    final_dataframe = get_hgnc_id(reshaped_dataframe)
    return final_dataframe
