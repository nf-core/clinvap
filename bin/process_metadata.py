#!/usr/bin/env python3

import json
import sys
import pandas as pd
import copy
import os
import math
from operator import itemgetter
import re


# PROCESS METADATA

def get_field(key, metadata):
    """Process metadata function."""
    val = metadata.get(key, "null")
    if not val:
        val = "null"
    return val


def give_na_diagnosis_tag(report):
    """Process metadata function."""
    copy_report = copy.deepcopy(report)
    copy_report["seen_in_diagnosis"] = "Not applicable"
    return copy_report


def in_vocabulary(vocab_str, column_name, result_list):
    """Process metadata function."""
    vocab_list = vocab_str.split(",")
    for v in vocab_list:
        if column_name == "icd10":
            if not vocab[vocab[column_name].str.contains(v, na=False)].empty:
                result_list.append(v)
    else:
        if not vocab[vocab[column_name] == v].empty:
            result_list.append(v)
    return result_list


def get_tumor_name(vocab_list, vocab_df, vocab_df_column_name):
    """Process metadata function.
    Function to obtain tumor type keywords as given databases from disease ontology id, icd10 codes or do_"""
    conv_list = []
    for v in vocab_list:
        cond = vocab_df[vocab_df_column_name].str.contains(v, na=False)
        indices = cond[cond].index
        conv_list = vocab_df.loc[indices, "tumor_type"].tolist()
    return conv_list


def tumor_type_match(section, keywords, match_type):
    """Process metadata function."""
    for d in section:
        if not keywords:
            d[match_type] = "null"
        else:
            d[match_type] = "False"
        for k in keywords:
            if k in d["db_tumor_repr"].lower():
                d[match_type] = "True"
                break
            else:
                continue
    return section


def get_matching_annotations(annotation_type):
    """Process metadata function."""
    match_list = []
    try:
        match = icd10_annotation_dict[icd10][annotation_type]
        if isinstance(match, str):
            cond = vocab[annotation_type] == match
            match_list = vocab[cond]["tumor_type"].tolist()
        return match_list
    except KeyError:
        print("Not a cancer ICD10 code. Exiting...")
        exit(0)


# TODO bunu metinde anlatmak lazim

def calculate_diagnosis_match_score(section):
    """Process metadata function."""
    for ann in section:
        if ann.get("seen_in_diagnosis") == "True":
            ann["score"] = 1.0
        else:
            diagnosis_annotation = [ann.get("system", 0), ann.get("organ", 0), ann.get(
                "histology_0", 0), ann.get("histology_1", 0), ann.get("histology_2", 0)]
            if set(diagnosis_annotation) == {0}:
                if ann.get("seen_in_diagnosis") == "False":
                    ann["score"] = 0.0
                elif ann.get("seen_in_diagnosis", "na") == "na":
                    ann["score"] = -1.0
            else:
                num_trues = diagnosis_annotation.count("True")
                num_falses = diagnosis_annotation.count("False")
                ann["score"] = round(num_trues / (num_trues + num_falses), 2)
    return section


def clean_up(section):
    """Process metadata function."""
    for ann in section:
        try:
            del ann["system"]
            del ann["organ"]
            del ann["histology_0"]
            del ann["histology_1"]
            del ann["histology_2"]
        except KeyError:
            pass
    return section


def filter_on_diagnosis(section):
    """Filter, sort, prioritize function"""
    filtered_section = [
        d for d in section if d["score"] == 1]
    return filtered_section


def sort_on_diagnosis(section):
    """Filter, sort, prioritize function"""
    sorted_list = sorted(section, key=itemgetter('score'), reverse=True)
    return sorted_list


def priori_on_diagnosis(section):
    """Filter, sort, prioritize function"""
    new_section = []
    for ann in section:
        evidence_level = re.sub("-.*", "", ann["match_level"])
        if evidence_level in ["D", "E"]:
            if ann["score"] == 1:
                new_section.append(ann)
            else:
                continue
        else:
            new_section.append(ann)
    sorted_new_section = sort_on_diagnosis(new_section)
    return sorted_new_section


main_json = sys.argv[1]
metadata_json = sys.argv[2]
output_json = sys.argv[3]
db_lookup = sys.argv[4]
icd10_lookup = sys.argv[5]
process_type = sys.argv[6]


if len(sys.argv) != 7:
    print("usage: {} <main_json> <metadata_json> <output_json> <disease_ontology_map> <icd10_lookup_table> <process_type>".format(
        sys.argv[0]))
    exit(1)


# check if process_type is correct
OPTIONS = ["sort", "filter", "prioritize,filter,sort"]
if not process_type in OPTIONS:
    raise ValueError(
        "Invalid process type input. Choose among 'sort', 'filter', 'prioritize, filter, sort'")


# Read metadata file
with open(metadata_json) as m:
    metadata = json.load(m)

# Read main report file
with open(main_json) as r:
    main_report = json.load(r)

# merge main report content with metadata
report = metadata.copy()
report.update(main_report)

# Load disease vocabulary TODO degisecek
dtypes = {"tumor_type": "str", "disease_name": "str", "doid": "str", "icd10": "str", "do_name": "str",
          "abbr": "str", "system": "str", "organ": "str", "histology_0": "str", "histology_1": "str", "histology_2": "str"}
columns = ["tumor_type", "doid", "icd10", "do_name", "system",
           "organ", "histology_0", "histology_1", "histology_2"]
vocab = pd.read_csv(db_lookup, sep="\t", dtype=dtypes)[columns]
vocab = vocab.apply(lambda x: x.str.lower())

# Load ICD10 lookup table
dtypes = {"code": "str", "long_desc": "str", "system": "str", "organ": "str",
          "histology_0": "str", "histology_1": "str", "histology_2": "str"}
columns = ["code", "system", "organ",
           "histology_0", "histology_1", "histology_2"]
icd10_vocab = pd.read_csv(icd10_lookup, sep="\t", dtype=dtypes)[columns]
icd10_vocab = icd10_vocab.apply(lambda x: x.str.lower())


do_name = get_field("do_name", report).lower()
doid = get_field("doid", report).lower()
icd10 = get_field("icd10", report).lower()


if all(val == "null" for val in [do_name, doid, icd10]):
    tagged_report = give_na_diagnosis_tag(report)
    with open(output_json, "w") as o:
        json.dump(tagged_report, o, indent=4)
    print("Empty metadata file.")
    sys.exit(0)


KEYS = ["driver_table", "direct_pharm_table", "pharm_table"]

vocab_doname = []
vocab_doid = []
vocab_icd10 = []


if do_name != "null":
    vocab_doname = in_vocabulary(do_name, "do_name", vocab_doname)
    vocab_doname_tn = get_tumor_name(vocab_doname, vocab, "do_name")
else:
    vocab_doname_tn = []


if doid != "null":
    vocab_doid = in_vocabulary(doid, "doid", vocab_doid)
    vocab_doid_tn = get_tumor_name(vocab_doid, vocab, "doid")
else:
    vocab_doid_tn = []


if icd10 != "null":
    vocab_icd10 = in_vocabulary(icd10, "icd10", vocab_icd10)
    vocab_icd10_tn = get_tumor_name(vocab_icd10, vocab, "icd10")
else:
    vocab_icd10_tn = []


keyword_lists = [vocab_doname_tn, vocab_doid_tn, vocab_icd10_tn]
if not all(result == [] for result in keyword_lists):
    keywords = []
    for x in keyword_lists:
        keywords.extend(x)
    keywords = list(set(keywords))
    report["driver_table"] = tumor_type_match(
        report["driver_table"], keywords, "seen_in_diagnosis")
    report["direct_pharm_table"] = tumor_type_match(
        report["direct_pharm_table"], keywords, "seen_in_diagnosis")
    report["pharm_table"] = tumor_type_match(
        report["pharm_table"], keywords, "seen_in_diagnosis")


if icd10:
    cond = icd10_vocab["code"] == icd10
    icd10_annotation_dict = icd10_vocab[cond].set_index(
        "code").to_dict("index")

    same_system_tumors_list = get_matching_annotations("system")
    same_organ_tumors_list = get_matching_annotations("organ")
    same_histology_0_tumors_list = get_matching_annotations("histology_0")
    same_histology_1_tumors_list = get_matching_annotations("histology_1")
    same_histology_2_tumors_list = get_matching_annotations("histology_2")

    categories = ["system", "organ", "histology_0",
                  "histology_1", "histology_2"]
    for word in KEYS:
        report[word] = tumor_type_match(
            report[word], same_system_tumors_list, "system")
        report[word] = tumor_type_match(
            report[word], same_organ_tumors_list, "organ")
        report[word] = tumor_type_match(
            report[word], same_histology_0_tumors_list, "histology_0")
        report[word] = tumor_type_match(
            report[word], same_histology_1_tumors_list, "histology_1")
        report[word] = tumor_type_match(
            report[word], same_histology_2_tumors_list, "histology_2")


for word in KEYS:
    report[word] = calculate_diagnosis_match_score(report[word])
    report[word] = clean_up(report[word])


# APPLY FILTERS


# Check if "seen in diagnosis": "not applicable" is in the report
# if so, do nothing, give report back
if report.get("seen_in_diagnosis") == "Not applicable":
    with open(output_json, "w") as o:
        json.dump(report, o, indent=4)
    os.remove(main_json)
    sys.exit(0)


if process_type == "sort":
    report["driver_table"] = sort_on_diagnosis(report["driver_table"])
    report["direct_pharm_table"] = sort_on_diagnosis(
        report["direct_pharm_table"])
    report["pharm_table"] = sort_on_diagnosis(report["pharm_table"])


elif process_type == "filter":
    report["driver_table"] = filter_on_diagnosis(report["driver_table"])
    report["direct_pharm_table"] = filter_on_diagnosis(
        report["direct_pharm_table"])
    report["pharm_table"] = filter_on_diagnosis(report["pharm_table"])


elif process_type == "prioritize,filter,sort":
    report["driver_table"] = sort_on_diagnosis(report["driver_table"])
    report["direct_pharm_table"] = priori_on_diagnosis(
        report["direct_pharm_table"])
    report["pharm_table"] = priori_on_diagnosis(report["pharm_table"])


else:
    print("Procees type option error.")
    sys.exit(1)


with open(output_json, "w") as o:
    json.dump(report, o, indent=4)
