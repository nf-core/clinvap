#!/usr/bin/env python3

import db_process.process_pharmacogenomics_effects as pharma
import db_process.process_mechanistic_drugs as mechanistic
import table_content.create_content as content
import table_content.prioritize as prior
import table_content.handling as handle
import db_process.process_driver as driver
from collections import defaultdict
import db_process.helpers as helper
import db_process.query as query
import input_process.snv as snv
import input_process.cnv as cnv
import numpy as np
from itertools import chain
import functools
from cyvcf2 import VCF
import pandas as pd
import argparse
import requests
import json
import sys
import time


# Process arguments
parser = argparse.ArgumentParser(
    prog="reporting", description='Reporting script arguments')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument(
    "-i", "--input", required=True, help="Path to input VCF")
requiredNamed.add_argument(
    "-o", "--output", required=True, help="Output report name")
requiredNamed.add_argument("-g", "--genome_assembly",
                           required=True, help="Genome assembly version")
requiredNamed.add_argument(
    "-k", "--knowledgebase", required=True, help="Path to knowledgebase JSON file")
parser.add_argument("-c", "--cnv", required=False,
                    help="Path to somatic CNV file")

args = parser.parse_args()
input_vcf = args.input
output_json = args.output
genome = args.genome_assembly
KNOWLEDGEBASE = query.read_knowledgebase(args.knowledgebase)
input_cnv = args.cnv


###############################################################################################################

# PROCESS VCF TO GET INFO THAT CNV ANALYSIS DEPEND ON

# Read in and process VCF file, parse vcf function is checking whether mvld is empty or not.
mvld = snv.parse_vcf(input_vcf).rename(columns={"SYMBOL": "gene"})


# Filter and divide dataframe
mvld_high_moderate = snv.get_high_moderate_effect(mvld)
mvld_modifier = snv.get_modifier_effect(mvld)

# Check if mvld_high_moderate dataframe is empty. If it is empty, use mvld_moderate
if not mvld_high_moderate.empty:
    main_mvld = mvld_high_moderate
elif not mvld_modifier.empty:
    handle.empty_high_impact_mvld()
    main_mvld = mvld_modifier

# HGVSp one letter representation
main_mvld = snv.one_letter_repr(main_mvld)


# list of genes
snv_gene_list = query.get_gene_list(main_mvld)[0]
snv_gene_dict = query.get_gene_list(main_mvld)[1]


# fields needed from mlvd:
from_mlvd = main_mvld[[
    "HGNC_ID", "gene", "one_letter_repr", "Consequence", "vaf", "var_type"]].rename(columns={"SYMBOL": "gene"})


###############################################################################################################

########################
### Read Process CNV ###
########################

# PROCESS CNV

if input_cnv:
    df_cnv = cnv.read_cnv(input_cnv)
    cnv_flag = True
else:
    # if cnv is empty, this is placeholder for filter_pair_variant function.
    cnv_gene_tuples = [()]
    cnv_flag = False


if cnv_flag and not df_cnv.empty:
    df_cnv_processed = cnv.parse_cnv(df_cnv)

    #  Get driver gene annotation for CNVs
    cnv_driver_ann = query.get_driver_annotation(
        df_cnv_processed, KNOWLEDGEBASE)[2]

    # Query database - Get mechanistic drug targets, cancer drugs
    cnv_mechanistic_drug_targets = query.get_mechanistic_drugs(
        df_cnv_processed, KNOWLEDGEBASE)

    # Query database - Get direct drug targets
    cnv_direct_variant_targets = query.ann_pharm_variant(
        df_cnv_processed, genome, KNOWLEDGEBASE, "CNA")

    # Filter direct drug targets for variant_class
    cnv_direct_variant_targets = query.filter_main_variant_class(
        cnv_direct_variant_targets, "CNA")

    # Query database - Get drugs targeting affected gene
    cnv_gene_targets = query.ann_pharm_gene(
        df_cnv_processed, KNOWLEDGEBASE, genome)

    # Filter drug targets for variant_class
    cnv_gene_targets = query.filter_main_variant_class(cnv_gene_targets, "CNA")

    # Get gene, type tuples reqired in pair variant filtering in process_pharmacogenomics part.
    cnv_gene_tuples = query.get_gene_list(df_cnv_processed)[2]

    # fields needed from cnv dataframe
    from_cnv = df_cnv_processed[["HGNC_ID",
                                 "gene", "type", "copy_number", "effect"]]

elif cnv_flag and df_cnv.empty:
    print("CNV file is empty. There will be no results related with copy number variants.")
    sys.exit(0)
    # cnv_flag = False
    # # we need empty cnv_gene_tuple list as a place holder for upcaming parts (filter_gene_pair function)
    # cnv_gene_tuple = [()]

if (cnv_flag) and (not cnv_mechanistic_drug_targets) and (not cnv_driver_ann) and (not cnv_gene_targets) and (not cnv_direct_variant_targets):
    print("Database queries did not return any results for CNVs. There will be no results related with copy number variants.")
    sys.exit(0)
    # cnv_flag = False

# Check if mechanistic drug targets are empty or not.
if cnv_mechanistic_drug_targets:
    cnv_mechanistic_flag = True
else:
    cnv_mechanistic_flag = False
    handle.empty_mechanistic("CNV")


# Check if driver gene annotation dictionary is empty or not.
if cnv_driver_ann:
    cnv_driver_flag = True
else:
    cnv_driver_flag = False


###############################################################################################################

############
### CNVs ###
############


# CNV DRIVER ANNOTATION


# add tumor string to driver dictionary and remove tumor_type list of dictionaries
for key in cnv_driver_ann:
    helper.add_tumor_string(key, cnv_driver_ann)

if cnv_driver_flag:
    df_cnv_driver_ann = driver.get_driver_df(cnv_driver_ann)
else:
    df_cnv_driver_ann = handle.empty_driver_annotation("CNVs")

# Start building up driver content for CNVs. References will be mapped later.

cnv_driver_ann = df_cnv_driver_ann[[
    "hgnc_id", "driver_role", "reference_id", "tumor_list", "db_tumor_repr", "var_type"]]


from_cnv = df_cnv_processed[["type", "copy_number",
                             "gene", "var_type", "HGNC_ID", "effect"]]

cnv_driver_content = pd.merge(from_cnv, cnv_driver_ann, how="right",
                              left_on=["HGNC_ID", "var_type"], right_on=["hgnc_id", "var_type"]).drop(columns=["HGNC_ID"])


# Get references
cnv_ref_driver = content.get_references(df_cnv_driver_ann)


##############################################################################################################

# CNV DIRECT VARIANT TARGETS

if cnv_direct_variant_targets:
    cnv_empty_variant_targets = False

    # DIRECT PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from variant targets
    cnv_direct_pharmacogenomics = pharma.get_subset_dict(
        cnv_direct_variant_targets, "pharmacogenomics_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        cnv_direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    # Direct pharmacogenomics effect dictionary to dataframe
    df_cnv_direct_pharmacogenomics = pharma.dict_to_dataframe(
        cnv_direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    if df_cnv_direct_pharmacogenomics.empty:
        df_cnv_direct_pharmacogenomics = handle.empty_dataframe_direct_pharmacogenomics(
            "CNVs")

    # Direct target combination will not return anything from the database, such a case for CNV does not exists. All the combination
    # cases for CNV as main variant has null in cordinate related fields.
    # Skipping CNV DIRECT PHARMACOGENOMICS FOR COMBINED VARIANTS. To be implemented later after a database update, if we receive such cases

    # Direct adverse combination effects will not return anything from the database, such a case for CNV does not exits.
    # Skipping CNV DIRECT ADVERSE EFFECTS COMBINED VARIANTS. To be implemented after a database update, if we receive such cases.

    # Direct adverse effects will not return anything from the database, such a case for CNV does not exits.
    # The adverse effects entries have only MUT, null and BIA as variant class.
    # Skipping CNV DIRECT ADVERSE EFFECTS. To be implemented after a database update, if we receive such cases.


else:
    df_cnv_direct_pharmacogenomics, df_cnv_direct_adverse = handle.empty_direct_variant_annotation(
        "CNV")
    cnv_empty_variant_targets = True

#############################################################################################################

# CNV GENE TARGETS

if cnv_gene_targets:
    cnv_empty_gene_targets = False

    # PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from gene targets
    cnv_pharmacogenomics_therapeutics = pharma.get_subset_dict(
        cnv_gene_targets, "pharmacogenomics_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        cnv_pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    # Pharmacogenomics effects to dataframe
    df_cnv_pharmacogenomics = pharma.dict_to_dataframe(
        cnv_pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    if df_cnv_pharmacogenomics.empty:
        df_cnv_pharmacogenomics = handle.empty_dataframe_pharmacogenomics(
            "CNVs")

    # PHARMACOGENOMICS COMBINED VARIANTS

    # seperate non-empty pharmacogenomics combined variant therapeutics from gene targets
    cnv_pharmacogenomics_combined_variants_therapeutics = pharma.get_subset_dict(
        cnv_gene_targets, "pharmacogenomics_combined_variants_therapeutics")

    # Filter based on the distrupted pair gene
    cnv_pharmacogenomics_combined_variants_therapeutics = pharma.filter_gene_pair(
        snv_gene_list, cnv_gene_tuples, cnv_pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(cnv_pharmacogenomics_combined_variants_therapeutics,
                                  "pharmacogenomics_combined_variants_therapeutics")

    df_cnv_pharmacogenomics_combined = pharma.dict_to_dataframe(
        cnv_pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")

    if df_cnv_pharmacogenomics_combined.empty:
        df_cnv_pharmacogenomics_combined = handle.empty_dataframe_pharmacogenomics_combined(
            "CNVs")

    # Direct adverse effects will not return anything from the database, such a case for CNV does not exits.
    # The adverse effects entries have only MUT, null and BIA as variant class.
    # Skipping CNV  ADVERSE EFFECTS. To be implemented after a database update, if we receive such cases.

    df_cnv_adverse = handle.empty_dataframe_adverse()

    # Direct adverse combination effects will not return anything from the database, such a case for CNV does not exits.
    # Skipping CNV ADVERSE EFFECTS COMBINED VARIANTS. To be implemented after a database update, if we receive such cases.

    # ADVERSE EFFECTS COMBINED VARIANTS
    df_cnv_adverse_combined = handle.empty_dataframe_adverse_combined()

else:
    df_cnv_pharmacogenomics, df_cnv_pharmacogenomics_combined, df_cnv_adverse, df_cnv_adverse_combined = handle.empty_variant_annotation(
        "CNV")
    cnv_empty_gene_targets = True


if cnv_empty_gene_targets and cnv_empty_variant_targets:
    cnv_skip_pharmacodynamics_content = True
else:
    cnv_skip_pharmacodynamics_content = False


#############################################################################################################

if not cnv_skip_pharmacodynamics_content:

    # CREATE CNV PHARMACOGENOMICS TABLES CONTENT

    COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source",
               "variant_drug_association", "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "info"]
    COMB_COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source", "variant_drug_association",
                    "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "info", "info_pair"]

    # Content from direct pharmacogenomics effect
    cnv_direct_pharm_content = content.cnv_get_content(
        df_cnv_direct_pharmacogenomics, from_cnv, COLUMNS, direct=1, pharm=1, notcombined=1)

    # Remove direct results covered by df_pharmacogenomics
    cnv_pharm = pd.concat(
        [df_cnv_pharmacogenomics, df_cnv_direct_pharmacogenomics]).drop_duplicates(keep=False)

    # Content from pharmacogenomics effect
    cnv_pharm_content = content.cnv_get_content(
        cnv_pharm, from_cnv, COLUMNS, pharm=1, indirect=1, notcombined=1)

    # Content from pharmacogenomics variant combination effect - Direct combination does not possible for CNV
    cnv_pharm_comb_content = content.cnv_get_content(
        df_cnv_pharmacogenomics_combined, from_cnv, COMB_COLUMNS, mvld=from_mlvd, pharm=1, indirect=1, combined=1)

    # move same variant for both genes, if there are any, content to direct pharm table.
    cond = (cnv_pharm_comb_content["info"] == "Same gene, same variant, same consequence.") & (
        cnv_pharm_comb_content["info_pair"] == "Same gene, same variant, same consequence.")
    rows = cnv_pharm_comb_content.loc[cond, :]
    if not rows.empty:
        cnv_direct_pharm_content = cnv_direct_pharm_content.append(
            rows, ignore_index=True)
        cnv_pharm_comb_content.drop(rows.index, inplace=True)

    # Add match level to see whether the results are same gene, same, variant,same consequence or whether the
    # variant and consequence differs in the input than the database result. For detais see documentation

    cnv_sorted_direct_pharm_content = prior.compatibility_sort(
        cnv_direct_pharm_content)

    # Join references of the same rows (collapsing same rows with different reference ids)
    cols = ["gene", "drug_name", "variant_drug_association",
            "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "match_level"]

    cnv_sorted_direct_pharm_content = prior.unnest_ref(
        cnv_sorted_direct_pharm_content, cols)

    cnv_pharm_content = pd.concat([cnv_pharm_content, cnv_pharm_comb_content])
    cnv_sorted_pharm_content = prior.compatibility_sort(cnv_pharm_content)
    cnv_sorted_pharm_content = prior.unnest_ref(cnv_sorted_pharm_content, cols)

    # Get references
    cnv_ref_direct_pharm = content.get_references(cnv_direct_pharm_content)

    cnv_ref_pharm = content.get_references(cnv_pharm_content)

else:
    cnv_sorted_direct_pharm_content = handle.empty_pharmacodynamics_place_holder()
    cnv_sorted_pharm_content = handle.empty_pharmacodynamics_place_holder()
    cnv_ref_pharm = handle.empty_reference_place_holder()
    cnv_ref_direct_pharm = handle.empty_reference_place_holder()

#############################################################################################################

if cnv_mechanistic_flag:

    # Update mechanistic drug targets dict
    cnv_mechanistic_drug_targets = mechanistic.reshape(
        cnv_mechanistic_drug_targets, {})

    # collapse info from sources into one representation
    temp = {}
    for key, annotation in cnv_mechanistic_drug_targets.items():
        temp[key] = []
        for db_id, info in annotation.items():
            if len(info) == 1:
                temp[key].append(info[0])
            elif db_id != "null":
                approval = set()
                drug = set()
                source = []
                source_pmid = []
                reference = []
                reference_source = set()
                for drug_info in info:
                    helper.populate_set(drug_info, "approval_status", approval)
                    helper.populate_set(drug_info, "drug_name", drug)
                    helper.populate_list(drug_info, "source_name", source)
                    helper.populate_list(drug_info, "source_pmid", source_pmid)
                    helper.populate_list(drug_info, "reference_id", reference)
                    helper.populate_set(
                        drug_info, "reference_source", reference_source)

                approval_status = mechanistic.assign_approval_status(approval)
                name = list(drug)[0]
                agg_source = "|".join(source)
                agg_pmid = "|".join(source_pmid)
                agg_ref = "|".join(reference)
                agg_ref_source = "|".join(list(reference_source))

                temp[key].append(
                    {"drugbank_id": db_id, "drug_name": name, "approval_status": approval_status,
                        "source_name": agg_source, "reference_id": agg_ref, "reference_source": agg_ref_source,
                        "source_pmid": agg_pmid})
            else:
                for drug_info in info:
                    temp[key].append(drug_info)
    cnv_mechanistic_drug_targets = temp

    # filter not-approved drugs from mechanistic drugs
    cnv_mechanistic_drugs_filtered = mechanistic.filter_approval_info(
        cnv_mechanistic_drug_targets)

    # seperate investigational and approved mechanistic drugs
    cnv_approved_mechanistic_drugs = mechanistic.divide_approved_investigational(
        cnv_mechanistic_drugs_filtered)[0]

    cnv_investigational_mechanistic_drugs = mechanistic.divide_approved_investigational(
        cnv_mechanistic_drugs_filtered)[1]

    # Mechanistic drug targets information to dataframe
    if cnv_approved_mechanistic_drugs:
        cnv_main_mechanistic = cnv_approved_mechanistic_drugs
    else:
        handle.empty_approved_mechanistic()
        cnv_main_mechanistic = cnv_investigational_mechanistic_drugs

    df_cnv_mechanistic_targets = mechanistic.dict_to_dataframe(
        cnv_main_mechanistic)
else:
    df_cnv_mechanistic_targets = pd.DataFrame(
        columns=['drugbank_id', 'drug_name', 'approval_status', 'source_name', 'reference_id', 'reference_source', 'source_pmid', 'hgnc_id'])


# CREATE MECHANISTIC DRUG TARGETS TABLE CONTENT

if not cnv_skip_pharmacodynamics_content:
    # remove the content found in pharmacogenomics therapeutics - bunlardan biri bossa bura patlar mi?
    cnv_all_pharm_drugs = pd.concat([cnv_direct_pharm_content[["drug_name", "hgnc_id"]], cnv_pharm_content[[
        "drug_name", "hgnc_id"]], cnv_pharm_comb_content[["drug_name", "hgnc_id"]]]).drop_duplicates()

    cond = df_cnv_mechanistic_targets["drug_name"].isin(
        cnv_all_pharm_drugs["drug_name"]) & df_cnv_mechanistic_targets["hgnc_id"].isin(cnv_all_pharm_drugs["hgnc_id"])
    cnv_mechanistic = df_cnv_mechanistic_targets.loc[~cond]
else:
    cnv_mechanistic = df_cnv_mechanistic_targets

MECH_COLS = ["gene", "drug_name", "approval_status",
             "reference_id", "reference_source", "hgnc_id"]
cnv_mechanistic_content = content.get_mechanistic_content(
    cnv_mechanistic, from_cnv, MECH_COLS)


cnv_ref_mechanistic = content.get_references(cnv_mechanistic_content)


##############################################################################################

# CREATE APPENDIX TABLES

# appendix table for all variants
cnv_appendix_content = content.appendix_variants(df_cnv_processed, "CNV")


# Appendix table for references


# concatanate reference dataframes
dataframes = [cnv_ref_driver, cnv_ref_direct_pharm,
              cnv_ref_pharm, cnv_ref_mechanistic]

cnv_ref_df = functools.reduce(
    lambda left, right: pd.concat([left, right], ignore_index=True), dataframes)  # concatanate dataframes to create one reference dataframe


# create pubmed query string for eutils
pubmed_ids, ref_map_df = content.id_string(cnv_ref_df)
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?retmode=json;db=pubmed;id={}".format(
    pubmed_ids)


# get API results, create dataframe with reference information
ref_info_df = content.get_ref_details(URL)

# merge reference map and reference info dataframes two get a complete reference map.
ref_map_complete = ref_map_df.merge(
    ref_info_df, how="left", left_on="reference_id", right_index=True)

cnv_reference_content, reference_map_dict = content.reference_table(
    ref_map_complete)

# Map reference ids to the dataframes and sort based on the number of the references they have

# get reference map ids for driver_content and sort based n number of references
# listing unknown type of drivers under known typed ones.
cnv_driver_content = content.reference_map_ids(
    cnv_driver_content, reference_map_dict)
cnv_driver_content = prior.driver_type_sort(cnv_driver_content)


# get reference map ids for sorted_direct_pharm_content
cnv_sorted_direct_pharm_content = content.reference_map_ids(
    cnv_sorted_direct_pharm_content, reference_map_dict)

# get reference map ids for sorted_pharm_content
cnv_sorted_pharm_content = content.reference_map_ids(
    cnv_sorted_pharm_content, reference_map_dict)

# get reference map ids for mechanistic_content and sort based n number of references
cnv_mechanistic_content = content.reference_map_ids(
    cnv_mechanistic_content, reference_map_dict)
cnv_mechanistic_content = prior.reference_sort(cnv_mechanistic_content)


##############################################################################################

# COUNT NUMBERS


# Get number of oncogenes
num_oncogene = content.oncogenes_count(cnv_driver_content)

# Get number of tumor suppressor genes
num_tsg = content.tsg_count(cnv_driver_content)

# Get number of genes, both onco and tumor suppressor genes.
num_onco_tsg = content.onco_tsg_count(cnv_driver_content)

# Get unknown type
num_unknown = content.unknown_count(cnv_driver_content)


##############################################################################################

# CONVERT CONENT INTO JSON

report_json = {}
driver_table_json = cnv_driver_content.to_dict("records")
direct_pharm_table = cnv_sorted_direct_pharm_content.to_dict("records")
pharm_table = cnv_sorted_pharm_content.to_dict("records")
mechanistic_table = cnv_mechanistic_content.to_dict("records")
appendix_variant_table = cnv_appendix_content.to_dict("records")
appendix_reference_table = cnv_reference_content.to_dict("records")

report_json["tag"] = "CNV"
report_json["num_oncogene"] = num_oncogene
report_json["num_tsg"] = num_tsg
report_json["num_onco_tsg"] = num_onco_tsg
report_json["num_unknown"] = num_unknown

report_json["driver_table"] = driver_table_json
report_json["direct_pharm_table"] = direct_pharm_table
report_json["pharm_table"] = pharm_table
report_json["mechanistic_drug_table"] = mechanistic_table
report_json["appendix_variant_table"] = appendix_variant_table
report_json["appendix_reference_table"] = appendix_reference_table


with open(output_json, "w") as f:
    json.dump(report_json, f, indent=4, ensure_ascii=False)
