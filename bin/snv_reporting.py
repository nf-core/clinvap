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
import argparse
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


########################
### Read Process SNV ###
########################

# PROCESS VCF

# Read in and process VCF file
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
    "HGNC_ID", "gene", "one_letter_repr", "Consequence", "vaf", "var_type"]]


#############################################################################################################

########################
### Read Process CNV ###
########################

# Process CNV to get info that SNV analysis depends on (in section combined pharmacogenomics,
# the combination partner might be a CNV)

if input_cnv:
    df_cnv = cnv.read_cnv(input_cnv)
    cnv_flag = True
else:
    # if cnv is empty, this is placeholder for filter_pair_variant function.
    cnv_gene_tuples = [()]
    cnv_flag = False


if cnv_flag and not df_cnv.empty:
    df_cnv_processed = cnv.parse_cnv(df_cnv)

    # Get gene, type tuples reqired in pair variant filtering in process_pharmacogenomics part.
    cnv_gene_tuples = query.get_gene_list(df_cnv_processed)[2]

    # fields needed from cnv dataframe
    from_cnv = df_cnv_processed[["HGNC_ID",
                                 "gene", "type", "copy_number", "effect"]]

elif cnv_flag and df_cnv.empty:
    print("CNV file is empty")
    cnv_flag = False
    # we need empty cnv_gene_tuple list as a place holder for upcaming parts (filter_gene_pair function)
    cnv_gene_tuples = [()]


#############################################################################################################

############
### SNVs ###
############


# QUERY DATABASE


#  Get driver gene annotations for snv and small indels
driver_ann_mutation = query.get_driver_annotation(main_mvld, KNOWLEDGEBASE)[0]
driver_ann_indel = query.get_driver_annotation(main_mvld, KNOWLEDGEBASE)[1]


# Query database - Get direct drug targets
snv_direct_variant_targets = query.ann_pharm_variant(
    main_mvld, genome, KNOWLEDGEBASE, "MUT")


# Filter direct drug targets for variant_class
mut_direct_variant_targets = query.filter_main_variant_class(
    snv_direct_variant_targets, "MUT")


# Query database - Get drugs targeting affected gene
snv_gene_targets = query.ann_pharm_gene(main_mvld, KNOWLEDGEBASE, genome)

# Filter direct drug targets for variant_class
mut_gene_targets = query.filter_main_variant_class(snv_gene_targets, "MUT")

# Query database - Get mechanistic drug targets, cancer drugs (later - remove the ones found in direct targets and indirect targets)
mechanistic_drug_targets = query.get_mechanistic_drugs(
    main_mvld, KNOWLEDGEBASE)

if (not mechanistic_drug_targets) and (not driver_ann_mutation) and (not driver_ann_indel) and (not mut_gene_targets) and (not mut_direct_variant_targets):
    handle.empty_database_result()


# Check if mechanistic drug targets are empty or not.
if mechanistic_drug_targets:
    mechanistic_flag = True
else:
    mechanistic_flag = False
    handle.empty_mechanistic("SNV")


# Check if driver gene annotation dictionary is empty or not.
if driver_ann_mutation:
    driver_mut_flag = True
else:
    driver_mut_flag = False

if driver_ann_indel:
    driver_indel_flag = True
else:
    driver_indel_flag = False


############################################################################

# PROCESS DATABASE RESULTS

############################################################################
# DRIVER ANNOTATION

# Update driver dict
# add tumor string to driver dictionary and remove tumor_type list of dictionaries (same function doing both)
for dt in [driver_ann_mutation, driver_ann_indel]:
    for key in dt:
        helper.add_tumor_string(key, dt)


if driver_mut_flag:
    df_driver_ann_mutation = driver.get_driver_df(driver_ann_mutation)
else:
    df_driver_ann_mutation = handle.empty_driver_annotation("SNVs")


if driver_indel_flag:
    df_driver_ann_indel = driver.get_driver_df(driver_ann_indel)
else:
    df_driver_ann_indel = handle.empty_driver_annotation("INDELs")

# Concat snv and indel driver annotations to get driver content on those.
df_driver_ann = pd.concat([df_driver_ann_mutation, df_driver_ann_indel],
                          ignore_index=True)
                        #   .drop(columns=["index"])

# Start building up driver content. References will be mapped later.
driver_ann = df_driver_ann[[
    "hgnc_id", "driver_role", "reference_id", "tumor_list", "db_tumor_repr", "var_type"]]

driver_content = pd.merge(from_mlvd, driver_ann, how="right",
                          left_on=["HGNC_ID", "var_type"], right_on=["hgnc_id", "var_type"]).drop(columns=["HGNC_ID"])

# Get references
ref_driver = content.get_references(df_driver_ann)

###############################################################################################################

# DIRECT VARIANT TARGETS MUT

if mut_direct_variant_targets:
    empty_variant_targets = False

    # DIRECT PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from variant targets
    direct_pharmacogenomics = pharma.get_subset_dict(
        mut_direct_variant_targets, "pharmacogenomics_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    # Direct pharmacogenomics effect dictionary to dataframe
    df_direct_pharmacogenomics = pharma.dict_to_dataframe(
        direct_pharmacogenomics, "pharmacogenomics_therapeutics")

    if df_direct_pharmacogenomics.empty:
        df_direct_pharmacogenomics = handle.empty_dataframe_direct_pharmacogenomics(
            "SNVs")

    # DIRECT PHARMACOGENOMICS FOR COMBINED VARIANTS

    # seperate non-empty pharmacogenomics combined variant therapeutics from variant targets
    direct_pharmacogenomics_combined = pharma.get_subset_dict(
        mut_direct_variant_targets, "pharmacogenomics_combined_variants_therapeutics")

    # Filter based on the distrupted pair gene
    direct_pharmacogenomics_combined = pharma.filter_gene_pair(
        snv_gene_list, cnv_gene_tuples, direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")

    # Direct pharmacogenomics combined variant effect dictionary to dataframe
    df_direct_pharmacogenomics_combined = pharma.dict_to_dataframe(
        direct_pharmacogenomics_combined, "pharmacogenomics_combined_variants_therapeutics")

    if df_direct_pharmacogenomics_combined.empty:
        pass

    # DIRECT ADVERSE EFFECTS

    # seperate non-empty adverse effect therapeutics from variant targets
    direct_adverse = pharma.get_subset_dict(
        mut_direct_variant_targets, "adverse_effect_therapeutics")

    # Direct adverse effect dictionary to dataframe
    df_direct_adverse = pharma.dict_to_dataframe(
        direct_adverse, "adverse_effect_therapeutics")

    if df_direct_adverse.empty:
        df_direct_adverse = handle.empty_dataframe_direct_adverse()
else:
    df_direct_pharmacogenomics, df_direct_adverse = handle.empty_direct_variant_annotation(
        "SNV")
    empty_variant_targets = True


#############################################################################################################

# GENE TARGETS MUT

if mut_gene_targets:
    empty_gene_targets = False

    # PHARMACOGENOMICS

    # seperate non-empty pharmacogenomics therapeutics from gene targets
    pharmacogenomics_therapeutics = pharma.get_subset_dict(
        mut_gene_targets, "pharmacogenomics_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(
        pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    # Pharmacogenomics effects to dataframe
    df_pharmacogenomics = pharma.dict_to_dataframe(
        pharmacogenomics_therapeutics, "pharmacogenomics_therapeutics")

    if df_pharmacogenomics.empty:
        df_pharmacogenomics = handle.empty_dataframe_pharmacogenomics("SNVs")

    # PHARMACOGENOMICS COMBINED VARIANTS

    # seperate non-empty pharmacogenomics combined variant therapeutics from gene targets
    pharmacogenomics_combined_variants_therapeutics = pharma.get_subset_dict(
        mut_gene_targets, "pharmacogenomics_combined_variants_therapeutics")

    # Filter based on the distrupted pair gene
    pharmacogenomics_combined_variants_therapeutics = pharma.filter_gene_pair(
        snv_gene_list, cnv_gene_tuples, pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")

    # merge tumor_type list into a list, remove "evidence_statement" and "rating"
    pharma.apply_add_tumor_string(pharmacogenomics_combined_variants_therapeutics,
                                  "pharmacogenomics_combined_variants_therapeutics")

    df_pharmacogenomics_combined = pharma.dict_to_dataframe(
        pharmacogenomics_combined_variants_therapeutics, "pharmacogenomics_combined_variants_therapeutics")

    if df_pharmacogenomics_combined.empty:
        df_pharmacogenomics_combined = handle.empty_dataframe_pharmacogenomics_combined(
            "SNVs")

    # ADVERSE EFFECTS

    # seperate non-empty adverse effect therapeutics from gene targets
    adverse_effect_therapeutics = pharma.get_subset_dict(
        mut_gene_targets, "adverse_effect_therapeutics")

    df_adverse = pharma.dict_to_dataframe(
        adverse_effect_therapeutics, "adverse_effect_therapeutics")

    if df_adverse.empty:
        df_adverse = handle.empty_dataframe_adverse()

    # ADVERSE EFFECTS COMBINED VARIANTS

    # seperate non-empty adverse effect combined variants from gene targets
    adverse_effect_combined_variants_therapeutics = pharma.get_subset_dict(
        mut_gene_targets, "adverse_effect_combined_variants_therapeutics")

    # filter based on the distrupted pair
    adverse_effect_combined_variants_therapeutics = pharma.filter_gene_pair(
        snv_gene_list, cnv_gene_tuples, adverse_effect_combined_variants_therapeutics, "adverse_effect_combined_variants_therapeutics")

    # Adverse effcts combination variants dictionary to dataframe
    df_adverse_combined = pharma.dict_to_dataframe(
        adverse_effect_combined_variants_therapeutics, "adverse_effect_combined_variants_therapeutics")

    if df_adverse_combined.empty:
        df_adverse_combined = handle.empty_dataframe_adverse_combined()
else:
    df_pharmacogenomics, df_pharmacogenomics_combined, df_adverse, df_adverse_combined = handle.empty_variant_annotation(
        "SNV")
    empty_gene_targets = True


if empty_gene_targets and empty_variant_targets:
    skip_pharmacodynamics_content = True
else:
    skip_pharmacodynamics_content = False


########################################################################################################
if not skip_pharmacodynamics_content:

    # CREATE PHARMACOGENOMICS TABLES CONTENT

    COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source",
               "variant_drug_association", "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "info"]
    COMB_COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source", "variant_drug_association",
                    "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "info", "info_pair"]

    # Content from direct pharmacogenomics effect
    direct_pharm_content = content.snv_get_content(
        df_direct_pharmacogenomics, from_mlvd, COLUMNS, direct=1, pharm=1, notcombined=1)

    # Remove direct results covered by df_pharmacogenomics
    pharm = pd.concat(
        [df_pharmacogenomics, df_direct_pharmacogenomics]).drop_duplicates(keep=False)

    # Content from pharmacogenomics effect
    pharm_content = content.snv_get_content(
        pharm, from_mlvd, COLUMNS, indirect=1, pharm=1, notcombined=1)

    # Content from pharmacogenomics variant combination effect
    pharm_comb_content = content.snv_get_content(
        df_pharmacogenomics_combined, from_mlvd, COMB_COLUMNS, combined=1, pharm=1, indirect=1)

    # move same variant for both genes, if there are any, content to direct pharm table.
    cond = (pharm_comb_content["info"] == "Same gene, same variant, same consequence.") & (
        pharm_comb_content["info_pair"] == "Same gene, same variant, same consequence.")
    rows = pharm_comb_content.loc[cond, :]
    if not rows.empty:
        direct_pharm_content = direct_pharm_content.append(
            rows, ignore_index=True)
        pharm_comb_content.drop(rows.index, inplace=True)

    # Add match level to see whether the results are same gene, same, variant,same consequence or whether the
    # variant and consequence differs in the inut than the database result. For detais see documentation

    sorted_direct_pharm_content = prior.compatibility_sort(
        direct_pharm_content)

    # Join references of the same rows (collapsing same rows with different reference ids)
    cols = ["gene", "drug_name", "variant_drug_association",
            "tumor_list", "db_tumor_repr", "hgnc_id", "variant", "variant_type", "match_level"]
    sorted_direct_pharm_content = prior.unnest_ref(
        sorted_direct_pharm_content, cols)

    pharm_content = pd.concat([pharm_content, pharm_comb_content])
    sorted_pharm_content = prior.compatibility_sort(pharm_content)
    sorted_pharm_content = prior.unnest_ref(sorted_pharm_content, cols)

    # Get references
    ref_direct_pharm = content.get_references(direct_pharm_content)

    ref_pharm = content.get_references(pharm_content)

    ##############################################################################################

    # CREATE ADVERSE EFFECT TABLE CONTENTS

    ADV_COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source",
                   "variant_drug_association", "hgnc_id", "variant", "variant_type", "info"]

    ADV_COMB_COLUMNS = ["gene", "drug_name", "evidence_level", "reference_id", "reference_source", "variant_drug_association",
                        "hgnc_id", "variant", "variant_type", "info", "info_pair"]

    # Content from direct adverse effects table
    direct_adverse_content = content.snv_get_content(
        df_direct_adverse, from_mlvd, ADV_COLUMNS, direct=1, adverse=1, notcombined=1)

    # remove direct results covered by df_adverse
    adverse = pd.concat([df_direct_adverse, df_adverse]
                        ).drop_duplicates(keep=False)

    # Content from adverse effects table
    adverse_content = content.snv_get_content(
        adverse, from_mlvd, ADV_COLUMNS, adverse=1, indirect=1, notcombined=1)

    # Content from adverse effects combined variants
    adverse_comb_content = content.snv_get_content(
        df_adverse_combined, from_mlvd, ADV_COMB_COLUMNS, combined=1, adverse=1, indirect=1)

    # merge direct_adverse_content, adverse_content, adverse_comb_content
    adverse_content = pd.concat(
        [direct_adverse_content, adverse_content, adverse_comb_content]).drop_duplicates()

    # add merge level and sort based on levels
    cols = ["gene", "drug_name", "variant_drug_association",
            "hgnc_id", "variant", "variant_type", "match_level"]
    sorted_adverse_content = prior.compatibility_sort(adverse_content)
    sorted_adverse_content = prior.unnest_ref(sorted_adverse_content, cols)

    ref_adverse = content.get_references(adverse_content)

else:
    sorted_direct_pharm_content = handle.empty_pharmacodynamics_place_holder()
    sorted_pharm_content = handle.empty_pharmacodynamics_place_holder()
    sorted_adverse_content = handle.empty_pharmacodynamics_place_holder(
        adverse=1)
    ref_pharm = handle.empty_reference_place_holder()
    ref_direct_pharm = handle.empty_reference_place_holder()
    ref_adverse = handle.empty_reference_place_holder()


#############################################################################################################


if mechanistic_flag:

    # Update mechanistic drug targets dict
    mechanistic_drug_targets = mechanistic.reshape(
        mechanistic_drug_targets, {})

    # collapse info from sources into one representation
    temp = {}
    for key, annotation in mechanistic_drug_targets.items():
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
    mechanistic_drug_targets = temp

    # filter not-approved drugs from mechanistic drugs
    mechanistic_drugs_filtered = mechanistic.filter_approval_info(
        mechanistic_drug_targets)

    # seperate investigational and approved mechanistic drugs
    approved_mechanistic_drugs = mechanistic.divide_approved_investigational(
        mechanistic_drugs_filtered)[0]

    investigational_mechanistic_drugs = mechanistic.divide_approved_investigational(
        mechanistic_drugs_filtered)[1]

    # Mechanistic drug targets information to dataframe
    if approved_mechanistic_drugs:
        main_mechanistic = approved_mechanistic_drugs
    else:
        handle.empty_approved_mechanistic()
        main_mechanistic = investigational_mechanistic_drugs

    df_mechanistic_targets = mechanistic.dict_to_dataframe(main_mechanistic)
else:
    df_mechanistic_targets = pd.DataFrame(
        columns=['drugbank_id', 'drug_name', 'approval_status', 'source_name', 'reference_id', 'reference_source', 'source_pmid', 'hgnc_id'])


# CREATE MECHANISTIC DRUG TARGETS TABLE CONTENT


if not skip_pharmacodynamics_content:
    # remove the content found in pharmacogenomics therapeutics - bunlardan biri bossa bura patlar mi?
    all_pharm_drugs = pd.concat([direct_pharm_content[["drug_name", "hgnc_id"]], pharm_content[[
        "drug_name", "hgnc_id"]], pharm_comb_content[["drug_name", "hgnc_id"]]]).drop_duplicates()

    cond = df_mechanistic_targets["drug_name"].isin(
        all_pharm_drugs["drug_name"]) & df_mechanistic_targets["hgnc_id"].isin(all_pharm_drugs["hgnc_id"])
    mechanistic = df_mechanistic_targets.loc[~cond]
else:
    mechanistic = df_mechanistic_targets

MECH_COLS = ["gene", "drug_name", "approval_status",
             "reference_id", "reference_source", "hgnc_id"]
mechanistic_content = content.get_mechanistic_content(
    mechanistic, from_mlvd, MECH_COLS)


ref_mechanistic = content.get_references(mechanistic_content)


##############################################################################################

# CREATE APPENDIX TABLES

# appendix table for all variants
appendix_content = content.appendix_variants(main_mvld, "SNV")


# Appendix table for references


# concatanate reference dataframes
dataframes = [ref_driver, ref_direct_pharm,
              ref_pharm, ref_mechanistic, ref_adverse]

ref_df = functools.reduce(
    lambda left, right: pd.concat([left, right], ignore_index=True), dataframes)  # concatanate dataframes to create one reference dataframe


# create pubmed query string for eutils
pubmed_ids, ref_map_df = content.id_string(ref_df)
URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?retmode=json;db=pubmed;id={}".format(
    pubmed_ids)


# get API results, create dataframe with reference information
ref_info_df = content.get_ref_details(URL)

# merge reference map and reference info dataframes two get a complete reference map.
ref_map_complete = ref_map_df.merge(
    ref_info_df, how="left", left_on="reference_id", right_index=True)

reference_content, reference_map_dict = content.reference_table(
    ref_map_complete)

# Map reference ids to the dataframes and sort based on the number of the references they have

# get reference map ids for driver_content and sort based n number of references
# listing unknown type of drivers under known typed ones.
driver_content = content.reference_map_ids(driver_content, reference_map_dict)
driver_content = prior.driver_type_sort(driver_content)


# get reference map ids for sorted_direct_pharm_content
sorted_direct_pharm_content = content.reference_map_ids(
    sorted_direct_pharm_content, reference_map_dict)

# get reference map ids for sorted_pharm_content
sorted_pharm_content = content.reference_map_ids(
    sorted_pharm_content, reference_map_dict)

# get reference map ids for mechanistic_content and sort based n number of references
mechanistic_content = content.reference_map_ids(
    mechanistic_content, reference_map_dict)
mechanistic_content = prior.reference_sort(mechanistic_content)

# get reference map ids for sorted_adverse_content
sorted_adverse_content = content.reference_map_ids(
    sorted_adverse_content, reference_map_dict)


##############################################################################################

# COUNT NUMBERS


# Get number of nonsynonymous genes
num_nonsynonymous = content.nonsynonymous_count(main_mvld)

# Get number of oncogenes
num_oncogene = content.oncogenes_count(driver_content)

# Get number of tumor suppressor genes
num_tsg = content.tsg_count(driver_content)

# Get number of genes, both onco and tumor suppressor genes.
num_onco_tsg = content.onco_tsg_count(driver_content)

# Get unknown type
num_unknown = content.unknown_count(driver_content)


##############################################################################################

# CONVERT CONENT INTO JSON


report_json = {}
driver_table_json = driver_content.to_dict("records")
direct_pharm_table = sorted_direct_pharm_content.to_dict("records")
pharm_table = sorted_pharm_content.to_dict("records")
mechanistic_table = mechanistic_content.to_dict("records")
adverse_table = sorted_adverse_content.to_dict("records")
appendix_variant_table = appendix_content.to_dict("records")
appendix_reference_table = reference_content.to_dict("records")

report_json["tag"] = "SNV"
report_json["num_nonsynonymous"] = num_nonsynonymous
report_json["num_oncogene"] = num_oncogene
report_json["num_tsg"] = num_tsg
report_json["num_onco_tsg"] = num_onco_tsg
report_json["num_unknown"] = num_unknown

report_json["driver_table"] = driver_table_json
report_json["direct_pharm_table"] = direct_pharm_table
report_json["pharm_table"] = pharm_table
report_json["mechanistic_drug_table"] = mechanistic_table
report_json["adverse_table"] = adverse_table
report_json["appendix_variant_table"] = appendix_variant_table
report_json["appendix_reference_table"] = appendix_reference_table


with open(output_json, "w") as f:
    json.dump(report_json, f, indent=4, ensure_ascii=False)
