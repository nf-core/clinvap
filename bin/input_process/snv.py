# VALID ONLY FOR SNVs

import table_content.handling as handle
from cyvcf2 import VCF

import sys
import pandas as pd
import re


def read_vcf(inputfile):
    """Read vcf"""
    vcf = VCF(inputfile)
    return vcf


def get_fields(vcf):
    """Get annotation field names - CSQ comes from Ensembl VEP. It might not work for SV or CNV"""
    try:
        header = vcf.get_header_type(
            "CSQ")["Description"][51:].strip('"').split("|")
        return header
    except KeyError:
        print("Input file is not annotated. Run Ensembl VEP on input file. Terminating...")
        sys.exit()


def _strip_allele_comma(annotation_list):
    allele = annotation_list[0]
    stripped_allele = allele.strip(',')
    annotation_list[0] = stripped_allele
    return annotation_list


def get_annotation(variant, header):
    """Get variant annotation info"""
    annotation = variant.INFO.get('CSQ').split('|')
    pubmed = annotation[-1]
    del annotation[-1]
    list_annotation = [annotation[x: x + len(header) - 1]
                    for x in range(0, len(annotation), len(header) - 1)]
    [ls.insert(39, pubmed) for ls in list_annotation if len(ls) == 38]
    [_strip_allele_comma(ls) for ls in list_annotation]
    return list_annotation


def ngs_source(vcf):  # takes vcf object as argument
    """Function to get NGS pipeline name from the VCF header."""
    if vcf.contains("source"):
        source_info = vcf.get_header_type("source")
        source_name = source_info["source"]
    else:
        source_name = "null"
    return source_name


# vaf formula: https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md
def strelka_vaf_calculation(variant, tumor_sample_index):
    if variant.var_type == "snp":
        ref_counts_field = "{}U".format(variant.REF)
        alt_counts_field = "{}U".format("".join(variant.ALT))
        tier1_ref_counts = variant.format(ref_counts_field)[
            tumor_sample_index][0]
        tier1_alt_counts = variant.format(alt_counts_field)[
            tumor_sample_index][0]
    elif variant.var_type == "indel":
        tier1_ref_counts = variant.format("TAR")[tumor_sample_index][0]
        tier1_alt_counts = variant.format("TIR")[tumor_sample_index][0]
    else:
        return "."
    vaf = tier1_alt_counts / (tier1_ref_counts + tier1_alt_counts)
    return "%.2f" % vaf


def stratacaller_vaf_identification(variant):
    vaf = float("".join(variant.format("FREQ")).strip("%"))/100
    return "%.2f" % vaf


def other_ngs_vaf(variant):
    vaf = "."
    return vaf


# ##SAMPLE filed and IsTumor assumption is made here. VCF has to be checked for those even if it comes from strelka
def get_tumor_sample(vcf):
    """Function to get tumor sample id from the VCF header."""
    sample_tumor_dict = {}
    header_list = vcf.raw_header.split("\n")
    info_sample = [h for h in header_list if h.startswith("##SAMPLE")]
    for s in info_sample:
        sample_name = re.search("(?<=##SAMPLE=<ID=).+?(?=,)", s).group(0)
        tumor_state = re.search("(?<=IsTumor=).+?(?=,)", s).group(0)
        if tumor_state == "yes":
            sample_tumor_dict[sample_name] = tumor_state
        else:
            continue
    return sample_tumor_dict


# Only one tumor sample assumption is made. If more, fix is needed here.
def tumor_sample_format_index(vcf, source):
    """Function to return index of tumor sample name. The required Format filed of the tumor sample can be reached with the index. 
    It returns index if ##Sample field if present in VCF header. It returns index if sample names are Tumor, Normal. It returns "null" otherwise.
    """
    if source == "strelka":
        tumor_sample = get_tumor_sample(vcf)
        samples = vcf.samples
        if tumor_sample:
            for s in samples:
                if tumor_sample.get(s) == "yes":
                    tumor_sample_index = samples.index(s)
        elif "NORMAL" in samples and "TUMOR" in samples:
            tumor_sample_index = samples.index("TUMOR")
        else:
            tumor_sample_index = "null"
        return tumor_sample_index
    else:
        return "null"


# Read VCF file and create dataframe minium variant level data (MVLD) according to
# Ritter et al. (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0367-z)

def parse_vcf(inputfile):
    """Parse and process VCF, create dataframe"""
    complete_line_list = list()
    vcf = read_vcf(inputfile)
    fields = get_fields(vcf)  # Get fields of VCF
    source = ngs_source(vcf)
    tumor_sample_index = tumor_sample_format_index(vcf, source)
    for variant in vcf:
        if source == "strelka":
            if tumor_sample_index != "null":
                vaf = strelka_vaf_calculation(variant, tumor_sample_index)
            else:
                vaf = other_ngs_vaf(variant)
        elif "StrataCaller" in source:
            vaf = stratacaller_vaf_identification(variant)
        else:
            vaf = other_ngs_vaf(variant)
        chrom = variant.CHROM.replace("chr", "")
        if chrom != "X":
            chrom = int(chrom)
        start = variant.POS
        stop = variant.end
        vid = variant.ID
        ref = variant.REF
        qual = variant.QUAL
        vfilter = variant.FILTER
        alt = variant.ALT
        var_type = variant.var_type
        alt = ''.join(variant.ALT)
        location = "{}:{}_{}/{}".format(chrom, start, ref, alt)
        line = [chrom, start, stop, location,
                vid, ref, alt, qual, vfilter, vaf, var_type]
        annotation = get_annotation(variant, fields)
        for ann in annotation:
            for ln in line:
                ann.append(ln)
            complete_line_list.append(ann)
    header = fields + ["chr", "start", "end", "location",
                       "vid", "ref", "alt", "qual", "filter", "vaf", "var_type"]  # create VCF dataframe
    dataframe = pd.DataFrame(complete_line_list, columns=header)

    dataframe["dbSNP"] = dataframe["Existing_variation"].str.findall(
        "(rs.\d+)").apply(",".join)
    dataframe["COSMIC"] = dataframe["Existing_variation"].str.findall(
        "(COSM.\d+)").apply(",".join)

    if dataframe.empty:
        handle.empty_mvld()
    else:
        return dataframe


def get_high_moderate_effect(dataframe):
    """Filter mvld, create mvld dataframe for high and moderate variant effects, works for Ensembl VEP output"""
    dataframe = dataframe[dataframe["PICK"] == "1"]

    if len(dataframe.index != 0):
        first_cond = dataframe[(dataframe["IMPACT"] == "HIGH") | (
            dataframe["IMPACT"] == "MODERATE")]
        second_cond = first_cond[~((first_cond["SIFT"] == "tolerated") & (
            first_cond["PolyPhen"] == "benign"))]
        third_cond = second_cond[~((second_cond["SIFT"] == "tolerated_low_confidence") & (
            second_cond["PolyPhen"] == "benign"))]
        return third_cond
    else:
        print("No variants found that passed default filters.")
        return pd.DataFrame()


def get_modifier_effect(dataframe):
    """Create mvld dataframe for variants that has modifying effect, works for Ensembl VEP output"""
    if len(dataframe.index != 0):
        cond = dataframe[(dataframe["IMPACT"] == "MODIFIER")]
        return cond
    else:
        print("No variants found that passed default filters.")
        return pd.DataFrame()


ONE_LETTER = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L",
              "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Asx": "B", "Glx": "Z", "Xaa": "X", "Ter": "*"}


def one_letter_repr(dataframe):
    """Function to represent HGVSp with one letter aminoacids"""
    dataframe["aa"] = dataframe["HGVSp"].str.replace("ENSP.*p\.", "")
    dataframe["aa_1"] = dataframe["aa"].str.extract("([a-zA-Z]{3}(?=\d))")
    dataframe["aa_number"] = dataframe["aa"].str.extract(
        "((?<=[a-zA-Z]{3}).*(?=[a-zA-Z]{3}))")
    dataframe["aa_2"] = dataframe["aa"].str.extract("((?<=\d)[a-zA-Z]{3})")
    dataframe["aa_1"] = dataframe["aa_1"].map(ONE_LETTER)
    dataframe["aa_2"] = dataframe["aa_2"].map(ONE_LETTER)
    cols = ["aa_1", "aa_number", "aa_2"]
    dataframe["one_letter_repr"] = dataframe[cols].apply(
        lambda row: "".join(row.values.astype(str)), axis=1)
    dataframe = dataframe.drop(["aa", "aa_1", "aa_2", "aa_number"], axis=1)
    return dataframe
