"""
Script to perform post-filtering using tumor sample and different time points
"""

import argparse
import numpy as np
import pandas as pd
import os.path
import io

# constants
header_names_gsvar = {
    "chr": "#chr",
    "start": "start",
    "end": "end",
    "ref": "ref",
    "obs": "obs",
    "pval": "p-value",
    "multi AF": "m_AF",
    "multi REF count": "m_REF",
    "multi ALT count": "m_ALT",
    "Strand": "Strand",
    "Homopolymer": "Homopolymer",
    "filter": "filter"
}
header_names_tsv = {
    "chr": "CHROM",
    "start": "POS",
    "ref": "REF",
    "obs": "ALT",
    "pval": "Pval",
    "multi AF": "Multi_AF",
    "multi REF count": "Multi_REF",
    "multi ALT count": "Multi_ALT",
    "Strand": "STRAND",
    "Homopolymer": "Homopolymer",
    "filter": "FILTER"
}


def parse_args():
    parser = argparse.ArgumentParser(description="Script to perform post-filtering using tumor sample and different time points.")
    parser.add_argument("input_files", type=str, help="Comma separated list of cfDNA GSvar/TSV files.")
    parser.add_argument("output_files", type=str, help="Comma separated list of cfDNA GSvar/TSV files.")
    parser.add_argument("--tumor_samples", type=str, help="Comma separated list of tumor-normal GSvar files.", default="")
    parser.add_argument("--keep_indels", action='store_true', help="Do not remove InDels from variant list.")
    parser.add_argument("--min_alt_count", type=int, help="Minimal number of hits in sample for a variant to be kept", default=2)
    parser.add_argument("--keep_homopolymers", action='store_true', help="Do not remove homopolymers from variant list.")
    parser.add_argument("--min_depth", type=int, help="Minimal depth in each sample for a variant to be kept", default=1000)
    parser.add_argument("--min_strand_count", type=int, help="Minimal AC on each strtand in each sample for a variant to be kept", default=2)
    parser.add_argument("--keep_tumor_filter", action='store_true', help="Keep variants which has entries in the tumor filter column.")
    parser.add_argument("--keep_outliers", action='store_true', help="Do not remove variants which differ more than 3 stddev from mean.")
    parser.add_argument("--log_file", type=str, help="File path to log file.", default="filter_variants.log")

    args = parser.parse_args()
    return args


def high_pop_filter(row):
    if row["gnomAD"] > 0.05 or row["NGSD_hom"] == "n/a (AF>5%)":
        return True
    return False


def off_target_filter(row, tumors):
    n_off_target = 0
    for tumor in tumors:
        if "off-target" in row[tumor] or row[tumor] == "NOT_CALLED":
            n_off_target += 1
    return n_off_target == len(tumors)


def load_sample_file(path):
    # skip comment lines
    skip_rows = 0
    comments = []
    with open(path, 'r') as file:
        for line in file:
            if line.startswith("##"):
                comments.append(line)
                skip_rows += 1
            else:
                break

    # read TSV file
    tsv = pd.read_csv(path, sep='\t', header=0, skiprows=skip_rows)

    # remove previous annotations
    if "Tumor_Filter" in tsv.columns:
        tsv = tsv.drop("Tumor_Filter", axis=1)
    if "Post_Filter" in tsv.columns:
        tsv = tsv.drop("Post_Filter", axis=1)

    return comments, tsv


def parse_tumor_filter_column(path, sample):
    # read tumor-normal file
    _, tumor_gsvar = load_sample_file(path)
    tumor_gsvar["variant"] = tumor_gsvar["#chr"] + ":" + tumor_gsvar["start"].astype(str) + " " + tumor_gsvar["ref"] + ">" + tumor_gsvar["obs"]
    tumor_gsvar = tumor_gsvar.set_index("variant")
    tumor_gsvar["filter"].fillna("PASS", inplace=True)
    tumor_gsvar.loc[tumor_gsvar.apply(high_pop_filter, axis=1), "filter"] = tumor_gsvar["filter"] + ";HIGH_POP_AF"
    tumor_gsvar = tumor_gsvar.rename(columns={"filter": sample})
    return tumor_gsvar[sample]


def main():
    # get parameter
    args = parse_args()

    file_extension = ""

    cfdna_input_file_paths = args.input_files.split(',')
    cfdna_output_file_paths = args.output_files.split(',')
    if len(cfdna_input_file_paths) != len(cfdna_output_file_paths):
        raise ValueError("Number  off input and output files differ! ")
    file_types = set([os.path.basename(path).split('.')[-1] for path in cfdna_input_file_paths])
    if len(file_types) > 1:
        raise ValueError("File types of input cfDNA files differ! (" + ", ".join(file_types) + ")")
    file_type = file_types.pop().strip().lower()
    if file_type == "tsv":
        file_extension = "tsv"
        header_names = header_names_tsv
    elif file_type == "gsvar":
        file_extension = "GSvar"
        header_names = header_names_gsvar
    else:
        raise ValueError("Invalid file type '" + file_type + "'!")
    ext_offset = -(len(file_extension) + 1)

    # check if file exists
    for path in cfdna_input_file_paths:
        if not os.path.exists(path) or not os.path.isfile(path):
            raise FileNotFoundError(file_extension.upper() + " file '" + path + "' does not exist!")

    cfdna_samples = [os.path.basename(path)[:ext_offset] for path in cfdna_input_file_paths]
    print("cfDNA samples: " + ",".join(cfdna_samples))

    combined_gsvar_filter = pd.DataFrame()
    if args.tumor_samples.strip() != "":
        tn_sample_file_paths = args.tumor_samples.split(',')
        tn_samples = [os.path.basename(path)[:ext_offset] for path in tn_sample_file_paths]
        # check if file exists
        for path in tn_sample_file_paths:
            if not os.path.exists(path) or not os.path.isfile(path):
                raise FileNotFoundError("Tumor-normal GSvar file '" + path + "' does not exist!")
        print("tumor-normal samples: " + ",".join(tn_samples))

        # parse tumor files

        for path, tumor_normal in zip(tn_sample_file_paths, tn_samples):
            tumor_gsvar = parse_tumor_filter_column(path, tumor_normal)
            combined_gsvar_filter = pd.concat([combined_gsvar_filter, tumor_gsvar], axis=1, sort=True)
        # collapse filter columns:
        combined_gsvar_filter.fillna("NOT_CALLED", inplace=True)
        combined_gsvar_filter = combined_gsvar_filter[~(combined_gsvar_filter == "PASS").any(axis=1)]
        combined_gsvar_filter["Tumor_Filter"] = combined_gsvar_filter.apply(";".join, axis=1)
        if combined_gsvar_filter.size > 0:
            # remove off-target filter, if not off-target in all tumor samples
            mask = combined_gsvar_filter.apply(lambda x: off_target_filter(x, tn_samples), axis=1)
            combined_gsvar_filter.loc[~mask, "Tumor_Filter"] = combined_gsvar_filter["Tumor_Filter"].str.replace("off-target;", "")
        # combined_gsvar_filter = combined_gsvar_filter["Tumor_Filter"]
    else:
        combined_gsvar_filter["Tumor_Filter"] = []

    # perform post filtering
    combined_dataset = pd.DataFrame()
    log = ""
    initialized = False
    single_sample_buffer = {}

    for sample, path, in zip(cfdna_samples, cfdna_input_file_paths):
        print(sample + ":")
        log += sample + "\t"

        # determine file path for combined output file (always first sample of time series)
        if not initialized:
            output_file_path_combined = cfdna_output_file_paths[0][:ext_offset] + "_combined." + file_extension
        if not os.path.isfile(path):
            print("Error: Monitoring file " + sample + "_monitoring.tsv is missing! skipping sample!")
            continue

        # load sample VC file:
        comments, data = load_sample_file(path)

        # add index column
        data["variant"] = data[header_names["chr"]] + ":" + data[header_names["start"]].astype(str) + " " + data[header_names["ref"]] + ">" + data[header_names["obs"]]
        data = data.set_index("variant")

        # add tumor filter column
        data = pd.merge(data, combined_gsvar_filter["Tumor_Filter"].to_frame(), left_index=True, right_index=True, how='left')

        log += str(len(data.index)) + "\t"

        # filter data

        # add post-filter column:
        data["Post_Filter"] = ""

        # remove 0AC-SNVs
        if len(data.index) > 0:
            n_rows = len(data.index)
            mask = (data[header_names["multi ALT count"]] < args.min_alt_count)
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "LOW_mALT_COUNT;"
            print("Filter variants with less than " + str(args.min_alt_count) + "... \t" + str(n_rows - len(data[mask].index))
                  + " variant(s) removed. ")

        # other filters (indel, homoploymer, depth, extreme strand bias)
        # remove indels
        if not args.keep_indels and len(data.index) > 0:
            n_rows = len(data.index)
            mask = ((data[header_names["ref"]].str.replace('-', '').str.len() != 1) | (data[header_names["obs"]].str.replace('-', '').str.len() != 1))
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "INDEL;"
            print("Remove InDels... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # remove homopolymers
        if not args.keep_homopolymers and len(data.index) > 0:
            n_rows = len(data.index)
            mask = ((data[header_names["Homopolymer"]] == 1) | (data[header_names["Homopolymer"]] == "true"))
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "HOMOPOLYMER;"
            # data = data[data["Homopolymer"] == 0]
            print("Remove homopolymers... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # filter by depth
        if len(data.index) > 0:
            n_rows = len(data.index)
            mask = ((data[header_names["multi REF count"]] + data[header_names["multi ALT count"]]) < args.min_depth)
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "LOW_DEPTH;"
            print("Filter low depth variants... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # filter by strand bias
        if len(data.index) > 0:
            n_rows = len(data.index)
            data[header_names["Strand"]] = data[header_names["Strand"]].astype("object").fillna("0-0-0-0")
            mask = ((data[header_names["Strand"]].str.split('-', expand=True)[0].astype(int) < args.min_strand_count)
                    | (data[header_names["Strand"]].str.split('-', expand=True)[1].astype(int) < args.min_strand_count))
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "STRAND_BIAS;"
            print("Filter strand biased variants... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # filter by tumor filter column
        if not args.keep_tumor_filter and len(data.index) > 0:
            n_rows = len(data.index)
            mask = (data["Tumor_Filter"].isna() == False)
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "TUMOR_FILTER;"
            # data = data[data["Tumor_Filter"].isna()]
            print("Remove varaints with entries in tumor filter column... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # remove outliers
        if not args.keep_outliers and len(data.index) > 0:
            n_rows = len(data.index)
            std = data[header_names["multi AF"]].std()
            mean = data[header_names["multi AF"]].mean()
            mask = ((data[header_names["multi AF"]] < (mean - 3 * std)) | (data[header_names["multi AF"]] > (mean + 3 * std)))
            data.loc[mask, "Post_Filter"] = data["Post_Filter"] + "OUTLIER;"
            print("Remove outliers... \t" + str(n_rows - len(data[mask].index)) + " variant(s) removed. ")

        # buffer output file (to allow post-filter by combined file)
        print("Buffer single sample output...")
        single_sample_buffer[sample] = data
        single_sample_buffer[sample + "_header"] = comments
        log += str(len(data[data["Post_Filter"] == ""].index)) + "\n"

        # combine with previous
        if not initialized:
            combined_dataset = pd.DataFrame({"variant": data.index})
            combined_dataset.set_index("variant", inplace=True)
            combined_dataset = pd.concat([combined_dataset, data[[]]], axis=1)
            initialized = True

        # remove not needed columns
        data = data[[header_names["multi REF count"], header_names["multi ALT count"], header_names["multi AF"], header_names["pval"], header_names["Strand"],
                     header_names["Homopolymer"], header_names["filter"]]]

        combined_dataset = pd.merge(combined_dataset, data.add_prefix(sample + "_"), left_index=True, right_index=True, how="outer")

        print("-------------------------------------------")

    # extract variant info from index
    tmp = combined_dataset.index.to_series().str.split(' ', 1, expand=True)
    pos = tmp[0]
    change = tmp[1]
    tmp = pos.str.split(':', 1, expand=True)
    chr = tmp[0]
    start = tmp[1]
    combined_dataset.insert(0, header_names["chr"], chr)
    combined_dataset.insert(1, header_names["start"], start)
    combined_dataset[header_names["start"]] = combined_dataset[header_names["start"]].astype("int64")
    tmp = change.str.split('>', 1, expand=True)
    ref = tmp[0]
    obs = tmp[1]
    combined_dataset.insert(2, header_names["ref"], ref)
    combined_dataset.insert(3, header_names["obs"], obs)

    if file_extension == "GSvar":
        end = combined_dataset[header_names["start"]] + combined_dataset[header_names["ref"]].str.len() - 1
        combined_dataset.insert(2, header_names["end"], end)

    # Apply filter based on combined table

    # check if combined column contains at least 3 informative variants
    alt_columns = []
    informative = {}
    for column in combined_dataset.columns:
        if column.endswith(header_names["multi ALT count"]):
            alt_columns.append(column)
    multi_alt_prod = combined_dataset.loc[:, alt_columns].fillna(0).prod(axis=1)
    n_informative_vars = multi_alt_prod[multi_alt_prod > 0].size

    if n_informative_vars >= 3:
        for sample in cfdna_samples:
            informative[sample] = True
    else:
        for sample in cfdna_samples:
            # check if sample contains 3 distinct variants
            buffer = single_sample_buffer[sample]
            buffer = buffer[(buffer["Post_Filter"] == "") & (buffer[header_names["multi ALT count"]] > 1e-6)]
            n_informative_vars = buffer.shape[0]
            informative[sample] = (n_informative_vars >= 3)

    # mark non-informative variants
    multi_alt_sum = combined_dataset.loc[:, alt_columns].fillna(0).sum(axis=1)
    for sample in cfdna_samples:
        if not informative[sample]:
            continue
        mask = (multi_alt_sum == 0)

        buffer = single_sample_buffer[sample]
        mask = mask[mask.index.intersection(buffer.index)]
        test = buffer.loc[mask, "Post_Filter"]
        buffer.loc[mask, "Post_Filter"] = buffer["Post_Filter"] + "NON-INFORMATIVE;"
        single_sample_buffer[sample] = buffer

    # check if all variants have at least 1000x multi-depth
    alt_columns = []
    ref_columns = []
    for column in combined_dataset.columns:
        if column.endswith(header_names["multi ALT count"]):
            alt_columns.append(column)
        if column.endswith(header_names["multi REF count"]):
            ref_columns.append(column)
    mask = []
    for index, row in combined_dataset.iterrows():
        for alt_column, ref_column in zip(alt_columns, ref_columns):
            if (row[alt_column] + row[ref_column]) < 1000:
                mask.append(index)
                break
    # apply filter to sample files
    for sample in cfdna_samples:
        buffer = single_sample_buffer[sample]
        indices = buffer.index.intersection(mask)
        buffer.loc[indices, "Post_Filter"] = buffer["Post_Filter"] + "LOW_MULTIDEPTH;"
        single_sample_buffer[sample] = buffer

    # calculate mean and median
    combined_dataset.loc["mean"] = combined_dataset.mean(numeric_only=True)
    combined_dataset.loc["mean", header_names["chr"]] = "mean"
    combined_dataset.loc["mean", header_names["start"]] = -1
    if file_extension == "GSvar":
        combined_dataset.loc["mean", header_names["end"]] = -1
    combined_dataset.loc["median"] = combined_dataset.median(numeric_only=True)
    combined_dataset.loc["median", header_names["chr"]] = "median"
    combined_dataset.loc["median", header_names["start"]] = -1
    if file_extension == "GSvar":
        combined_dataset.loc["median", header_names["end"]] = -1

    # store filtered monitoring file:
    print("Writing output (single sample)...")
    for sample, output_file_path in zip(cfdna_samples, cfdna_output_file_paths):
        # skip if output file name is empty
        if output_file_path.strip() == "":
            continue
        with open(output_file_path, "w") as file:
            file.write("".join(single_sample_buffer[sample + "_header"]))
        single_sample_buffer[sample].to_csv(output_file_path, float_format="%.8g", sep="\t", index=False, mode='a')

    # store in file
    print("Writing output (combined)...")
    if cfdna_output_file_paths[0].strip() != "":
        with open(output_file_path_combined, "w") as file:
            file.write("".join(single_sample_buffer[cfdna_samples[0] + "_header"]))
        combined_dataset.to_csv(output_file_path_combined, float_format="%.8g", sep="\t", index=False, mode='a')

    with open(args.log_file, 'a') as file:
        file.write(log)

    print("finished!")


if __name__ == '__main__':
    main()

