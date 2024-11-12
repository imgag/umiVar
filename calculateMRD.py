"""
Alternative MRD calculation
Calculates the MRD based on the difference between monitoring AF and background AF
"""
import argparse
import scipy
from scipy import stats
import numpy as np
import pandas as pd
import os
from functools import cmp_to_key


def parse_args():
    parser = argparse.ArgumentParser(description="Script to calculate MRD based on background noise.")
    parser.add_argument("folder", type=str, help="umiVar folder containing dedup TSV files.")
    parser.add_argument("panel", type=str, help="VCF containing the cfDNA panel.(commas seperated if multiple)")
    parser.add_argument("output", type=str, help="Output TSV containing MRD value.")
    parser.add_argument("variants", type=str, help="Output TSV containing read counts and p-values for each monitoring variants.")
    parser.add_argument("--max_af", type=float, help="Maximum allele frequency for variants to include in the computation. (>1 to disable)", default=0.2, required=False)
    parser.add_argument("--keep_gonosomes", action="store_true", help="Do not remove gonosomes")
    parser.add_argument("--blacklist", type=str, help="VCF of variants which should be excluded for the background error rate", default="")
    parser.add_argument("--remove_off_target", action="store_true", help="Remove tumor off-target variants")
    parser.add_argument("--keep_indels", action="store_true", help="Do not remove InDels.")

    args = parser.parse_args()
    return args


def parse_table(file_path):
    table = pd.read_csv(file_path, delimiter='\t')

    # set index
    table["index"] = table["CHROM"] + ":" + table["POS"].astype(str)
    table = table.set_index("index")

    # cleanup
    table = table[["REF", "DP", "DP_HQ", "A", "C", "T", "G", "INS", "DEL"]]

    return table


def parse_panel(file_paths):
    monitoring_snps = []
    id_snps = []
    off_target = []
    for file_path in file_paths.split(','):
        with open(file_path, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                if line.strip() == "":
                    continue
                split_line = line.split("\t")
                if split_line[2] == "ID":
                    id_snps.append(split_line[0] + ":" + split_line[1] + " " + split_line[3] + ">" + split_line[4])
                elif split_line[2] == "M" or split_line[2] == ".":
                    if "off-target" in split_line[6]:
                        off_target.append(split_line[0] + ":" + split_line[1] + " " + split_line[3] + ">" + split_line[4])
                    else:
                        monitoring_snps.append(split_line[0] + ":" + split_line[1] + " " + split_line[3] + ">" + split_line[4])
                else:
                    raise ValueError("Invalid variant type '" + split_line[2] + "'!")

    return set(monitoring_snps), set(id_snps), set(off_target)


def parse_blacklist_positions(file_path):
    blacklist_pos = []
    with open(file_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            if line.strip() == "":
                continue
            split_line = line.split("\t")
            blacklist_pos.append(split_line[0] + ":" + split_line[1])

    return set(blacklist_pos)


def get_alt_count(row, include_indel=False):
    if include_indel:
        alt = {"A", "C", "G", "T", "INS", "DEL"}
    else:
        alt = {"A", "C", "G", "T"}
    alt.remove(row["REF"])
    return row[list(alt)].sum()


def chr_compare(line1, line2):
    chr_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                 "chr20", "chr21", "chr22", "chrX", "chrY", "chrMT"]
    chr_pos1 = chr_order.index(line1["chr"])
    chr_pos2 = chr_order.index(line2["chr"])
    if chr_pos1 == chr_pos2:
        return int(line1["pos"]) - int(line2["pos"])
    else:
        return chr_pos1 - chr_pos2


def main():
    args = parse_args()

    # parse dedup files
    if not os.path.exists(args.folder):
        raise FileNotFoundError("umiVar folder '" + args.folder + "'not found!")
    dedup_tables = {}
    dedup_tables["1-fold"] = parse_table(os.path.join(args.folder, "dedup_DP1.tsv"))
    dedup_tables["2-fold"] = parse_table(os.path.join(args.folder, "dedup_DP2.tsv"))
    dedup_tables["3-fold"] = parse_table(os.path.join(args.folder, "dedup_DP3.tsv"))
    dedup_tables["4-fold"] = parse_table(os.path.join(args.folder, "dedup_DP4.tsv"))

    # parse cfDNA panel
    monitoring_snps, id_snps, off_target_snps = parse_panel(args.panel)

    # parse blacklist
    blacklist = set()
    if args.blacklist != "":
        blacklist = parse_blacklist_positions(args.blacklist)

    # output file buffer
    output = [["#MRD_log10", "MRD_pval", "SUM_DP", "SUM_ALT", "Mean_AF", "Stddev_AF", "Median_AF", "duplication", "error_rate", "BG_REF", "BG_other", "monitoring_count_pre_filter",
               "monitoring_count_post_filter", "monitoring_count_post_filter_with_counts", "proportion_Z_statistic", "minimal_detectable_af"]]

    monitoring_variant_counts = []
    for snp in monitoring_snps:
        snp_line = {}
        tmp1 = snp.split(' ')
        snp_line["chr"] = tmp1[0].split(':')[0]
        snp_line["pos"] = tmp1[0].split(':')[1]
        snp_line["REF"] = tmp1[1].split('>')[0]
        snp_line["ALT"] = tmp1[1].split('>')[1]

        # ignore indels
        if not args.keep_indels:
            if len(snp_line["ALT"]) != 1 or len(snp_line["REF"]) != 1:
                continue
        monitoring_variant_counts.append(snp_line)

    if not args.remove_off_target:
        for snp in off_target_snps:
            snp_line = {}
            tmp1 = snp.split(' ')
            snp_line["chr"] = tmp1[0].split(':')[0]
            snp_line["pos"] = tmp1[0].split(':')[1]
            snp_line["REF"] = tmp1[1].split('>')[0]
            snp_line["ALT"] = tmp1[1].split('>')[1]
            # ignore indels
            if len(snp_line["ALT"]) != 1 or len(snp_line["REF"]) != 1:
                continue
            monitoring_variant_counts.append(snp_line)

    # sort table
    monitoring_variant_counts = sorted(monitoring_variant_counts, key=cmp_to_key(chr_compare))

    # calculate MRD for 2-fold to 4-fold
    dup_levels = ["1-fold", "2-fold", "3-fold", "4-fold"]
    combined_tables = {}
    for dup_level in dup_levels:
        combined_tables[dup_level] = dedup_tables[dup_level][["A", "C", "G", "T", "DEL", "INS", "DP", "DP_HQ"]]
        included_levels = [dup_level]
        print(dup_level)
        for idx in range(dup_levels.index(dup_level) + 1, len(dup_levels)):
            combined_tables[dup_level] = combined_tables[dup_level].add(dedup_tables[dup_levels[idx]][["A", "C", "G", "T", "DEL", "INS", "DP", "DP_HQ"]], fill_value=0)
            included_levels.append(dup_levels[idx])
            print(" - " + dup_levels[idx])

        # get combined ref column
        ref_column = pd.concat([dedup_tables[dedup]["REF"] for dedup in included_levels])
        ref_column = ref_column[~ref_column.index.duplicated(keep='first')]
        combined_tables[dup_level]["REF"] = ref_column

        # use alias to simplify code
        combined = combined_tables[dup_level]

        # calculate ALT count and AF
        combined["ALT_count"] = combined.apply(get_alt_count, args=(args.keep_indels,), axis=1)
        combined["AF"] = combined["ALT_count"] / combined["DP_HQ"]

        # get monitoring variant count (pre-filter)
        n_mon_var_count_pre_filter = 0
        for snp in monitoring_variant_counts:
            idx_string = snp["chr"] + ":" + snp["pos"]
            if idx_string in combined.index:
                n_mon_var_count_pre_filter += 1

        # filter highAF variants
        if args.max_af < 1.0:
            combined = combined[combined["AF"] < args.max_af]

        # remove gonosomes
        if not args.keep_gonosomes:
            combined = combined[~combined.index.str.startswith("chrX:", "chrY:")]

        # filter table for monitoring SNPs
        monitoring_indices = list(set([snp.split(' ')[0] for snp in monitoring_snps]) & set(combined.index))
        # monitoring_table = combined.loc[monitoring_indices, ]

        # filter table for background data
        panel_indices = list(set(monitoring_indices + [snp.split(' ')[0] for snp in id_snps] + [snp.split(' ')[0] for snp in off_target_snps]) & set(combined.index))
        background_table = combined.loc[~combined.index.isin(panel_indices), ]  # remove panel variants
        background_table = background_table.loc[~background_table.index.isin(blacklist), ]  # remove blacklist variants

        # get monitoring counts for fisher test
        alt_counts_bg = int(background_table["ALT_count"].sum())
        ref_counts_bg = int(background_table["DP_HQ"].sum() - alt_counts_bg)

        # calculate p-value and AF for each monitoring variant:

        # first run calculate single variant values to determine distribution
        mon_afs = []
        monitoring_variant_counts_filtered = []
        for snp in monitoring_variant_counts:
            idx_string = snp["chr"] + ":" + snp["pos"]
            if idx_string in combined.index:
                snp["m_REF_" + dup_level] = int(combined.loc[idx_string, snp["REF"][0]])
                if len(snp["ALT"]) > 1:
                    snp["m_ALT_" + dup_level] = int(combined.loc[idx_string, "INS"])
                elif len(snp["REF"]) > 1:
                    snp["m_ALT_" + dup_level] = int(combined.loc[idx_string, "DEL"])
                else:
                    snp["m_ALT_" + dup_level] = int(combined.loc[idx_string, snp["ALT"]])
                snp["m_OTHER_" + dup_level] = max(0, int(get_alt_count(combined.loc[idx_string, :]) - snp["m_ALT_" + dup_level]))
                if snp["m_REF_" + dup_level] > 0:
                    snp["m_AF_" + dup_level] = snp["m_ALT_" + dup_level] / (snp["m_REF_" + dup_level] + snp["m_ALT_" + dup_level])
                else:
                    snp["m_AF_" + dup_level] = np.nan
                mon_afs.append(snp["m_AF_" + dup_level])
                # perform Fisher Exact Test for each variant
                contingency_table = np.array([[snp["m_ALT_" + dup_level], alt_counts_bg], [snp["m_REF_" + dup_level], ref_counts_bg]])
                res = scipy.stats.fisher_exact(contingency_table, alternative="greater")
                snp["p_value_" + dup_level] = res[1]

                monitoring_variant_counts_filtered.append(snp)
        # apply filter to list
        monitoring_variant_counts = monitoring_variant_counts_filtered

        # get statistics for outlier detection
        if len(mon_afs) > 0:
            mon_af_mean = np.mean(mon_afs)
            mon_af_median = np.median(mon_afs)
            mon_af_median_wo0 = np.median([af for af in mon_afs if af > 0.0])
            mon_af_std = np.std(mon_afs)
            mon_af_upper_limit = mon_af_mean + 3 * mon_af_std
        else:
            mon_af_mean = np.nan
            mon_af_median = np.nan
            mon_af_median_wo0 = np.nan
            mon_af_std = np.nan
            mon_af_upper_limit = np.nan

        # second run: remove outliers and sum up stats
        n_mon_var_count_post_filter = 0
        n_mon_var_count_post_filter_with_counts = 0
        alt_counts_mon = 0
        ref_counts_mon = 0
        alt_counts_mon_raw = 0
        ref_counts_mon_raw = 0
        n_monitoring = len(monitoring_variant_counts)
        n_background = background_table.shape[0]
        monitoring_variant_counts_filtered = []
        for snp in monitoring_variant_counts:
            ref_counts_mon_raw += snp["m_REF_" + dup_level]
            alt_counts_mon_raw += snp["m_ALT_" + dup_level]
            if snp["m_AF_" + dup_level] <= mon_af_upper_limit:
                snp["outlier_" + dup_level] = "0"
                # aggregate data
                ref_counts_mon += snp["m_REF_" + dup_level]
                alt_counts_mon += snp["m_ALT_" + dup_level]
                # count variants
                n_mon_var_count_post_filter += 1
                if snp["m_ALT_" + dup_level] > 0:
                    n_mon_var_count_post_filter_with_counts += 1
            else:
                snp["outlier_" + dup_level] = "1"

            # compute fold-change
            snp["fold-change_" + dup_level] = snp["m_AF_" + dup_level] / mon_af_mean - 1

            # perform Proportion Z Test for each variant
            # based on mean
            p_0 = mon_af_mean
            p = snp["m_AF_" + dup_level]
            n = len(mon_afs)
            snp["z_stats_mean-" + dup_level] = (p - p_0) / np.sqrt((p_0 * (1 - p_0)) / n)
            # based on median
            p_0 = mon_af_median
            if p_0 > 0.0:
                snp["z_stats_median-" + dup_level] = (p - p_0) / np.sqrt((p_0 * (1 - p_0)) / n)
            else:
                snp["z_stats_median-" + dup_level] = np.nan
            # based on median of >0
            p_0 = mon_af_median_wo0
            n = len([af for af in mon_afs if af > 0.0])
            if p_0 > 0.0 and n > 0:
                snp["z_stats_median_wo0-" + dup_level] = (p - p_0) / np.sqrt((p_0 * (1 - p_0)) / n)
            else:
                snp["z_stats_median_wo0-" + dup_level] = np.nan

        depth_mon = alt_counts_mon + ref_counts_mon
        depth_bg = alt_counts_bg + ref_counts_bg

        mon_af = np.nan
        if depth_mon > 0:
            mon_af = alt_counts_mon / depth_mon
        median_af = np.median([(snp["m_ALT_" + dup_level] / (snp["m_REF_" + dup_level] + snp["m_ALT_" + dup_level])) if (snp["m_REF_" + dup_level] + snp["m_ALT_" + dup_level]) > 0
                               else np.nan for snp in monitoring_variant_counts if ("m_ALT_" + dup_level) in snp.keys()])

        # print info
        print("Background AF: " + str(alt_counts_bg/depth_bg))
        if depth_mon > 0:
            print("Monitoring AF: " + str(alt_counts_mon/depth_mon))
        else:
            print("Monitoring AF: nan")

        # perform Proportion Z Test
        p_0 = alt_counts_bg/ref_counts_bg
        p = alt_counts_mon_raw/ref_counts_mon_raw
        n = n_monitoring + n_background

        z = (p-p_0)/np.sqrt((p_0*(1-p_0))/n)
        print("proportion Z statistic: " + str(z))


        # perform Fisher Exact Test
        contingency_table = np.array([[alt_counts_mon, alt_counts_bg], [ref_counts_mon, ref_counts_bg]])
        res = scipy.stats.fisher_exact(contingency_table, alternative="greater")
        print("fisher statistic: " + str(res[0]))
        print("fisher pvalue: " + str(res[1]))

        # dertermine minimal detectable AF:
        for i in range(alt_counts_bg):
            table = np.array([[i+1, alt_counts_bg], [ref_counts_mon, ref_counts_bg]])
            p_value = scipy.stats.fisher_exact(table, alternative="greater")[1]
            if p_value <= 0.05:
                break
        min_detectable_af = np.nan
        if p_value <= 0.05 and depth_mon > 0:
            min_detectable_af = i / depth_mon
        print("minimal detectable AF: " + str(min_detectable_af))

        # store output
        output.append(["{:.3f}".format(np.log10(max(res[1], 1e-20))), #MRD_log10
                       "{:.5f}".format(res[1]), #MRD_pval
                       str(ref_counts_mon), #SUM_DP
                       str(alt_counts_mon), #SUM_ALT
                       "{:.5f}".format(mon_af), #Mean_AF
                       "{:.5f}".format(mon_af_std), #Stddev_AF
                       "{:.5f}".format(median_af), #Median_AF
                       dup_level, #duplication
                       "{:.8f}".format(alt_counts_bg/(depth_bg)), #error_rate
                       str(ref_counts_bg), #BG_REF
                       str(alt_counts_bg), #BG_other
                       str(n_mon_var_count_pre_filter), #monitoring_count_pre_filter
                       str(n_mon_var_count_post_filter), #monitoring_count_post_filter
                       str(n_mon_var_count_post_filter_with_counts), #monitoring_count_post_filter_with_counts
                       "{:.5f}".format(z), #proportion_Z_statistic
                       "{:.8f}".format(min_detectable_af) #minimal_detectable_af
                       ])

    # write output files
    with open(args.output, 'w') as mrd_file:
        for line in output:
            mrd_file.write("\t".join(line))
            mrd_file.write("\n")

    # format p-value to precision of 8
    monitoring_variant_counts_df = pd.DataFrame(monitoring_variant_counts)
    monitoring_variant_counts_df["p_value_1-fold"] = monitoring_variant_counts_df["p_value_1-fold"].map(lambda x: '{0:.8}'.format(x))
    monitoring_variant_counts_df["p_value_2-fold"] = monitoring_variant_counts_df["p_value_2-fold"].map(lambda x: '{0:.8}'.format(x))
    monitoring_variant_counts_df["p_value_3-fold"] = monitoring_variant_counts_df["p_value_3-fold"].map(lambda x: '{0:.8}'.format(x))
    monitoring_variant_counts_df["p_value_4-fold"] = monitoring_variant_counts_df["p_value_4-fold"].map(lambda x: '{0:.8}'.format(x))

    monitoring_variant_counts_df.to_csv(args.variants, sep='\t', index=False)

    print("finished")

if __name__ == '__main__':
    main()


