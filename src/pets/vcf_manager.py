from cyvcf2 import VCF, Writer
import pandas as pd
import json
from collections import Counter
from statistics import mean, median


import json
import math
import re

from collections import Counter
from statistics import mean, median

import pandas as pd
from cyvcf2 import VCF, Writer
import numpy as np
import pandas as pd


class VCFManager:
    """
    Manage normalized, biallelic, single-sample VCF files.

    Supported allele-depth formats:
        - GATK:     AD = REF,ALT
        - DRAGEN:   AD = REF,ALT (Number=R)
        - VarScan2: RD = REF, AD = ALT

    Provides:
        - VCF loading and FORMAT normalization
        - Variant statistics
        - Variant filtering
        - CSV, JSON, and VCF export
    """

    REPORT_METRICS = ("QUAL", "DP", "GQ", "VAF", "AD_REF", "AD_ALT")

    def __init__(self, vcf_path):
        self.vcf_path = str(vcf_path)

        self.data = {
            "metadata": {},
            "variants": [],
            "summary": {}
        }

        self.filtered_variants = None

    # ==========================================
    # FORMAT PARSING HELPERS
    # ==========================================

    @staticmethod
    def _valid_number(value):
        """
        Convert a numerical FORMAT value to a non-negative integer.

        Returns None for missing, invalid, or sentinel values.
        """
        if value is None:
            return None

        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError):
            return None

        if not math.isfinite(number):
            return None

        # Handle missing and vector-end sentinels from htslib.
        if number < 0 or number >= 2147483645:
            return None

        if not number.is_integer():
            return None

        return int(number)

    @staticmethod
    def _valid_float(value):
        """
        Convert a numerical value to a finite float.

        Missing and invalid values return None.
        """
        if value is None:
            return None

        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError):
            return None

        return number if math.isfinite(number) else None

    @classmethod
    def _format_values(cls, variant, field):
        """
        Extract FORMAT values from a single-sample VCF record.

        Examples:
            AD=80,20 -> [80, 20]
            AD=20    -> [20]
            DP=100   -> [100]

        Returns None when the field is absent.
        """
        if field not in variant.FORMAT:
            return None

        try:
            raw = variant.format(field)
        except (KeyError, ValueError):
            return None

        if raw is None or len(raw) == 0:
            return None

        values = raw[0]

        if hasattr(values, "tolist"):
            values = values.tolist()

        if not isinstance(values, (list, tuple)):
            values = [values]

        return [
            cls._valid_number(value)
            for value in values
    ]

    @classmethod
    def _scalar_format(cls, variant, field):
        """
        Extract a single integer FORMAT value.

        Returns None when the field is missing or invalid.
        """
        values = cls._format_values(variant, field)

        if values is None or len(values) != 1:
            return None

        return values[0]

    @staticmethod
    def _get_ad_number(vcf):
        """
        Read the Number attribute of FORMAT/AD from the VCF header.

        Expected values:
            R -> REF + ALT allele depths
            1 -> Single allele-depth value
        """
        match = re.search(
            r"##FORMAT=<ID=AD,Number=([^,>]+)",
            vcf.raw_header
        )

        return match.group(1) if match else None

    @classmethod
    def _parse_allelic_depths(cls, variant, ad_number):
        """
        Normalize allele depths to (AD_REF, AD_ALT).

        Supported representations:
            GATK/DRAGEN: AD=80,20
            VarScan2:    RD=80, AD=20

        The declared FORMAT/AD Number takes precedence
        over automatic detection.
        """
        ad = cls._format_values(variant, "AD")
        rd = cls._scalar_format(variant, "RD")

        if ad is None:
            return rd, None

        # Standard allele-depth representation.
        if ad_number == "R":
            if len(ad) == 2:
                return ad[0], ad[1]

            return None, None

        # VarScan2 allele-depth representation.
        if ad_number == "1":
            if len(ad) == 1 and rd is not None:
                return rd, ad[0]

            return rd, None

        # Fallback for files without a usable AD declaration.
        if ad_number is None:
            if len(ad) == 2:
                return ad[0], ad[1]

            if len(ad) == 1 and rd is not None:
                return rd, ad[0]

        # Unsupported or ambiguous representation.
        return rd, None

    @staticmethod
    def _extract_info(variant):
        """
        Extract INFO attributes into a Python dictionary.
        """
        result = {}

        for key, value in variant.INFO:
            if hasattr(value, "tolist"):
                value = value.tolist()

            if isinstance(value, tuple):
                value = list(value)

            result[key] = value

        return result

    # ==========================================
    # LOADING
    # ==========================================

    def load(self):
        """
        Load a normalized, biallelic, single-sample VCF.

        Allele depths are normalized regardless of caller.

        Internal representation:
            AD_REF: Reference-supporting reads
            AD_ALT: Alternative-supporting reads
            VAF:    Alternative allele fraction (0.0 to 1.0)
        """
        vcf = VCF(self.vcf_path)

        try:
            if len(vcf.samples) != 1:
                raise ValueError(
                    "VCF must contain exactly one sample."
                )

            ad_number = self._get_ad_number(vcf)

            self.data["metadata"] = {
                "file": self.vcf_path,
                "sample": vcf.samples[0],
                "ad_number": ad_number
            }

            self.data["variants"] = []
            self.data["summary"] = {}
            self.filtered_variants = None

            for variant in vcf:

                if (
                    len(variant.ALT) != 1
                    or variant.ALT[0] is None
                ):
                    raise ValueError(
                        "Non-biallelic variant found at "
                        f"{variant.CHROM}:{variant.POS}"
                    )

                # Skip non-variant and symbolic ALT records.
                alt = variant.ALT[0]

                if (
                    alt == "."
                    or alt == "*"
                    or alt.startswith("<")
                    or "[" in alt
                    or "]" in alt
                ):
                    raise ValueError(
                        "Unsupported ALT allele at "
                        f"{variant.CHROM}:{variant.POS}: {alt}"
                    )

                # Extract genotype information.
                gt_data = variant.genotypes[0]

                gt = [
                    int(allele)
                    for allele in gt_data[:-1]
                ]

                phased = bool(gt_data[-1])

                # Extract common FORMAT fields.
                dp = self._scalar_format(variant, "DP")
                gq = self._scalar_format(variant, "GQ")

                # Normalize allele depths.
                ad_ref, ad_alt = self._parse_allelic_depths(
                    variant,
                    ad_number
                )

                # Calculate VAF from informative allele depths.
                if ad_ref is not None and ad_alt is not None:
                    total_depth = ad_ref + ad_alt

                    vaf = (
                        ad_alt / total_depth
                        if total_depth > 0
                        else None
                    )
                else:
                    vaf = None

                info = self._extract_info(variant)

                record = {
                    "CHROM": variant.CHROM,
                    "POS": variant.POS,
                    "ID": variant.ID,
                    "REF": variant.REF,
                    "ALT": alt,
                    "QUAL": self._valid_float(variant.QUAL),
                    "FILTER": variant.FILTER,
                    "GT": gt,
                    "PHASED": phased,
                    "AD_REF": ad_ref,
                    "AD_ALT": ad_alt,
                    "DP": dp,
                    "GQ": gq,
                    "VAF": vaf,
                    "INFO": info,
                    "GENE": info.get("GENE"),
                    "CONSEQUENCE": info.get("CONSEQUENCE"),
                    "IMPACT": info.get("IMPACT")
                }

                self.data["variants"].append(record)

        finally:
            vcf.close()

        return self.data["variants"]

    # ==========================================
    # GETTERS
    # ==========================================

    def get_variants(self):
        """Return currently selected variants."""
        if self.filtered_variants is not None:
            return self.filtered_variants

        return self.data["variants"]

    def get_summary(self):
        """Return calculated statistics."""
        return self.data["summary"]

    def get_metadata(self):
        """Return VCF metadata."""
        return self.data["metadata"]

    # ==========================================
    # STATISTICS
    # ==========================================

    def calculate_all_stats(self):
        """
        Calculate statistics for currently selected variants.
        """
        self.data["summary"] = {}

        self.calculate_basic_stats()
        self.calculate_quality_stats()
        self.calculate_vaf_stats()
        self.calculate_genotype_stats()
        self.calculate_variant_type_stats()

        return self.data["summary"]

    def calculate_basic_stats(self):
        """Calculate general variant statistics."""
        variants = self.get_variants()

        self.data["summary"]["basic"] = {
            "total_variants": len(variants),
            "variants_per_chromosome": dict(
                Counter(v["CHROM"] for v in variants)
            ),
            "pass_variants": sum(
                v["FILTER"] in (None, "PASS")
                for v in variants
            ),
            "filtered_variants": sum(
                v["FILTER"] not in (None, "PASS", ".")
                for v in variants
            )
        }

    def calculate_quality_stats(self):
        """
        Calculate QUAL, DP, and GQ statistics.
        """
        variants = self.get_variants()
        stats = {}

        for field in ["QUAL", "DP", "GQ"]:

            values = [
                v[field]
                for v in variants
                if (
                    v[field] is not None
                    and isinstance(v[field], (int, float))
                    and math.isfinite(v[field])
                    and v[field] >= 0
                )
            ]

            stats[field] = {
                "mean": mean(values) if values else None,
                "median": median(values) if values else None,
                "min": min(values) if values else None,
                "max": max(values) if values else None
            }

        self.data["summary"]["quality"] = stats

    def calculate_vaf_stats(self):
        """Calculate variant allele fraction statistics."""
        vafs = [
            v["VAF"]
            for v in self.get_variants()
            if v["VAF"] is not None
        ]

        self.data["summary"]["vaf"] = {
            "mean": mean(vafs) if vafs else None,
            "median": median(vafs) if vafs else None,
            "min": min(vafs) if vafs else None,
            "max": max(vafs) if vafs else None,
            "count": len(vafs)
        }

    def calculate_genotype_stats(self):
        """Calculate genotype distribution statistics."""
        variants = self.get_variants()

        stats = Counter()

        for variant in variants:
            gt = variant["GT"]

            if not gt or any(allele < 0 for allele in gt):
                stats["missing"] += 1

            elif all(allele == 0 for allele in gt):
                stats["homozygous_ref"] += 1

            elif all(allele == 1 for allele in gt):
                stats["homozygous_alt"] += 1

            elif 0 in gt and 1 in gt:
                stats["heterozygous"] += 1

            else:
                stats["other"] += 1

        self.data["summary"]["genotypes"] = dict(stats)

    @staticmethod
    def _get_variant_type(variant):
        """Classify a small variant by REF and ALT."""
        ref = variant["REF"].upper()
        alt = variant["ALT"].upper()

        if len(ref) == 1 and len(alt) == 1:
            return "SNP" if ref != alt else "OTHER"

        if len(ref) < len(alt):
            return "INDEL"

        if len(ref) > len(alt):
            return "INDEL"

        if len(ref) == len(alt) and ref != alt:
            return "MNP"

        return "OTHER"

    def calculate_variant_type_stats(self):
        """
        Calculate variant types and transition/transversion ratio.
        """
        variants = self.get_variants()

        stats = Counter()

        transitions = {
            "AG", "GA", "CT", "TC"
        }

        for variant in variants:
            ref = variant["REF"].upper()
            alt = variant["ALT"].upper()

            variant_type = self._get_variant_type(variant)

            if variant_type == "SNP":
                stats["SNP"] += 1

                pair = ref + alt

                if pair in transitions:
                    stats["transitions"] += 1

                elif ref in "ACGT" and alt in "ACGT":
                    stats["transversions"] += 1

            elif variant_type == "INDEL":
                if len(ref) < len(alt):
                    stats["insertion"] += 1
                else:
                    stats["deletion"] += 1

            elif variant_type == "MNP":
                stats["MNP"] += 1

            else:
                stats["other"] += 1

        ti = stats["transitions"]
        tv = stats["transversions"]

        stats["ti_tv_ratio"] = (
            ti / tv if tv > 0 else None
        )

        self.data["summary"]["variant_types"] = dict(stats)

    # ==========================================
    # FILTERING
    # ==========================================

    def reset_filters(self):
        """Reset the active variant selection."""
        self.filtered_variants = None
        self.data["summary"] = {}

        return self.data["variants"]

    def _apply_filter(self, condition):
        """
        Apply a filter to the current variant selection.
        """
        variants = self.get_variants()

        self.filtered_variants = [
            variant
            for variant in variants
            if condition(variant)
        ]

        self.data["summary"] = {}

        return self.filtered_variants

    def apply_filters(self, filters):
        """
        Reset and apply an ordered sequence of filters.

        Example:
            [
                ("quality", {"min_qual": 30}),
                ("depth", {"min_dp": 20}),
                ("vaf", {"min_vaf": 0.1})
            ]
        """
        self.reset_filters()

        for filter_type, value in filters:
            method = getattr(
                self,
                f"filter_by_{filter_type}",
                None
            )

            if method is None:
                raise ValueError(
                    f"Unknown filter type: {filter_type}"
                )

            if isinstance(value, dict):
                method(**value)
            else:
                method(value)

        return self.get_variants()

    def filter_variants(self, **filters):
        """Apply multiple filters to the current selection."""
        for name, params in filters.items():
            method = getattr(
                self,
                f"filter_by_{name}",
                None
            )

            if method is None:
                raise ValueError(
                    f"Unknown filter type: {name}"
                )

            if isinstance(params, dict):
                method(**params)
            else:
                method(params)

        return self.get_variants()

    def filter_by_chromosome(self, chromosomes):
        """Filter by chromosome."""
        if isinstance(chromosomes, str):
            chromosomes = [chromosomes]

        return self._apply_filter(
            lambda v: v["CHROM"] in chromosomes
        )

    def filter_by_position(self, chrom, start, end):
        """Filter by genomic position."""
        return self._apply_filter(
            lambda v:
                v["CHROM"] == chrom
                and start <= v["POS"] <= end
        )

    def filter_by_variant_type(self, variant_types):
        """Filter by SNP, INDEL, MNP, or OTHER."""
        if isinstance(variant_types, str):
            variant_types = [variant_types]

        return self._apply_filter(
            lambda v:
                self._get_variant_type(v) in variant_types
        )

    def filter_by_quality(self, min_qual=None, max_qual=None):
        """Filter by QUAL."""
        return self._apply_filter(
            lambda v:
                v["QUAL"] is not None
                and (
                    min_qual is None
                    or v["QUAL"] >= min_qual
                )
                and (
                    max_qual is None
                    or v["QUAL"] <= max_qual
                )
        )

    def filter_by_filter_status(self, statuses):
        """Filter by VCF FILTER status."""
        if isinstance(statuses, str):
            statuses = [statuses]

        def normalize_status(status):
            if status is None:
                return "PASS"
            return status

        return self._apply_filter(
            lambda v:
                normalize_status(v["FILTER"]) in statuses
        )

    def filter_by_genotype(self, genotypes):
        """
        Filter by genotype.

        Supports:
            0/0, 0/1, 1/0, 1/1, 0|1, 1|0, etc.
        """
        if isinstance(genotypes, str):
            genotypes = [genotypes]

        def genotype_string(variant):
            separator = "|" if variant["PHASED"] else "/"

            return separator.join(
                str(allele) if allele >= 0 else "."
                for allele in variant["GT"]
            )

        return self._apply_filter(
            lambda v: genotype_string(v) in genotypes
        )

    def filter_by_depth(self, min_dp=None, max_dp=None):
        """Filter by read depth."""
        return self._apply_filter(
            lambda v:
                v["DP"] is not None
                and (
                    min_dp is None
                    or v["DP"] >= min_dp
                )
                and (
                    max_dp is None
                    or v["DP"] <= max_dp
                )
        )

    def filter_by_genotype_quality(
        self,
        min_gq=None,
        max_gq=None
    ):
        """Filter by genotype quality."""
        return self._apply_filter(
            lambda v:
                v["GQ"] is not None
                and (
                    min_gq is None
                    or v["GQ"] >= min_gq
                )
                and (
                    max_gq is None
                    or v["GQ"] <= max_gq
                )
        )

    def filter_by_vaf(self, min_vaf=None, max_vaf=None):
        """Filter by variant allele fraction (0 to 1)."""
        return self._apply_filter(
            lambda v:
                v["VAF"] is not None
                and (
                    min_vaf is None
                    or v["VAF"] >= min_vaf
                )
                and (
                    max_vaf is None
                    or v["VAF"] <= max_vaf
                )
        )

    def filter_by_gene(self, genes):
        """Filter by GENE annotation."""
        if isinstance(genes, str):
            genes = [genes]

        return self._apply_filter(
            lambda v: v.get("GENE") in genes
        )

    def filter_by_consequence(self, consequences):
        """Filter by CONSEQUENCE annotation."""
        if isinstance(consequences, str):
            consequences = [consequences]

        return self._apply_filter(
            lambda v:
                v.get("CONSEQUENCE") in consequences
        )

    def filter_by_impact(self, impacts):
        """Filter by IMPACT annotation."""
        if isinstance(impacts, str):
            impacts = [impacts]

        return self._apply_filter(
            lambda v: v.get("IMPACT") in impacts
        )

    def filter_by_info(self, field, value):
        """Filter by an INFO attribute."""
        return self._apply_filter(
            lambda v:
                v.get("INFO", {}).get(field) == value
        )

    # ==========================================
    # EXPORT
    # ==========================================

    def to_dataframe(self):
        """Convert selected variants to a pandas DataFrame."""
        return pd.DataFrame(self.get_variants())

    def to_csv(self, output_path):
        """Export selected variants to CSV."""
        self.to_dataframe().to_csv(
            output_path,
            index=False
        )

    def to_json(self, output_path):
        """
        Export metadata, selected variants, and summary to JSON.
        """
        data = {
            "metadata": self.get_metadata(),
            "variants": self.get_variants(),
            "summary": self.get_summary()
        }

        with open(output_path, "w", encoding="utf-8") as file:
            json.dump(
                data,
                file,
                indent=4,
                allow_nan=False
            )

    def to_vcf(self, output_path):
        """
        Export selected variants to a VCF file.

        Original FORMAT fields and the original header
        are preserved. This method does not convert
        caller-specific FORMAT fields.
        """
        variants = self.get_variants()

        selected = {
            (
                v["CHROM"],
                v["POS"],
                v["REF"],
                v["ALT"]
            )
            for v in variants
        }

        vcf = VCF(self.vcf_path)
        writer = None

        try:
            writer = Writer(str(output_path), vcf)

            for variant in vcf:
                if len(variant.ALT) != 1:
                    continue

                key = (
                    variant.CHROM,
                    variant.POS,
                    variant.REF,
                    variant.ALT[0]
                )

                if key in selected:
                    writer.write_record(variant)

        finally:
            if writer is not None:
                writer.close()

            vcf.close()

    @staticmethod
    def _chromosome_key(chrom):
        value = str(chrom)
        label = value[3:] if value.lower().startswith("chr") else value
        label = label.upper()
        if label.isdigit():
            return (0, int(label), "")
        return (1, {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(label, 1000), label)

    @staticmethod
    def _matrix(frame, column_keys):
        """Return header-first, rectangular data for report_html heatmaps."""
        rows = [["Chromosome"] + [str(key) for key in column_keys]]
        for chrom, series in frame.iterrows():
            rows.append([
                str(chrom),
                *[
                    float(series[key]) if pd.notna(series[key]) else None
                    for key in column_keys
                ],
            ])
        return rows

    def _prepare_frame(self):
        frame = self.to_dataframe().copy()
        required = {"CHROM", "POS", "REF", "ALT"}
        if not required.issubset(frame.columns):
            raise ValueError("VCF variant data must include CHROM, POS, REF and ALT")
        if frame.empty:
            raise ValueError("Cannot generate a report from an empty variant selection")

        frame["CHROM"] = frame["CHROM"].astype(str)
        frame["POS"] = pd.to_numeric(frame["POS"], errors="coerce")
        frame = frame.loc[frame["POS"].notna() & (frame["POS"] >= 1)].copy()
        if frame.empty:
            raise ValueError("No valid genomic positions were found")
        frame["POS"] = frame["POS"].astype("int64")

        for metric in self.REPORT_METRICS:
            if metric not in frame.columns:
                frame[metric] = np.nan
            frame[metric] = pd.to_numeric(frame[metric], errors="coerce")
            frame.loc[~np.isfinite(frame[metric]), metric] = np.nan
        return frame

    def _build_summary_data(self, frame):
        summary = self.calculate_all_stats()
        basic = summary.get("basic", {})
        result = [
            ["Statistic", "Value"],
            ["Sample", str(self.get_metadata().get("sample", "Unknown"))],
            ["Variants", int(len(frame))],
            ["Chromosomes with variants", int(frame["CHROM"].nunique())],
            ["PASS variants", basic.get("pass_variants")],
            ["Filtered variants", basic.get("filtered_variants")],
            ["Ti/Tv ratio", summary.get("variant_types", {}).get("ti_tv_ratio")],
        ]
        for metric in ("QUAL", "DP", "GQ"):
            stats = summary.get("quality", {}).get(metric, {})
            result += [[f"{metric} mean", stats.get("mean")],
                       [f"{metric} median", stats.get("median")]]
        result.append(["VAF median", summary.get("vaf", {}).get("median")])
        return result

    def _build_chromosome_data(self, frame, chromosomes):
        rows = [["Chromosome", "Variants", "Mean DP", "Mean GQ", "Mean VAF"]]
        for chrom in chromosomes:
            subset = frame.loc[frame["CHROM"] == chrom]

            def mean_or_none(metric):
                mean_value = subset[metric].mean()
                return round(float(mean_value), 4) if pd.notna(mean_value) else None

            rows.append([chrom, len(subset), mean_or_none("DP"),
                         mean_or_none("GQ"), mean_or_none("VAF")])
        return rows

    def _build_genomic_profiles(self, frame, chromosomes, window_bp):
        frame = frame.copy()
        frame["WINDOW"] = (frame["POS"] - 1) // window_bp
        windows = sorted(frame["WINDOW"].unique())
        if len(windows) > 1500:
            raise ValueError("More than 1500 occupied genomic windows; increase window_bp")

        counts = frame.groupby(["CHROM", "WINDOW"]).size().unstack("WINDOW")
        counts = counts.reindex(index=chromosomes, columns=windows).fillna(0)
        profiles = {"variant_counts": self._matrix(counts, windows)}
        for metric in self.REPORT_METRICS:
            medians = frame.groupby(["CHROM", "WINDOW"])[metric].median().unstack("WINDOW")
            medians = medians.reindex(index=chromosomes, columns=windows)
            profiles[metric] = self._matrix(medians, windows)
        return profiles, [
            f"{w * window_bp / 1e6:g}-{(w + 1) * window_bp / 1e6:g} Mb"
            for w in windows
        ]

    def _build_distributions(self, frame, chromosomes, histogram_bins):
        distributions = {}
        for metric in self.REPORT_METRICS:
            subset = frame.loc[frame[metric].notna(), ["CHROM", metric]]
            if subset.empty:
                continue
            values = subset[metric].to_numpy(dtype=float)
            minimum, maximum = float(values.min()), float(values.max())
            if minimum == maximum:
                minimum -= 0.5
                maximum += 0.5
            edges = np.linspace(minimum, maximum, histogram_bins + 1)
            labels = [f"{edges[i]:.4g}–{edges[i + 1]:.4g}"
                      for i in range(histogram_bins)]
            matrix = [["Chromosome"] + labels]
            for chrom in ["Genome"] + chromosomes:
                selected = (values if chrom == "Genome" else
                            subset.loc[subset["CHROM"] == chrom, metric].to_numpy(dtype=float))
                counts, _ = np.histogram(selected, bins=edges)
                probabilities = counts / counts.sum() if counts.sum() else counts.astype(float)
                matrix.append([chrom] + probabilities.tolist())
            distributions[metric] = matrix
        return distributions

    @staticmethod
    def _build_histograms(distributions):
        """Build CanvasXpress barplot tables from the exact heatmap bins.

        Genome and chromosome profiles are empirical relative frequencies.
        No independently recalculated bin boundaries are introduced.
        """
        histograms = {}
        for metric, matrix in distributions.items():
            labels = matrix[0][1:]
            names = [row[0] for row in matrix[1:]]
            if not labels or not names or names[0] != "Genome":
                continue

            # Column 0 contains categorical bin names (row names).
            # Each remaining column is a numerical series.
            genome = [["Bin", "Genome"]] + [
                [label, float(matrix[1][index + 1])]
                for index, label in enumerate(labels)
            ]
            by_chromosome = [["Bin"] + names[1:]] + [
                [label] + [float(row[index + 1]) for row in matrix[2:]]
                for index, label in enumerate(labels)
            ]
            histograms[metric] = {
                "genome": genome,
                "chromosomes": by_chromosome,
            }
        return histograms

    def build_report_data(self, window_bp=5_000_000, histogram_bins=24):
        """Return a presentation-ready dictionary from the active variant selection.

        This method does not depend on report_html or Mako.
        """
        if not isinstance(window_bp, int) or window_bp <= 0:
            raise ValueError("window_bp must be a positive integer")
        if not isinstance(histogram_bins, int) or histogram_bins < 2:
            raise ValueError("histogram_bins must be an integer >= 2")
        frame = self._prepare_frame()
        chromosomes = sorted(frame["CHROM"].unique(), key=self._chromosome_key)
        profiles, window_labels = self._build_genomic_profiles(frame, chromosomes, window_bp)
        distributions = self._build_distributions(frame, chromosomes, histogram_bins)
        return {
            "sample": str(self.get_metadata().get("sample", "Unknown")),
            "chromosomes": chromosomes,
            "summary": self._build_summary_data(frame),
            "chrom_table": self._build_chromosome_data(frame, chromosomes),
            "position": profiles,
            "window_labels": window_labels,
            "window_bp": window_bp,
            "distributions": distributions,
            "histograms": self._build_histograms(distributions),
            "metrics": list(self.REPORT_METRICS),
            "statistical_note": (
                "Empirical probability mass is the fraction of observed variants in a value bin. "
                "It is not a pathogenicity probability or a per-base mutation rate. "
                "Missing values are excluded separately for each metric."
            ),
        }

