from cyvcf2 import VCF, Writer
import pandas as pd
import json
from collections import Counter
from statistics import mean, median

class VCFManager:

    def __init__(self, vcf_path):
        self.vcf_path = vcf_path

        self.data = {
            "metadata": {},
            "variants": [],
            "summary": {}
        }

    # ==========================================
    # LOADING
    # ==========================================

    def load(self):
        """Load a single-sample VCF file."""
        vcf = VCF(self.vcf_path)

        if len(vcf.samples) != 1:
            raise ValueError("VCF must contain exactly one sample.")

        self.data["metadata"] = {
            "file": self.vcf_path,
            "sample": vcf.samples[0]
        }

        self.data["variants"] = []

        for variant in vcf:
            print("the variant is:", variant.FORMAT)
            rd = variant.format("RD") if "RD" in variant.FORMAT else None
            ad = variant.format("AD") if "AD" in variant.FORMAT else None
            dp = variant.format("DP") if "DP" in variant.FORMAT else None
            gq = variant.format("GQ") if "GQ" in variant.FORMAT else None

            print(f"Processing variant: {variant.CHROM}:{variant.POS} {variant.REF}>{variant.ALT} with RD = {rd}, AD={ad}, DP={dp}, GQ={gq}")
            for i, alt in enumerate(variant.ALT, start=1):

                ad_ref = int(ad[0][0]) if ad is not None else None
                ad_alt = int(ad[0][i]) if ad is not None else None

                if ad_ref is not None and ad_ref >= 0 and ad_alt is not None and ad_alt >= 0:
                    total = sum(int(x) for x in ad[0] if x >= 0)
                    vaf = ad_alt / total if total > 0 else None
                else:
                    vaf = None

                self.data["variants"].append({
                    "CHROM": variant.CHROM,
                    "POS": variant.POS,
                    "ID": variant.ID,
                    "REF": variant.REF,
                    "ALT": alt,
                    "QUAL": variant.QUAL,
                    "FILTER": variant.FILTER,
                    "GT": variant.genotypes[0][:-1],
                    "AD_REF": ad_ref,
                    "AD_ALT": ad_alt,
                    "DP": int(dp[0][0]) if dp is not None else None,
                    "GQ": int(gq[0][0]) if gq is not None else None,
                    "VAF": vaf
                })

        vcf.close()

    # ==========================================
    # GETTERS
    # ==========================================

    def get_variants(self):
        """Devuelve todas las variantes."""
        return self.data["variants"]

    def get_summary(self):
        """Devuelve las estadísticas generales."""
        return self.data["summary"]

    def get_metadata(self):
        """Devuelve los metadatos del VCF."""
        return self.data["metadata"]

    # ==========================================
    # STATISTICS
    # ==========================================

    def calculate_all_stats(self):
        """Calculate all VCF statistics."""
        self.calculate_basic_stats()
        self.calculate_quality_stats()
        self.calculate_vaf_stats()
        self.calculate_genotype_stats()
        self.calculate_variant_type_stats()

        return self.data["summary"]

    def calculate_basic_stats(self):
        """Calculate general VCF statistics."""
        variants = self.data["variants"]

        self.data["summary"]["basic"] = {
            "total_variants": len(variants),
            "variants_per_chromosome": dict(Counter(v["CHROM"] for v in variants)),
            "pass_variants": sum(v["FILTER"] is None for v in variants),
            "filtered_variants": sum(v["FILTER"] is not None for v in variants)
        }


    def calculate_quality_stats(self):
        """Calculate QUAL, DP and GQ statistics."""
        variants = self.data["variants"]
        stats = {}

        for field in ["QUAL", "DP", "GQ"]:
            values = [
                v[field] for v in variants
                if v[field] is not None and 0 <= v[field] < 2147483645
            ]

            stats[field] = {
                "mean": mean(values) if values else None,
                "median": median(values) if values else None,
                "min": min(values) if values else None,
                "max": max(values) if values else None
            }

        self.data["summary"]["quality"] = stats


    def calculate_vaf_stats(self):
        """Calculate VAF statistics."""
        vafs = [
            v["VAF"] for v in self.data["variants"]
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
        """Calculate genotype statistics."""
        genotypes = [
            v["GT"] for v in self.data["variants"]
        ]

        stats = Counter()

        for gt in genotypes:
            if any(allele < 0 for allele in gt):
                stats["missing"] += 1
            elif all(allele == 0 for allele in gt):
                stats["homozygous_ref"] += 1
            elif len(set(gt)) == 1:
                stats["homozygous_alt"] += 1
            else:
                stats["heterozygous"] += 1

        self.data["summary"]["genotypes"] = dict(stats)


    def calculate_variant_type_stats(self):
        """Calculate variant types and Ti/Tv ratio."""
        variants = self.data["variants"]

        stats = Counter()
        transitions = {"AG", "GA", "CT", "TC"}

        for v in variants:
            ref, alt = v["REF"], v["ALT"]

            if len(ref) == 1 and len(alt) == 1:
                stats["SNP"] += 1

                if ref + alt in transitions:
                    stats["transitions"] += 1
                elif ref != alt and ref in "ACGT" and alt in "ACGT":
                    stats["transversions"] += 1

            elif len(ref) < len(alt):
                stats["insertion"] += 1
            elif len(ref) > len(alt):
                stats["deletion"] += 1
            elif len(ref) == len(alt):
                stats["MNP"] += 1
            else:
                stats["other"] += 1

        ti = stats["transitions"]
        tv = stats["transversions"]

        stats["ti_tv_ratio"] = ti / tv if tv > 0 else None

        self.data["summary"]["variant_types"] = dict(stats)

    # ==========================================
    # FILTERING
    # ==========================================

    def apply_filters(self, filters):
        """Apply filters sequentially in the given order."""
        self.reset_filters()

        for filter_type, value in filters:
            method = getattr(self, f"filter_by_{filter_type}")

            if isinstance(value, dict):
                method(**value)
            else:
                method(value)

        return self.get_variants()

    def _apply_filter(self, condition):
        """Apply a filter to the current variants."""
        variants = (
            self.data["variants"]
            if self.filtered_variants is None
            else self.filtered_variants
        )

        self.filtered_variants = [v for v in variants if condition(v)]
        return self.filtered_variants


    def filter_by_chromosome(self, chromosomes):
        """Filter variants by chromosome."""
        if isinstance(chromosomes, str):
            chromosomes = [chromosomes]

        return self._apply_filter(lambda v: v["CHROM"] in chromosomes)


    def filter_by_position(self, chrom, start, end):
        """Filter variants by genomic position."""
        return self._apply_filter(
            lambda v: v["CHROM"] == chrom and start <= v["POS"] <= end
        )


    def filter_by_variant_type(self, variant_types):
        """Filter variants by type."""
        if isinstance(variant_types, str):
            variant_types = [variant_types]

        def get_type(v):
            ref, alt = v["REF"], v["ALT"]

            if len(ref) == len(alt) == 1:
                return "SNP"
            if len(ref) != len(alt):
                return "INDEL"
            if len(ref) == len(alt):
                return "MNP"

            return "OTHER"

        return self._apply_filter(
            lambda v: get_type(v) in variant_types
        )


    def filter_by_quality(self, min_qual=None, max_qual=None):
        """Filter variants by QUAL."""
        return self._apply_filter(
            lambda v: v["QUAL"] is not None
            and (min_qual is None or v["QUAL"] >= min_qual)
            and (max_qual is None or v["QUAL"] <= max_qual)
        )


    def filter_by_filter_status(self, statuses):
        """Filter variants by FILTER status."""
        if isinstance(statuses, str):
            statuses = [statuses]

        return self._apply_filter(
            lambda v: (v["FILTER"] or "PASS") in statuses
        )


    def filter_by_genotype(self, genotypes):
        """Filter variants by genotype."""
        if isinstance(genotypes, str):
            genotypes = [genotypes]

        return self._apply_filter(
            lambda v: "/".join(map(str, v["GT"])) in genotypes
        )


    def filter_by_depth(self, min_dp=None, max_dp=None):
        """Filter variants by DP."""
        return self._apply_filter(
            lambda v: v["DP"] is not None
            and (min_dp is None or v["DP"] >= min_dp)
            and (max_dp is None or v["DP"] <= max_dp)
        )


    def filter_by_genotype_quality(self, min_gq=None, max_gq=None):
        """Filter variants by GQ."""
        return self._apply_filter(
            lambda v: v["GQ"] is not None
            and (min_gq is None or v["GQ"] >= min_gq)
            and (max_gq is None or v["GQ"] <= max_gq)
        )


    def filter_by_vaf(self, min_vaf=None, max_vaf=None):
        """Filter variants by VAF."""
        return self._apply_filter(
            lambda v: v["VAF"] is not None
            and (min_vaf is None or v["VAF"] >= min_vaf)
            and (max_vaf is None or v["VAF"] <= max_vaf)
        )


    def filter_by_gene(self, genes):
        """Filter variants by gene."""
        if isinstance(genes, str):
            genes = [genes]

        return self._apply_filter(
            lambda v: v.get("GENE") in genes
        )


    def filter_by_consequence(self, consequences):
        """Filter variants by consequence."""
        if isinstance(consequences, str):
            consequences = [consequences]

        return self._apply_filter(
            lambda v: v.get("CONSEQUENCE") in consequences
        )


    def filter_by_impact(self, impacts):
        """Filter variants by impact."""
        if isinstance(impacts, str):
            impacts = [impacts]

        return self._apply_filter(
            lambda v: v.get("IMPACT") in impacts
        )


    def filter_by_info(self, field, value):
        """Filter variants by INFO attribute."""
        return self._apply_filter(
            lambda v: v.get("INFO", {}).get(field) == value
        )


    def filter_variants(self, **filters):
        """Apply multiple filters."""
        for name, params in filters.items():
            method = getattr(self, f"filter_by_{name}")

            if isinstance(params, dict):
                method(**params)
            else:
                method(params)

        return self.filtered_variants


    def reset_filters(self):
        """Reset all applied filters."""
        self.filtered_variants = None
        return self.data["variants"]
    
    # ==========================================
    # EXPORT
    # ==========================================

    def to_dataframe(self):
        """Convert current variants to a DataFrame."""
        return pd.DataFrame(self.get_variants())


    def to_csv(self, output_path):
        """Export current variants to CSV."""
        self.to_dataframe().to_csv(output_path, index=False)


    def to_json(self, output_path):
        """Export metadata, variants and summary to JSON."""
        data = {
            "metadata": self.data["metadata"],
            "variants": self.get_variants(),
            "summary": self.data["summary"]
        }

        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)


    def to_vcf(self, output_path):
        """Export current variants to a VCF file."""
        variants = self.get_variants()

        selected = {
            (v["CHROM"], v["POS"], v["REF"], v["ALT"])
            for v in variants
        }

        vcf = VCF(self.vcf_path)
        writer = Writer(output_path, vcf)

        for variant in vcf:
            if any(
                (variant.CHROM, variant.POS, variant.REF, alt) in selected
                for alt in variant.ALT
            ):
                writer.write_record(variant)

        writer.close()
        vcf.close()