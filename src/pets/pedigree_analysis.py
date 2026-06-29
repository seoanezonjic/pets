

from dataclasses import fields
import gzip
import numpy as np


class PedigreeAnalyzer:
    def __init__(self, options=None):
        self.options = options or {}

        self.patient2attributes = {}
        self.vcf_data = {}

        self.patient_ids = []
        self.variant_ids = []
        self.V = None  # V ∈ Z3^{n x m}

        self.filtered_variants = None

    def load_data(self):
        self.load_pedigree(self.options["pedigree"])
        self.load_vcfs(self.options["vcfs"])
        self.build_matrixes()

    def load_pedigree(self, pedigree_file):
        with open(pedigree_file) as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                fam_id, patient, father, mother, sex, phenotype = line.split("\t")

                self.patient2attributes[patient] = {
                    "fam_id": fam_id,
                    "father": None if father == "0" else father,
                    "mother": None if mother == "0" else mother,
                    "sex": int(sex),
                    "affected": phenotype == "2"
                }
                print(self.patient2attributes[patient])
                print("-----------------")

    def load_vcf_merged(self, path):
        self.vcf_data = {}

        opener = gzip.open if path.endswith(".gz") else open
        mode = "rt" if path.endswith(".gz") else "r"

        sample_ids = None

        with opener(path, mode) as f:
            for line in f:
                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    fields = line.strip().split("\t")
                    sample_ids = fields[9:]

                    for patient_id in sample_ids:
                        self.vcf_data[patient_id] = {}

                    continue

                if sample_ids is None:
                    continue

                parsed = self.parse_merged_vcf_line(line, sample_ids)

                if parsed is None:
                    continue

                var_id, patient2zygosity = parsed

                for patient_id, zygosity in patient2zygosity.items():
                    self.vcf_data[patient_id][var_id] = zygosity

        for patient_id in self.vcf_data:
            print(
                f"Loaded {len(self.vcf_data[patient_id])} variants for patient {patient_id}"
            )

        return self.vcf_data

    def parse_merged_vcf_line(self, line, sample_ids):
        if line.startswith("#"):
            return None

        fields = line.strip().split("\t")

        chrom = fields[0].replace("chr", "")
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        var_id = f"{chrom}_{pos}_{ref}_{alt}"

        format_fields = fields[8].split(":")
        sample_columns = fields[9:]

        patient2zygosity = {}

        for patient_id, sample_value in zip(sample_ids, sample_columns):
            if sample_value in {".", "./.", ".|."}:
                continue

            sample_fields = sample_value.split(":")
            fmt = dict(zip(format_fields, sample_fields))

            gt = fmt.get("GT")

            if gt is None or gt in {"./.", ".|.", "."}:
                continue

            gt = gt.replace("|", "/")

            if gt == "1/1":
                zygosity = 2
            elif gt in {"0/1", "1/0"}:
                zygosity = 1
            elif gt == "0/0":
                zygosity = 0
            else:
                continue

            patient2zygosity[patient_id] = zygosity

        return var_id, patient2zygosity

    def load_vcfs(self, vcfs):
        self.vcf_data = {}

        for patient_id, vcf_path in vcfs.items():
            self.vcf_data[patient_id] = self.load_vcf_genotypes(vcf_path)
            print(f"Loaded {len(self.vcf_data[patient_id])} variants for patient {patient_id}")
            #print(self.vcf_data[patient_id])  # Print first 10 variants for each patient

    def load_vcf_genotypes(self, path):
        variants = {}

        opener = gzip.open if path.endswith(".gz") else open
        mode = "rt" if path.endswith(".gz") else "r"

        with opener(path, mode) as f:
            for line in f:
                parsed = self.parse_vcf_line(line)

                if parsed is not None:
                    var_id, zygosity = parsed
                    variants[var_id] = zygosity

        print("the number of variants in the vcf file is", len(variants))

        return variants


    def normalize_variant(self, chrom, pos, ref, alt):
        pos = int(pos)

        # eliminar prefijo común manteniendo una base de anclaje
        while (
            len(ref) > 1
            and len(alt) > 1
            and ref[0] == alt[0]
        ):
            ref = ref[1:]
            alt = alt[1:]
            pos += 1

        # eliminar sufijo común manteniendo una base de anclaje
        while (
            len(ref) > 1
            and len(alt) > 1
            and ref[-1] == alt[-1]
        ):
            ref = ref[:-1]
            alt = alt[:-1]

        return chrom, pos, ref, alt

    def parse_vcf_line(self, line):
        if line.startswith("#"):
            return None

        fields = line.strip().split("\t")

        chrom = fields[0].replace("chr", "")
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        # chrom, pos, ref, alt = self.normalize_variant(
        #     chrom,
        #     pos,
        #     ref,
        #     alt
        # )

        format_fields = fields[8].split(":")
        sample_fields = fields[9].split(":")

        fmt = dict(zip(format_fields, sample_fields))
        gt = fmt.get("GT")

        if gt is None:
            print(f"Warning: GT field not found in VCF line: {line.strip()}")
            return None

        gt = gt.replace("|", "/")

        if gt == "1/1":
            zygosity = 2
        elif gt in {"0/1", "1/0"}:
            zygosity = 1
        else:
            print(f"Warning: Unrecognized genotype format: {line.strip()}")
            return None

        var_id = f"{chrom}_{pos}_{ref}_{alt}"

        return var_id, zygosity
    
    def load_vcf_genotypes_ref(self, path):
        variants = {}

        opener = gzip.open if path.endswith(".gz") else open
        mode = "rt" if path.endswith(".gz") else "r"

        with opener(path, mode) as f:
            for line in f:
                parsed = self.parse_vcf_line_ref(line)

                if parsed is not None:
                    var_id, zygosity = parsed
                    variants[var_id] = zygosity

        print("the number of variants in the vcf file is", len(variants))

        return variants

    def parse_vcf_line_ref(self, line):
        if line.startswith("#"):
            return None

        fields = line.strip().split("\t")

        chrom = fields[0].replace("chr", "")
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        format_fields = fields[8].split(":")
        sample_fields = fields[9].split(":")

        fmt = dict(zip(format_fields, sample_fields))
        gt = fmt.get("GT")

        if gt is None:
            print(f"Warning: GT field not found in VCF line: {line.strip()}")
            return None

        gt = gt.replace("|", "/")

        if gt == "1/1":
            zygosity = 2
        elif gt in {"0/1", "1/0"}:
            zygosity = 1
        else:
            zygosity = 0  # Reference genotype

        var_id = f"{chrom}_{pos}_{ref}_{alt}"

        return var_id, zygosity

    def build_matrixes(self):
        self.build_variant_matrix()
        self.affected_col = self.build_affected_column()
        print("Affected column:", self.affected_col)

    def build_variant_matrix(self):
        self.patient_ids = list(self.vcf_data.keys())

        all_variants = set()
        for patient_vars in self.vcf_data.values():
            all_variants.update(patient_vars.keys())

        self.variant_ids = sorted(all_variants)

        V = np.zeros(
            (len(self.variant_ids), len(self.patient_ids)),
            dtype=np.int8
        )

        var2idx = {v: i for i, v in enumerate(self.variant_ids)}
        pat2idx = {p: j for j, p in enumerate(self.patient_ids)}

        for patient_id, patient_vars in self.vcf_data.items():
            j = pat2idx[patient_id]

            for var_id, zygosity in patient_vars.items():
                i = var2idx[var_id]
                V[i, j] = zygosity

        self.V = V
        # print(self.variant_ids)
        print(self.patient_ids)
        np.save("variant_matrix.npy", self.V)
        with open("variant_ids.txt", "w") as f:
            for var_id in self.variant_ids:
                chrom, pos, ref, alt = var_id.split("_")
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

    def build_affected_column(self):
        affected_col = np.zeros(len(self.patient_ids), dtype=np.int8)

        for j, patient_id in enumerate(self.patient_ids):
            affected_col[j] = int(self.patient2attributes[patient_id]["affected"])

        return affected_col
    
    def filter_by_inheritance(self, inheritance_pattern):
        if inheritance_pattern == "AD":
            self.filtered_variants = self.filter_autosomal_dominant()
        elif inheritance_pattern == "AR":
            self.filtered_variants = self.filter_autosomal_recessive()
        else:
            raise ValueError(f"Unknown inheritance pattern: {inheritance_pattern}")
        
    # def filter_autosomal_dominant(self):
    #     A = self.affected_col.astype(np.int8)

    #     B0V = self.B(self.V)

    #     n_affected = A.sum()

    #     affected_score = B0V @ A
    #     unaffected_score = B0V @ (1 - A)
    #     rho = (affected_score == n_affected) & (unaffected_score == 0)

    #     print(len(self.variant_ids), "variants before filtering")

    #     candidate_variants = [
    #         var_id
    #         for var_id, keep in zip(self.variant_ids, rho)
    #         if keep
    #     ]

    #     print("Number of candidate variants", len(set(candidate_variants)))
    #     print("Candidate variants:", candidate_variants[:10])

    #     return candidate_variants
    
    def filter_autosomal_dominant(self):
        A = self.affected_col.astype(bool)

        V_aff = self.V[:, A]
        V_unaff = self.V[:, ~A]

        # todos los afectados son heterocigotos
        affected_ok = np.all(V_aff == 1, axis=1)

        # todos los sanos son referencia
        unaffected_ok = np.all(V_unaff == 0, axis=1)

        rho = affected_ok & unaffected_ok

        candidate_variants = [
            var_id
            for var_id, keep in zip(self.variant_ids, rho)
            if keep
        ]
        print(len(self.variant_ids), "variants before filtering")
        print("Number of candidate variants", len(set(candidate_variants)))
        print("Candidate variants:", candidate_variants[:10])

        with open("candidate_variants.txt", "w") as f:
            for var_id in candidate_variants:
                chrom, pos, ref, alt = var_id.split("_")
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

        return candidate_variants


    @staticmethod
    def B(x):
        return (x > 0).astype(np.int8)

