

from dataclasses import fields
import gzip
import numpy as np
from cyvcf2 import VCF


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
        self.load_vcf_merged(self.options["merged_vcf"])
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

        vcf = VCF(path)
        sample_ids = vcf.samples

        for patient_id in sample_ids:
            self.vcf_data[patient_id] = {}

        for variant in vcf:
            chrom = variant.CHROM.replace("chr", "")
            pos = variant.POS
            ref = variant.REF

            # Si hay varios ALT, esto genera ALT1,ALT2...
            alt = ",".join(variant.ALT)

            var_id = f"{chrom}_{pos}_{ref}_{alt}"

            for patient_id, gt in zip(sample_ids, variant.genotypes):
                allele1, allele2, phased = gt

                # Missing: ./.
                # Here is very important to notice that is equivalent to say
                # this is 0, but this could not be the case for example if the variant is not present in the reference genome,
                #  so we will consider this as a missing value and we will not consider this variant for this patient
                # TODO FGC: See what to do in this cases so we select non called variants with a 
                # correct default that is the more tolerant (dont filter out variants that are not called in the reference genome)
                # Modify this EVEN in AD where we sohould condier . -> 0 for unaf but .-> 1 for aff 
                # And for AR something like the same.
                # No information is not information of "no variant" but "no information" and we should not consider this as a
                #  0 for the unaffected patients, because this could be a variant that is not called in the reference genome,
                #  but is present in the affected patients, so we should not filter this out.
                if allele1 == -1 or allele2 == -1:
                    continue

                if allele1 == 0 and allele2 == 0:
                    zygosity = 0

                elif {allele1, allele2} == {0, 1}:
                    zygosity = 1

                elif allele1 == 1 and allele2 == 1:
                    zygosity = 2

                else:
                    # Casos multialélicos: 1/2, 2/2, 0/2, etc.
                    continue

                self.vcf_data[patient_id][var_id] = zygosity

        for patient_id in self.vcf_data:
            print(
                f"Loaded {len(self.vcf_data[patient_id])} variants for patient {patient_id}"
            )

        return self.vcf_data

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
    
    def filter_autosomal_dominant(self):
        A = self.affected_col.astype(bool)

        V_aff = self.V[:, A]
        V_unaff = self.V[:, ~A]

        # todos los afectados son heterocigotos
        print("V_aff shape:", V_aff.shape)
        print("The affected patients are:", A)
        affected_ok = np.all(V_aff >= 1, axis=1) # This is different from jannovar documentation, where
        # They select just HET variants directly in the analysis. This is -> V_aff == 1 

        # todos los sanos son referencia
        print("V_unaff shape:", V_unaff.shape)
        print("The unaffected patients are:", A)
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
    
    def filter_denovo(self):
        # Just filter de novo variants which are "surely" de novo.
        # Nowadays the process of denovo variant calling is not very reliable, so this filter must be applied just with caution.
        # And being sure that the variants are really de novo, and not just a false positive.
        # https://academic.oup.com/bib/article/26/6/bbaf543/8315883?login=true

        # The rules are, 
        variants_in_paretns = set()
        for patient_id, attributes in self.patient2attributes.items():
            father_id = attributes["father"]
            mother_id = attributes["mother"]

            if father_id is None or mother_id is None:
                continue

            father_variants = self.vcf_data.get(father_id, {})
            mother_variants = self.vcf_data.get(mother_id, {})

            variants_in_paretns.update(father_variants.keys())
            variants_in_paretns.update(mother_variants.keys())

    def filter_homo_lethality(self):
        # This is just homo in the variant level and not the gene level.
        pass

