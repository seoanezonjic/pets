

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
        self.var2attrs = {}

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
                # print(self.patient2attributes[patient])
                # print("-----------------")

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
        #print("Affected column:", self.affected_col)

    def build_variant_matrix(self):
        self.patient_ids = list(self.vcf_data.keys())

        all_variants = set()
        for patient_vars in self.vcf_data.values():
            all_variants.update(patient_vars.keys())

        self.variant_ids = sorted(all_variants)

        V = np.full(
            (len(self.variant_ids), len(self.patient_ids)),
            fill_value=-1,
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
        # print(self.patient_ids)
        np.save("variant_matrix.npy", self.V)
        with open("variant_ids.txt", "w") as f:
            for var_id in self.variant_ids:
                chrom, pos, ref, alt = var_id.split("_")
                f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

    def build_affected_column(self):
        affected_col = np.full(len(self.patient_ids), 0, dtype=np.int8)

        for j, patient_id in enumerate(self.patient_ids):
            affected_col[j] = int(self.patient2attributes[patient_id]["affected"])

        return affected_col
    
    def filter_by_inheritance(self, inheritance_pattern, patient_id=None):
        if inheritance_pattern == "AD":
            self.filtered_variants = self.filter_autosomal_dominant()
        elif inheritance_pattern == "AR":
            self.filtered_variants = self.filter_autosomal_recessive()
        elif inheritance_pattern == "ARC":
            self.filtered_variants = self.filter_autosomal_recessive_compound_het(patient_id=patient_id)
        else:
            raise ValueError(f"Unknown inheritance pattern: {inheritance_pattern}")
    
    def filter_autosomal_dominant(self):
        A = self.affected_col.astype(bool)

        # Pass -1 to01 in affected and in unaffected
        V_aff = self.V[:, A]
        V_aff[V_aff == -1] = 0 # TODO FGC: This is a very important point, because we are 
        #considering that the missing values in the affected patients are 0,
        # but this could not be the case, because this could be a variant that is not 
        # called in the reference genome, but is present in the affected patients, so we should not filter this out.
        # When using missing, selection was going from 3871 to 10880!
        V_unaff = self.V[:, ~A]
        V_unaff[V_unaff == -1] = 0

        affected_ok = np.all(V_aff >= 1, axis=1) # I am supponsing in the model that we can have heterozygous
        # or homozygous variants in the affected patients, so we will consider both cases as valid for the AD model.
        # That is not the case for jannovar rools, but difference in minimum and we obtain more variants (and all variants from jannovar)
        #unaffected_ok = np.all(V_unaff == 0, axis=1)
        unaffected_ok = np.all(V_unaff == 0, axis=1)

        rho = affected_ok & unaffected_ok

        # print(len(self.variant_ids), "variants before filtering")
        # print("Number of candidate variants:", int(np.sum(rho)))

        return self.apply_variant_filter(
            rho
        )
    
    def filter_autosomal_recessive(self):
        A = self.affected_col.astype(bool)

        # Pass -1 to 1 in affected and to 0 in unaffected, to consider missing values as 0 in unaffected patients.
        V_aff = self.V[:, A]
        V_aff[V_aff == -1] = 2
        V_unaff = self.V[:, ~A]
        V_unaff[V_unaff == -1] = 0

        affected_ok = np.all(V_aff == 2, axis=1)
        unaffected_ok = np.all(V_unaff == 0, axis=1)

        rho = affected_ok & unaffected_ok

        # print(len(self.variant_ids), "variants before filtering")
        # print("Number of candidate variants:", int(np.sum(rho)))

        return self.apply_variant_filter(
            rho
        )

    def filter_autosomal_recessive_compound_het(self, patient_id=None):

        father_id = self.patient2attributes[patient_id]["father"]
        father_variant = self.V[:, self.patient_ids.index(father_id)]
        mother_id = self.patient2attributes[patient_id]["mother"]
        mother_variant = self.V[:, self.patient_ids.index(mother_id)]
        A = (father_variant >= 1) * (mother_variant <= 1)
        B = (mother_variant >= 1) * (father_variant <= 1)
        vPMv = np.dot(A[:, np.newaxis], B[np.newaxis, :]) 
        # vPMv = vPMv + np.dot(father_variant[:, np.newaxis] == -1, mother_variant[np.newaxis, :] >= 1)
        # vPMv = vPMv + np.dot(father_variant[:, np.newaxis] >= 1, mother_variant[np.newaxis, :] == -1)
        print("vPMv shape:", vPMv.shape)
        print(self.patient2attributes)


        A = self.affected_col.astype(bool)
        # Pass -1 to 1 in affected and to 0 in unaffected, to consider missing values as 0 in unaffected patients.
        V_aff = self.V[:, A]
        V_aff[V_aff == -1] = 1
        V_unaff = self.V[:, ~A]
        V_unaff[V_unaff == -1] = 0

        # pass from variant to gene
        self.build_var2attr_matrix()
        # get parents V info
        # We create a matrix of variant pairs where 1 if two variants are in the same gene and 0 otherwise.
        # load in sparse
        # This ca
        from scipy.sparse import csr_matrix
        csr_var2attrs_matrix = csr_matrix(self.var2attrs_matrix)
        vGv = np.dot(csr_var2attrs_matrix, csr_var2attrs_matrix.T)
        # Set diagonal to 0, because we are not interested in the same variant, but in different variants in the same gene.
        vGv.setdiag(0)
        print("vGv shape:", vGv.shape)
        # We create a matrix of variant pairs where 1 if one variant comes from a parent and  the other for the mother.
        # Now we create a pointwise multiplication of the two matrices to get a matrix of variant pairs where 1 if two variants are in the same gene and one comes from the father and the other from the mother.
        vPGMv = vGv * vPMv
        v_unaff = np.zeros((V_unaff.shape[0], V_unaff.shape[0]), dtype=np.int8)
        for i in range(V_unaff.shape[1]):
            v_unaff = V_unaff[:, i][:, np.newaxis] * V_unaff[:, i][np.newaxis, :] + v_unaff
            print(f"v_unaff shape after patient {i}: {v_unaff.shape}")
        
        v_aff = np.ones((V_aff.shape[0], V_aff.shape[0]), dtype=np.int8)
        for i in range(V_aff.shape[1]):
            v_aff = (V_aff[:, i][:, np.newaxis] * V_aff[:, i][np.newaxis, :])
            print(f"v_aff shape after patient {i}: {v_aff.shape}")
        vPGMv = vPGMv * (v_aff != 0) * (v_unaff == 0)
        rho = np.any(vPGMv, axis=1) * (self.V[:, self.patient_ids.index(patient_id)] != 0) #* (self.V[:, self.patient_ids.index(patient_id)] != -1)

        # Now we have a matrix of variants x attributes (genes) with 1 if the variant is in the gene and 0 otherwise.
        # Now we can filter the variants by gene, and keep only the genes that 
        # have at least 2 variants in the affected patients and 0 variants in the unaffected patients.
        # For this we can use the matrix multiplication of the variant matrix with the var2attr matrix, 
        # and then filter the genes that have at least 2 variants in the affected patients and 0 variants in the unaffected patients.
        # V_aff_gene = np.dot(V_aff.T, self.var2attrs_matrix)
        # V_aff_gene = V_aff_gene.T
        # np.save("V_aff_gene.npy", V_aff_gene)
        # V_unaff_gene = np.dot(V_unaff.T, self.var2attrs_matrix)
        # V_unaff_gene = V_unaff_gene.T
        # np.save("V_unaff_gene.npy", V_unaff_gene)
        # affected_ok_gene = np.all(V_aff_gene >= 2, axis=1)
        # unaffected_ok_gene = np.all(V_unaff_gene <= 1, axis=1)

        # # Now we select those variants that are in the genes that have at least 2 variants in the affected patients and 0 variants in the unaffected patients.
        # selected_genes = np.where(affected_ok_gene & unaffected_ok_gene)[0]
        # selected_variants = np.any(self.var2attrs_matrix[:, selected_genes], axis=1)

        # # Now rho
        # rho = selected_variants

        #rho = affected_ok_gene & unaffected_ok_gene

        print(len(self.variant_ids), "variants before filtering")
        print("Number of candidate variants:", int(np.sum(rho)))

        return self.apply_variant_filter(
            rho
        )


    @staticmethod
    def B(x):
        return (x > 0).astype(np.int8)

    def filter_variants_by_type(self, variant_type):
        if variant_type == "de_novo":
            self.filter_denovo()
        elif variant_type == "homozygous":
            self.filter_homo_lethality()
        elif variant_type == "compound_het":
            self.filter_compound_het()

    def filter_denovo(self):
        # Just filter de novo variants which are "surely" de novo.
        # Nowadays the process of denovo variant calling is not very reliable, so this filter must be applied just with caution.
        # And being sure that the variants are really de novo, and not just a false positive.
        # https://academic.oup.com/bib/article/26/6/bbaf543/8315883?login=true
        keep_variants = set()

        for child_id, attributes in self.patient2attributes.items():
            father_id = attributes["father"]
            mother_id = attributes["mother"]

            if father_id is None or mother_id is None:
                continue

            if (
                child_id not in self.patient_ids
                or father_id not in self.patient_ids
                or mother_id not in self.patient_ids
            ):
                continue

            child_idx = self.patient_ids.index(child_id)
            father_idx = self.patient_ids.index(father_id)
            mother_idx = self.patient_ids.index(mother_id)

            child_gt = self.V[:, child_idx]
            father_gt = self.V[:, father_idx]
            mother_gt = self.V[:, mother_idx]

            rho_child = (
                (father_gt == 0)
                & (mother_gt == 0)
                & (child_gt == 1)
            )

            for var_id, keep in zip(self.variant_ids, rho_child):
                if keep:
                    keep_variants.add(var_id)

        
        rho = np.array(
            [var_id in keep_variants for var_id in self.variant_ids],
            dtype=bool
        )

        # print(len(self.variant_ids), "variants before filtering")
        # print("Number of de novo candidate variants:", int(np.sum(rho)))

        self.apply_variant_filter(rho)

    def apply_variant_filter(self, rho):
        """
        Apply a boolean variant filter to self.V and self.variant_ids.

        Parameters
        ----------
        rho : np.ndarray
            Boolean array with one value per variant.
            True means keep the variant.
        """

        rho = np.asarray(rho, dtype=bool)

        if rho.shape[0] != len(self.variant_ids):
            raise ValueError(
                f"Filter length {rho.shape[0]} does not match "
                f"number of variants {len(self.variant_ids)}"
            )

        self.V = self.V[rho, :]
        self.variant_ids = [
            var_id
            for var_id, keep in zip(self.variant_ids, rho)
            if keep
        ]

    def filter_homo_lethality(self):
        homo_alt_present = np.any(self.V == 2, axis=1)
        rho = ~homo_alt_present

        # print(len(self.variant_ids), "variants before filtering")
        # print("Number of homozygous lethal candidate variants:", int(np.sum(rho)))

        return self.apply_variant_filter(
            rho,
            output_file="candidate_homo_lethal_variants.txt"
        )

    def filter_compound_het(self):
        # Implementation for filtering compound heterozygous variants
        pass

    def load_variant_attributes(self, path2var_attr):
        # TODO FGC: Ensure compatibility of ids between varcode and this system.
        var2attrs = {}
        with open(path2var_attr, 'r') as f:
            for idx, line in enumerate(f):
                attr_name2idx = {}
                parts = line.strip().split('\t')
                if idx == 0:
                    for idx_col, col in enumerate(parts):
                        if idx_col > 3:
                            attr_name2idx[col] = idx_col
                var_id = "_".join(parts[0:4])
                # print(var_id)
                attrs = {attr_name2idx.get(i, f"attr_{i}"): part.split(",") for i, part in enumerate(parts[4:])}
                var2attrs[var_id] = attrs
        self.var2attrs = var2attrs

    def build_var2attr_matrix(self):
        # generate columns for each attribute in var2attrs
        all_attr_values = set()
        for attrs in self.var2attrs.values():
            #print(list(attrs.values()))
            for element in list(attrs.values())[0]:
                all_attr_values.add(element)
        all_attr_values = sorted(all_attr_values)
        v2attr= np.zeros(
            (len(self.variant_ids), len(all_attr_values)),
            dtype=np.int8
        )
        attr_value2idx = {v: i for i, v in enumerate(all_attr_values)}
        count_vars_in_attrs = 0
        for idx, variant in enumerate(self.variant_ids):
            if variant in self.var2attrs:
                count_vars_in_attrs += 1
                # print("I am in variants", variant)
                attrs = self.var2attrs[variant]
                for attr_value in list(attrs.values())[0]:  # Assuming each attribute has a list of values
                    v2attr[idx, attr_value2idx[attr_value]] = 1
        print(f"Number of variants with attributes: {count_vars_in_attrs} out of {len(self.variant_ids)}")
        self.var2attrs_matrix = v2attr
        np.save("var2attr_matrix.npy", self.var2attrs_matrix)
        self.attr_ids = all_attr_values
        with open("attr_ids.txt", "w") as f:
            for attr_id in self.attr_ids:
                f.write(f"{attr_id}\n")
