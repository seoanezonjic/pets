import os, sys, glob, gzip
from pets.genomic_features import Genomic_Feature

def load_variants(variant_folder):
  variants = {}
  for pattern in ['*.tab', '*.vcf', '*.vcf.gz']:
    for path in glob.glob(os.path.join(variant_folder, pattern)):
      profile_id, ext = os.path.basename(path).split(".", maxsplit=1)
      if ext == 'tab' or ext == 'txt':
        vars = load_tabular_vars(path)
      elif ext == 'vcf' or ext == 'vcf.gz':
        vars = load_vcf(path, ext)
      variants[profile_id] = Genomic_Feature(vars)
  return variants


def load_tabular_vars(path):
  vars = []
  with open(path) as f:
    for line in f:
      fields = line.strip().split("\t")
      chr = fields[0].replace('chr','')
      start = int(fields[1])
      vars.append([chr, start, start])
  return vars


def load_vcf(path, ext): # Some compressed files are fragmented internally. If so, VCFfile only reads first fragment
  vars = []             # Use zcat original.vcf.gz | gzip > new.vcf.gz to obtain a contigous file
  if ext == 'vcf.gz':
    with gzip.open(path, 'rb') as gz:
      while True:
        line = gz.readline()
        if not line: break
        decoded_line = line.decode("ascii")
        if decoded_line.startswith("#"): continue
        chr, start, *rest = decoded_line.strip().split("\t")
        chr = chr.replace("chr", "")
        start = int(start)
        vars.append([chr, start, start])
  elif ext == 'vcf':
    with open(path) as f:
      for line in f:
        if line.startswith("#"): continue
        chr, start, *rest = line.strip().split("\t")
        chr = chr.replace("chr", "")
        start = int(start)
        vars.append([chr, start, start])
  print(len(vars))
  return vars


def load_evidences(evidences_path, hpo):
  genomic_coordinates = {}
  coord_files = glob.glob(os.path.join(evidences_path, '*.coords'))
  for cd_f in coord_files:
    entity = os.path.basename(cd_f).replace('.coords','')
    coordinates = load_coordinates(cd_f)
    genomic_coordinates[entity] = coordinates

  evidences = {}
  evidence_files = glob.glob(os.path.join(evidences_path, '*_HP.txt'))
  for e_f in evidence_files:
    pair = os.path.basename(e_f,).replace('.txt','')
    profiles, id2label = load_evidence_profiles(e_f, hpo)
    evidences[pair] = {"prof": profiles, "id2lab": id2label}
  return evidences, genomic_coordinates


def load_coordinates(file_path):
  coordinates = {}
  header = True
  with open(file_path) as f:
    for line in f:
      fields = line.strip().split("\t")
      if header:
        header = False
      else:
        entity, chrm, strand, start, stop = fields
        if chrm == 'NA':
          sys.stderr.write(f"Warning: Record {fields.__repr__} is undefined")
          continue
        coordinates[entity] = [chrm, int(start), int(stop), strand]
  return coordinates


def load_evidence_profiles(file_path, hpo):
  profiles = {}
  id2label = {}
  #count = 0
  with open(file_path) as f:
    for line in f:
      id, label, profile = line.strip().split("\t")
      hpos = profile.split(',')
      hpos, rejected_hpos = hpo.check_ids(hpos)
      if hpos and len(hpos) > 0:
        hpos = hpo.clean_profile(hpos)
        if hpos and len(hpos) > 0: profiles[id] = hpos
        id2label[id] = label
  return profiles, id2label


def load_hpo_ci_values(information_coefficient_file):
  hpos_ci_values = {}
  with open(information_coefficient_file) as f:
    for line in f:
      hpo_code, ci = line.rstrip().split("\t")
      hpos_ci_values[hpo_code] = float(ci)
  return hpos_ci_values


def write_tabulated_data(data, file, header = None):
  with open(file, 'w') as f:
    if header != None: f.write("\t".join(header))
    for row in data:
      f.write("\t".join(map(lambda x: str(x), row)) + "\n")


def load_profiles(file_path, hpo):
  profiles = {}
  with open(file_path) as f:
    for line in f:
      id, profile = line.rstrip().split("\t")
      hpos = profile.split(',')
      hpos, rejected_hpos = hpo.check_ids(hpos)
      if len(hpos) > 0:
        hpos = hpo.clean_profile(hpos)
        if len(hpos) > 0 : profiles[id] = hpos
  return profiles