import os, copy, tempfile, subprocess, time


def goto_dir(new_dir):
  if not os.path.isdir(new_dir):
    os.makedirs(new_dir)
  os.chdir(new_dir)


def insert_path(path, insert):
  if path.startswith('/'):
    result = path
  else:
    (head, tail) = os.path.split(path)
    result = os.path.join(insert, head, tail)
  return result


def temp_fname(suffix=''):
  fd, fname = tempfile.mkstemp(suffix, 'tmp-', '.')
  f = os.fdopen(fd, 'w')
  f.close()
  os.unlink(fname)
  return fname


def fname_variant(fname):
  root, ext = os.path.splitext(fname)
  i = 1
  new_fname = "%s-%d%s" % (root, i, ext)
  while os.path.isfile(new_fname):
    i += 1
    new_fname = "%s-%d%s" % (root, i, ext)
  return new_fname


def clean_fname(fname):
  try:
    os.remove(fname)
  except:
    pass


def run_with_output(cmd):
  p = subprocess.Popen(cmd, shell=True,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
  return p.stdout.read()


def is_same_dict_in_file(parms, fname):
  try:
    saved_parms = eval(open(fname).read())
    result = saved_parms == parms
  except:
    result = False
  return result


def replace_dict(t, d, max_len = 0):
  """Replace $token in string t using a dictionary

  By default, converts floats to %.3f format:
    - t: string whose tokens will be replace_dict
    - d: dictionary of tokens (string keys) and their replacements
         (values, may be of any type)
    - max_len: maximum length (in characters) of replacement values."""

  s = copy.deepcopy(t)
  for (k, v) in d.iteritems():
    token = "$%s" % k
    if max_len == 0:
      if type(v) is float:
        s = s.replace(token, "%.3f" % v)
      else:
        s = s.replace(token, str(v))
    else:
      if type(v) is float:
        s = s.replace(token, ("%.3f" % v).rjust(max_len))
      else:
        s = s.replace(token, str(v).rjust(max_len))
  return s


def write_dict(d, fname):
  """Outputs the dictionary to a file in repr format."""
  open(fname, "w").write(format_dict(d))


def format_dict(d):
  """Makes repr() comaptible string of a dictionary"""
  s = "{ "
  n = len(d)
  items = sorted(d.items())
  max_len = max([len(repr(k)) for k in d.keys()])
  for i in range(n):
    if i > 0:
      s += "  "
    key, value = items[i]
    s += "%s : %s" % (repr(key).ljust(max_len), repr(value))
    if i < n-1:
      s += ", \n"
  s += " }"
  return s


def words_in_file(fname):
  result = []
  for line in open(fname).readlines():
    result.extend(line.split())
  return result


def strip_lines(pdb_txt, tag_func):
  new_lines = []
  for line in pdb_txt.splitlines():
    if tag_func(line):
      continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def sort_file_by_line(fname, key_func):
  lines = open(fname, 'r').read().splitlines()
  pairs = [(key_func(l), l) for l in lines]
  new_lines = [l for (v, l) in sorted(pairs)]
  open(fname, 'w').write('\n'.join(new_lines))


def replace_txt_in_file(fname, tag, replacement):
  f = open(fname, 'r')
  txt = f.read().replace(tag, replacement)
  f.close()
  f = open(fname, 'w')
  f.write(txt)
  f.close()


def read_parameters(fname):
  class DataHolder: pass
  f = open(fname, 'r')
  result = DataHolder()
  result.__dict__ = eval(f.read())
  f.close()
  return result


class Timer:
  def __init__(self):
    self._elapsed = 0;
    self._start = time.time()

  def start(self):
    self._start = time.time()

  def stop(self):
    self._elapsed = time.time() - self._start

  def elapsed(self):
    if self._elapsed == 0:
      return time.time() - self._start
    else:
      return self._elapsed

  def str(self):
    elapsed_time = self.elapsed()
    min = elapsed_time / 60.0
    sec = elapsed_time % 60.0
    s = ""
    if min >= 1.0:
      s += "%.f:" % min
    if sec < 0.01:
      s += "%07.4fs" % sec
    else:
      s += "%05.2fs" % sec
    return s

  def __str__(self):
    return self.str()


def strip_hydrogens(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    if line.startswith("ATOM"):
      atom_type = line[12:16]
      if atom_type[0] == "H" or atom_type[1] == "H":
        continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def renumber_residues(pdb_txt):

  get_res_tag = lambda line: line[17:27]

  sorted_res_tags = []
  lines = pdb_txt.splitlines()
  for line in lines:
    if line.startswith('ATOM') or line.startswith('HETATM'):
      tag = get_res_tag(line)
      if tag not in sorted_res_tags:
        sorted_res_tags.append(tag)
  
  res_tag_to_new_resnum = {}
  for i, tag in enumerate(sorted_res_tags):
    res_tag_to_new_resnum[tag] = "%3d " % (i+1)

  new_lines = []
  for line in lines:
    new_line = line
    for start_tag in ['ATOM', 'ANISOU', 'HETATM']:
      if line.startswith(start_tag):
        tag = get_res_tag(line)
        resnum = res_tag_to_new_resnum[tag]
        new_line = line[:23] + resnum + line[27:]
    new_lines.append(new_line)
  return '\n'.join(new_lines)


def strip_other_nmr_models(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    new_lines.append(line)
    if line.startswith("ENDMDL"):
      break
  return '\n'.join(new_lines)


def strip_alternative_atoms(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    new_line = line
    if line.startswith('ATOM'):
      alt_loc = line[16]
      if not alt_loc in [' ']:
        if alt_loc in ['A', 'a']:
          new_line = line[:16] + ' ' + line[17:]
        else:
          continue
    new_lines.append(new_line)
  return '\n'.join(new_lines)


def clean_pdb(in_pdb, out_pdb):
  txt = open(in_pdb, 'r').read()
  txt = strip_other_nmr_models(txt)
  txt = strip_lines(txt, lambda l: l.startswith('HETATM'))
  txt = strip_lines(txt, lambda l: l.startswith('ANISOU'))
  txt = strip_lines(txt, lambda l: l.startswith('CONECT'))
  txt = strip_lines(txt, lambda l: l.startswith('MASTER'))
  txt = strip_lines(txt, lambda l: l.startswith("ATOM") and \
                                   l[17:20] in ['WAT', 'TIP', 'HOH'])
  txt = strip_alternative_atoms(txt)
  txt = strip_hydrogens(txt)
  txt = renumber_residues(txt)
  open(out_pdb, 'w').write(txt)
