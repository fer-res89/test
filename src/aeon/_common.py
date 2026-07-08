"""Shared infrastructure used across aeon subcommands: the embedded
CODEML/MCMCTree control-file and WAG matrix templates, subprocess/ctl-file
helpers used by `process`, `pipe-all`, and `run-dirs`.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

DUMMY_CTL_EMBED = '''          seed = -1
       seqfile = OG0014449_edit.fas
      treefile = species_tree_t_d.tre
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       seqtype = 2  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '>3.16<3.224'

         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 20 1   * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1    * gamma prior for sigma^2     (for clock=2 or 3)

      finetune = 1: 0.1  0.1  0.1  0.1  0.1  0.1  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 2000
      sampfreq = 10
       nsample = 20000
'''

WAG_DAT_EMBED = '''0.551571
0.509848  0.635346
0.738998  0.147304  5.429420
1.027040  0.528191  0.265256  0.0302949
0.908598  3.035500  1.543640  0.616783  0.0988179
1.582850  0.439157  0.947198  6.174160  0.021352  5.469470
1.416720  0.584665  1.125560  0.865584  0.306674  0.330052  0.567717
0.316954  2.137150  3.956290  0.930676  0.248972  4.294110  0.570025  0.249410
0.193335  0.186979  0.554236  0.039437  0.170135  0.113917  0.127395  0.0304501 0.138190
0.397915  0.497671  0.131528  0.0848047 0.384287  0.869489  0.154263  0.0613037 0.499462  3.170970
0.906265  5.351420  3.012010  0.479855  0.0740339 3.894900  2.584430  0.373558  0.890432  0.323832  0.257555
0.893496  0.683162  0.198221  0.103754  0.390482  1.545260  0.315124  0.174100  0.404141  4.257460  4.854020  0.934276
0.210494  0.102711  0.0961621 0.0467304 0.398020  0.0999208 0.0811339 0.049931  0.679371  1.059470  2.115170  0.088836  1.190630
1.438550  0.679489  0.195081  0.423984  0.109404  0.933372  0.682355  0.243570  0.696198  0.0999288 0.415844  0.556896  0.171329  0.161444
3.370790  1.224190  3.974230  1.071760  1.407660  1.028870  0.704939  1.341820  0.740169  0.319440  0.344739  0.967130  0.493905  0.545931  1.613280
2.121110  0.554413  2.030060  0.374866  0.512984  0.857928  0.822765  0.225833  0.473307  1.458160  0.326622  1.386980  1.516120  0.171903  0.795384  4.378020
0.113133  1.163920  0.0719167 0.129767  0.717070  0.215737  0.156557  0.336983  0.262569  0.212483  0.665309  0.137505  0.515706  1.529640  0.139405  0.523742  0.110864
0.240735  0.381533  1.086000  0.325711  0.543833  0.227710  0.196303  0.103604  3.873440  0.420170  0.398618  0.133264  0.428437  6.454280  0.216046  0.786993  0.291148  2.485390
2.006010  0.251849  0.196246  0.152335  1.002140  0.301281  0.588731  0.187247  0.118358  7.821300  1.800340  0.305434  2.058450  0.649892  0.314887  0.232739  1.388230  0.365369  0.314730

0.0866279 0.043972  0.0390894 0.0570451 0.0193078 0.0367281 0.0580589 0.0832518 0.0244313 0.048466  0.086209  0.0620286  0.0195027  0.0384319  0.0457631  0.0695179  0.0610127  0.0143859  0.0352742  0.0708956
'''

# =========================
# Utilities
# =========================


def run_cmd(cmd, cwd=None):
    msg = ' '.join(cmd) + (f" (cwd={cwd})" if cwd else "")
    print(f"Running: {msg}")
    subprocess.run(cmd, cwd=cwd, check=True)

def generate_ctl(template_text: str, seqfile: str, treefile: str, dest_path: Path):
    txt = re.sub(r'(?m)^(\s*seqfile\s*=).*$', rf"\1 {seqfile}", template_text)
    txt = re.sub(r'(?m)^(\s*treefile\s*=).*$', rf"\1 {treefile}", txt)
    dest_path.write_text(txt)


def _write_wag(dest_dir: Path, wag_file: str | None):
    wag_path = dest_dir / "wag.dat"
    if wag_file:
        # Copy user-provided WAG into working dir for predictable path
        shutil.copy(wag_file, wag_path)
    else:
        wag_path.write_text(WAG_DAT_EMBED)
    return wag_path

def _load_ctl_template(ctl_template: str | None) -> str:
    if ctl_template:
        return Path(ctl_template).read_text()
    return DUMMY_CTL_EMBED

