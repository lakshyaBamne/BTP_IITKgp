from Schemes.LF import run_lf
from Schemes.LLF import run_llf

PROBLEMS = ["MCW", "SCW", "BLW", "TORO-1", "TORO-2", "TORO-3", "TORO-4", "LAX", "SDW", "SEW"]
BCS = ["FREE", "REFLECTIVE"]
ORDERES = ["1", "2"]

# for PRB in PROBLEMS:
#     run_lf(PRB, "1")

# run_llf("MCW", "1")
# run_llf("SCW", "1")
# run_llf("BLW", "1")
# run_llf("LAX", "1")
run_llf("TORO-2", "1")

