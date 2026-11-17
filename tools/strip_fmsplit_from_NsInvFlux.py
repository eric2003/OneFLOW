"""Remove FMSplit*/FPSplit5 definitions from NsInvFlux.cpp (keep calls)."""
import re, sys
from pathlib import Path
p = Path(sys.argv[1] if len(sys.argv) > 1 else "codes/ns/src/NsInvFlux.cpp")
t = p.read_text(encoding="utf-8", errors="replace")
t2, n = re.subn(
    r"\nReal FMSplit1\( const Real & mach, const Real & ipn \)\s*\{.*?\n\}\s*"
    r"\nReal FMSplit2\( const Real & mach, const Real & ipn \)\s*\{.*?\n\}\s*"
    r"\nReal FMSplit4\( const Real & mach, const Real & beta, const Real & ipn \)\s*\{.*?\n\}\s*"
    r"\nReal FPSplit5\( const Real & mach, const Real & alpha, const Real & ipn \)\s*\{.*?\n\}\s*",
    "\n",
    t,
    count=1,
    flags=re.S,
)
if n != 1:
    # try without leading newline / alternate spacing
    t2, n = re.subn(
        r"Real FMSplit1\( const Real & mach, const Real & ipn \)\s*\{.*?\n\}\s*"
        r"Real FMSplit2\( const Real & mach, const Real & ipn \)\s*\{.*?\n\}\s*"
        r"Real FMSplit4\( const Real & mach, const Real & beta, const Real & ipn \)\s*\{.*?\n\}\s*"
        r"Real FPSplit5\( const Real & mach, const Real & alpha, const Real & ipn \)\s*\{.*?\n\}\s*",
        "",
        t,
        count=1,
        flags=re.S,
    )
print("removed", n, "block(s) from", p)
if n == 0:
    sys.exit("No FMSplit definitions found - open NsInvFlux.cpp and delete Real FMSplit1..FPSplit5 bodies manually")
if 'NsInvFluxPolicy.h' not in t2:
    t2 = t2.replace('#include "NsInvFlux.h"\n', '#include "NsInvFlux.h"\n#include "NsInvFluxPolicy.h"\n')
p.write_text(t2, encoding="utf-8")
