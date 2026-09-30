import contextlib, io, sys
from pathlib import Path
sys.stdout.reconfigure(encoding="utf-8")
class Sink(io.StringIO):
    def reconfigure(self, **kwargs):
        return None
with contextlib.redirect_stdout(Sink()):
    import balance_test as bt

corrected_tail=["ZZ_glc_transport","ZZ_glc_EX","ZZ_glycogen_EX","ZZ_g1p_transport","ZZ_g1p_EX","ZZ_Tre_amy","ZZ_Tre_cga","ZZ_Tre_EX","ZZ_Tre_transport","ZZ_TreT1_ADP","ZZ_TreT2_UDP","ZZ_TreYZ"]
corrected=bt.base_ids[:-12]+corrected_tail
print("sheet\tscenario\tbefore_max_abs_Sv\tafter_max_abs_Sv\tafter_top_metabolite\tafter_top_residual")
for sheet in ("Table S3","Table S4"):
    before=bt.calc(sheet, bt.original_ids, False)
    after=bt.calc(sheet, corrected, False)
    for b,a in zip(before,after):
        print(f"{sheet}\t{a[0]}\t{b[1]:.16g}\t{a[1]:.16g}\t{a[2][0][1]}\t{a[2][0][2]:.16g}")
print("GLOBAL_AFTER_MAX", max(a[1] for sheet in ("Table S3","Table S4") for a in bt.calc(sheet, corrected, False)))
