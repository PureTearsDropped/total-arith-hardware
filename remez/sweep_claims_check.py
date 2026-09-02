#!/usr/bin/env python3
# ⚠️ AI-assisted; verify. / 生成AI使用・要検証（research ブランチ・未監査）
"""sweep_claims_check — 掃引結果 (sweep/f32_*.json) を **現在の** claims() で 再判定する。

  掃引は 各 (fn, mode) の max_err（全入力での 最悪 ulp）と max_cross（hi の二進位跨ぎ）を記録しているので、
  主張の式だけを直した場合は 再掃引せずに 「max_err ≤ 新主張」で 全入力の合否が決まる（sup の単調性）。
  片側性（ge-lies / le-lies）と sunk / sign / no-sat / exact-claim は 主張の式に依らないので 元の判定をそのまま使う。

  使い方: python sweep_claims_check.py [sweep/f32_exp.json ...]
"""
import sys, os, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from funcs_spec import Spec

TOL = 2.0 ** -20                                       # sweep_f32.check と同じ余裕


def main(paths):
    spec = Spec("f32"); cl = spec.claims()
    bad_total = 0
    for p in paths:
        d = json.load(open(p))
        print(f"[{os.path.basename(p)}] group={d['group']} quick={d['quick']} jobs={d['jobs']} {d['seconds']}s")
        for key, a in sorted(d["stats"].items()):
            fn, mode = key.split("/")
            c = cl[fn]
            lim = c["near_ulp"] if mode == "near" else (c["side_ulp_hi"] if mode == "hi" else c["side_ulp_lo"])
            ok_err = a["max_err"][0] <= lim + TOL
            ok_cross = (mode != "hi") or a["max_cross"][0] <= 2 * lim + TOL
            other = [v for v in d["confirmed_violations"].get(key, []) if not str(v[0]).startswith("ulp")]
            ok = ok_err and ok_cross and not other
            bad_total += (not ok)
            print(f"  {key:12s} 最悪ulp={a['max_err'][0]:.6f} 主張={lim:.6f} {'✓' if ok_err else '✗'}"
                  + (f"  跨ぎ={a['max_cross'][0]:.6f} ≤ {2*lim:.6f} {'✓' if ok_cross else '✗'}" if mode == "hi" else "")
                  + f"  片側余裕={a['min_margin'][0]}  非最近接={a['n_not_cr']}  他の違反={len(other)}")
    print("PASS" if bad_total == 0 else f"FAIL ({bad_total})")
    return bad_total == 0


if __name__ == "__main__":
    paths = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, "sweep", "f32_*.json")))
    paths = [p for p in paths if not p.endswith("_quick.json")]
    sys.exit(0 if main(paths) else 1)
