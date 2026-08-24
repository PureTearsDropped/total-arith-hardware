#!/usr/bin/env python3
"""4 値の 3:2 圧縮器 — 表を **導出**する（推測しない）。

三値は 既存の compress3（18 ゲート・(p,n) レールに 全加算器 2 個）。
四値 `?`=(1,1) を 足したとき、`low + 2·high` が 真の 和の 集合を 含むような
最小の (low, high) を **64 通り 全列挙で 求める**。そこから ゲートの 形を 読む。
"""
import itertools
QS = {-1:(-1,), 0:(0,), 1:(1,), '?':(-1,0,1)}
NAMES = (-1, 0, 1, '?')

def repr_set(low, high):
    return frozenset(l + 2*h for l in QS[low] for h in QS[high])

def tightest(target):
    best = None
    for lo in NAMES:
        for hi in NAMES:
            s = repr_set(lo, hi)
            if target <= s:
                if best is None or len(s) < len(best[0]):
                    best = (s, lo, hi)
    return best

print("4 値 3:2 圧縮器 — 導出した表")
print(f"{'a':>4}{'b':>4}{'c':>4} | {'真の和':>18} | {'low':>4}{'high':>4} | {'余分':>5}")
rows = []
nq_stat = {}
for a, b, c in itertools.product(NAMES, repeat=3):
    t = frozenset(x+y+z for x in QS[a] for y in QS[b] for z in QS[c])
    best = tightest(t)
    assert best, f"表せない: {a},{b},{c} → {sorted(t)}"
    s, lo, hi = best
    nq = sum(1 for d in (a,b,c) if d == '?')
    rows.append((a,b,c,lo,hi,len(s)-len(t),nq))
    nq_stat.setdefault(nq, set()).add((lo,hi))
# 代表だけ 表示
seen = set()
for a,b,c,lo,hi,ex,nq in rows:
    if nq in seen and nq != 0: continue
    if nq == 0 and (a,b,c) not in ((0,0,0),(1,1,1),(1,1,-1)): continue
    seen.add(nq)
    t = frozenset(x+y+z for x in QS[a] for y in QS[b] for z in QS[c])
    print(f"{str(a):>4}{str(b):>4}{str(c):>4} | {str(sorted(t)):>18} | {str(lo):>4}{str(hi):>4} | {ex:>5}")
print()
print("? の 個数ごとに 出る (low,high) の 組")
for nq in sorted(nq_stat):
    print(f"  ? が {nq} 個: {sorted(nq_stat[nq], key=str)}")
print()
ex = [r[5] for r in rows]
print(f"余分な 値の 個数: 合計 {sum(ex)} / 64 通り中 余分ゼロ {sum(1 for e in ex if e==0)} 通り")
