#!/usr/bin/env python3
"""rescale は 0 でない値を 0 にするか。（＝ ±MIN が塞いだ穴が、ここにも在るか）"""
import numpy as np
from sedenion_tensor_logic import (to_trits, from_trits, rescale, sedenion_mult,
                                   sedenion_add, tensor_general_fp, M)

print("=" * 78)
print("① rescale は何をするか")
print("=" * 78)
f = 6
for v in (3**7, 3**6, 3**5, 3**3, 1, -1, 0):
    w = to_trits(v, 18)
    out = from_trits(rescale(w, f))
    tag = "  **0 でないのに 0 になった**" if (v != 0 and out == 0) else ""
    print(f"  {v:>8} を {f} 桁落とす → {out:>6}{tag}")

print()
print("=" * 78)
print("② 実際の当てはめで、何回起きるか（tensor_general_fp、8 引き）")
print("=" * 78)
fake = live = 0
for seed in range(8):
    rng = np.random.default_rng(seed)
    f, W = 6, 18
    q3 = lambda v: int(round(float(v) * 3**f))
    G = [[rng.uniform(-1, 1, M) for _ in range(2)] for _ in range(2)]
    Gw = [[[to_trits(q3(c), W) for c in e] for e in row] for row in G]
    xs = [rng.uniform(-1, 1, M) for _ in range(4)]
    xw = [[to_trits(q3(c), W) for c in q] for q in xs]
    st = {'g9': 0, 'g27': 0}
    # rescale を見張る
    import sedenion_tensor_logic as S
    orig = S.rescale
    def watched(word, ff):
        global fake, live
        before = from_trits(word)
        out = orig(word, ff)
        after = from_trits(out)
        if before != 0:
            live += 1
            if after == 0: fake += 1
        return out
    S.rescale = watched
    S.tensor_general_fp(Gw, xw, 2, f, st)
    S.rescale = orig
print(f"  rescale が 0 でない値に呼ばれた回数 : **{live}**")
print(f"  そのうち 0 になった回数            : **{fake}**  （**{100*fake/max(live,1):.1f}%**）")

print()
print("=" * 78)
print("③ それは何を壊すか — 0 になった成分は、後で 0 除算されうる")
print("=" * 78)
print("""  この回路には割り算がありません（`no division anywhere` と書いてある）。
  だから **今は壊れません。**

  でも 0 の意味は、既に二つに割れています:
      **本物の 0**        （e1+e10 と e4-e15 の積。零因子。**厳密に 0**）
      **潰れて 0**        （rescale が低い桁を捨てた）
  この二つが **同じ 16 トリットのベクトル**になっていて、見分けがつきません。""")
