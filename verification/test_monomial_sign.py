"""Exact regression tests; fixtures are algebra, not simulated coverage."""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from monomial_sign_test import (REQUIRED_ASSUMPTIONS, binomial_upper_tail,
    evaluate, projection_pair, square_bounds)


def options():
    return dict(degree=1, direction=[1, 0], theta="1/2", n=5,
                enabled=True, assumptions=REQUIRED_ASSUMPTIONS)


def row_at(index, value=1):
    row = [[0, 0] for _ in range(12)]
    row[index] = [value, value]
    return row


class MonomialSignTests(unittest.TestCase):
    def test_exact_tail_by_exhaustive_enumeration(self):
        for n in range(1, 9):
            counts = [sum(bits) for bits in itertools.product([0, 1], repeat=n)]
            for k in range(n + 2):
                self.assertEqual(binomial_upper_tail(n, k), F(sum(c >= k for c in counts), 2**n))

    def test_thresholds(self):
        for n, k in [(100, 59), (500, 269), (1000, 527)]:
            self.assertLessEqual(binomial_upper_tail(n, k), F(1, 20))
            self.assertGreater(binomial_upper_tail(n, k - 1), F(1, 20))

    def test_all_degrees_use_the_expected_raw_blocks(self):
        for i in range(3):
            opts = options() | {"degree": i}
            self.assertEqual(evaluate([row_at(4*i+2)]*5, **opts)["p_exact"], "1/32")
            self.assertFalse(evaluate([row_at(4*i)]*5, **opts)["reject"])

    def test_nontrivial_channel_direction(self):
        row = row_at(6); row[7] = [2, 2]
        self.assertEqual(evaluate([row]*5, **(options() | {"direction": [2, -1]}))["certified_negative"], 0)

    def test_ties_and_ambiguity_stay_in_n(self):
        result = evaluate([row_at(6)]*4 + [row_at(6, 0)], **options())
        self.assertEqual(result["n"], 5)
        self.assertEqual(result["p_exact"], "3/16")
        self.assertFalse(result["reject"])
        result = evaluate([[[-1, 1]]*12]*5, **options())
        self.assertEqual(result["ambiguous"], 5)
        self.assertEqual(result["p_exact"], "1")

    def test_known_replacements(self):
        self.assertEqual(evaluate([row_at(6)]*5, **(options() | {"max_replaced": 5}))["p_exact"], "1")

    def test_interval_square_enclosure(self):
        for lo in range(-3, 4):
            for hi in range(lo, 4):
                a, b = square_bounds(F(lo), F(hi))
                for j in range(9):
                    x = F(lo) + F(j, 8)*(hi-lo)
                    self.assertLessEqual(a, x*x); self.assertLessEqual(x*x, b)

    def test_guards(self):
        for change in [dict(enabled=False), dict(assumptions=[]), dict(n=4), dict(n=True),
                       dict(degree=3), dict(direction=[1]*6), dict(direction=[0,0]),
                       dict(theta=.5), dict(theta=True), dict(max_replaced=6), dict(alpha=1)]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                evaluate([row_at(6)]*5, **(options() | change))
        for row in [[[0,0]]*18, [[1,0]]*12, [[0.,0.]]*12, [[False,0]]*12]:
            with self.assertRaises(ValueError): evaluate([row]*5, **options())

    def test_full_covariance_identity(self):
        # Every symmetric12x12 basis entry, all degrees, nontrivial channels.
        for degree in range(3):
            c = [F(2,3), F(-5,7)]; a,b = projection_pair(degree,c)
            for u in range(12):
                for v in range(u,12):
                    K=[[F(0)]*12 for _ in range(12)];K[u][v]=K[v][u]=1
                    quad=lambda w:sum(w[i]*K[i][j]*w[j] for i in range(12) for j in range(12))
                    block=lambda t:sum(c[i]*K[2*t+i][2*t+j]*c[j] for i in range(2) for j in range(2))
                    self.assertEqual(F(9,16)*quad(a)-quad(b),F(9,16)*block(2*degree)-block(2*degree+1))

    def test_polynomial_cross_block_counterexample(self):
        theta=F(9,16); K=[[F(0)]*12 for _ in range(12)]
        for k in range(6):K[2*k][2*k]=K[2*k+1][2*k+1]=theta**k
        K[2][6]=K[6][2]=theta**2/2
        for i in range(3):
            for j in range(3):
                self.assertEqual(theta*K[2*(i+j)][2*(i+j)]-K[2*(i+j+1)][2*(i+j+1)],0)
        self.assertGreater(K[2][2]*K[6][6]-K[2][6]**2,0)
        self.assertEqual(theta*(K[0][0]+K[4][4])-K[2][2]-K[6][6]-2*K[2][6],-theta**2)

    def test_cli_opt_in_and_float_rejection(self):
        script=Path(__file__).resolve().parents[1]/"monomial_sign_test.py"
        payload=options();payload.pop("enabled");payload["assumptions"]=sorted(REQUIRED_ASSUMPTIONS)
        payload["rows"]=[row_at(6)]*5
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"input.json";p.write_text(json.dumps(payload))
            run=lambda *flags:subprocess.run([sys.executable,str(script),str(p),*flags],capture_output=True,text=True)
            self.assertEqual(run().returncode,2)
            r=run("--enable");self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(json.loads(r.stdout)["p_exact"],"1/32")
            payload["theta"]=.5;p.write_text(json.dumps(payload))
            self.assertEqual(run("--enable").returncode,2)


if __name__ == "__main__":
    unittest.main()
