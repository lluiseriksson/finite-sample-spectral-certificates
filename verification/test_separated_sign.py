"""Exact algebra and interface regressions, not a coverage simulation."""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'verification'))
from separated_sign_test import CERTIFIED_PAIRS, REQUIRED_ASSUMPTIONS, evaluate, projection_pair
from check_separated_sos import checks


def opts(pair=(0,3)):
    return dict(time_pair=list(pair),direction=[1,0],theta='1/2',n=5,
                enabled=True,assumptions=REQUIRED_ASSUMPTIONS)

def row_at(k,x=1):
    row=[[0,0] for _ in range(12)]; row[k]=[x,x]; return row


class SeparatedSignTests(unittest.TestCase):
    def test_symbolic_sos_full_covariance_and_counterexamples(self):
        cert=json.loads((ROOT/'verification/separated_time_certificate.json').read_text())
        report=checks(cert)
        self.assertEqual(report['symbolic_identities'],3)
        self.assertEqual(report['full_symmetric_covariance_basis_checks'],2106)
        self.assertEqual(report['altered_certificates_rejected'],3)

    def test_actual_api_projections_match_full_covariance_basis(self):
        v=[F(2,3),F(-5,7)]; theta=F(9,16)
        for pair in CERTIFIED_PAIRS:
            a,b=projection_pair(pair,v)
            for i in range(12):
                for j in range(i,12):
                    K=[[F(0)]*12 for _ in range(12)]; K[i][j]=K[j][i]=1
                    quad=lambda w:sum(w[r]*K[r][s]*w[s] for r in range(12) for s in range(12))
                    block=lambda k:sum(v[r]*K[2*k+r][2*k+s]*v[s] for r in range(2) for s in range(2))
                    weight=theta**(pair[1]-pair[0])
                    self.assertEqual(weight*quad(a)-quad(b),weight*block(pair[0])-block(pair[1]))

    def test_all_pairs_use_correct_exponent(self):
        for pair in CERTIFIED_PAIRS:
            row=row_at(2*pair[0],1); row[2*pair[1]]=[F(1,2),F(1,2)]
            r=evaluate([row]*5,**opts(pair))
            self.assertEqual(r['p_exact'],'1/32')
            self.assertEqual(r['time_pair'],list(pair)); self.assertEqual(r['theta'],'1/2')
            self.assertEqual(F(r['quadratic_weight']),F(1,2)**(pair[1]-pair[0]))
            self.assertNotIn('degree',r)
            self.assertFalse(evaluate([row_at(2*pair[0])]*5,**opts(pair))['reject'])

    def test_channel_cancellation_and_unused_coordinates(self):
        row=row_at(6); row[7]=[2,2]
        self.assertEqual(evaluate([row]*5,**(opts()|{'direction':[2,-1]}))['certified_negative'],0)
        row[11]=['bad','bad']
        with self.assertRaises(ValueError): evaluate([row]*5,**opts())

    def test_ties_ambiguity_and_replacements(self):
        self.assertEqual(evaluate([row_at(6)]*4+[row_at(6,0)],**opts())['p_exact'],'3/16')
        r=evaluate([[[-1,1]]*12]*5,**opts()); self.assertEqual(r['ambiguous'],5)
        self.assertEqual(r['p_exact'],'1')
        self.assertEqual(evaluate([row_at(6)]*5,**(opts()|{'max_replaced':5}))['p_exact'],'1')

    def test_projected_interval_count_is_conservative(self):
        for pair in CERTIFIED_PAIRS:
            for xl,xh in [(-2,-1),(-1,1),(1,2)]:
                for yl,yh in [(-2,-1),(-1,1),(1,2)]:
                    row=row_at(0,0);row[2*pair[0]]=[xl,xh];row[2*pair[1]]=[yl,yh]
                    r=evaluate([row],**(opts(pair)|{'n':1}))
                    weight=F(1,2)**(pair[1]-pair[0])
                    for x,y in itertools.product([F(xl),F(xl+xh,2),F(xh)],
                                                 [F(yl),F(yl+yh,2),F(yh)]):
                        self.assertLessEqual(r['certified_negative'],int(weight*x*x-y*y<0))

    def test_exact_boundary_stays_in_n(self):
        row=row_at(0,1);row[6]=[F(1,8),F(1,8)]
        r=evaluate([row]*5,**(opts()|{'theta':'1/4'}))
        self.assertEqual(r['certified_negative'],0);self.assertEqual(r['p_exact'],'1')

    def test_single_pair_and_model_guards(self):
        for change in [dict(time_pair=[[0,3],[0,5]]),dict(time_pair=[0,1]),dict(time_pair=[3,0]),
                       dict(time_pair=[0,5,2]),dict(time_pair=[False,3]),dict(time_pair=[0.,3]),
                       dict(time_pair='03'),dict(enabled=False),dict(assumptions=[]),dict(n=4),
                       dict(theta=.5),dict(theta=0),dict(theta=2),dict(direction=[1]*6),
                       dict(direction=[0,0]),dict(max_replaced=6),dict(alpha=True)]:
            with self.subTest(change=change),self.assertRaises(ValueError):
                evaluate([row_at(6)]*5,**(opts()|change))
        for name,value in [('degree',0),('pairs',[[0,3],[0,5]]),('scan',True)]:
            with self.assertRaises(TypeError): evaluate([row_at(6)]*5,**(opts()|{name:value}))

    def test_cli_separate_opt_in_and_invalid_modes(self):
        payload=opts();payload.pop('enabled');payload['assumptions']=sorted(REQUIRED_ASSUMPTIONS)
        payload['rows']=[row_at(6)]*5
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'input.json';p.write_text(json.dumps(payload))
            run=lambda *flags:subprocess.run([sys.executable,str(ROOT/'separated_sign_test.py'),str(p),*flags],capture_output=True,text=True)
            self.assertEqual(run().returncode,2)
            r=run('--enable');self.assertEqual(r.returncode,0,r.stderr)
            self.assertEqual(json.loads(r.stdout)['time_pair'],[0,3])
            payload['time_pair']=[[0,3],[0,5]];p.write_text(json.dumps(payload))
            self.assertEqual(run('--enable').returncode,2)
            payload['time_pair']=[0,3];payload['theta']=.5;p.write_text(json.dumps(payload))
            self.assertEqual(run('--enable').returncode,2)


if __name__=='__main__': unittest.main()
