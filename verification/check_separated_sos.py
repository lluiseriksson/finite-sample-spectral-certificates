"""Exact certificate replay; Python standard library only, no sampled data."""
from fractions import Fraction as F
from pathlib import Path
import json, copy, hashlib

ROOT = Path(__file__).resolve().parent

def add(p,q):
    r=p.copy()
    for k,v in q.items(): r[k]=r.get(k,F(0))+v
    return {k:v for k,v in r.items() if v}

def scale(p,c): return {k:v*c for k,v in p.items() if v*c}

def mul(p,q):
    r={}
    for (i,j),u in p.items():
        for (k,l),v in q.items():
            key=(i+k,j+l); r[key]=r.get(key,F(0))+u*v
    return {k:v for k,v in r.items() if v}

def parse(terms):
    p={}
    for i,j,v in terms:
        assert type(i) is int and type(j) is int and i>=0 and 0<=j<=2
        assert type(v) is str
        p=add(p,{(i,j):F(v)})
    return p

def verify_cert(cert):
    assert cert['variables']==['theta','t']
    assert [(r['a'],r['b']) for r in cert['pairs']]==[(0,3),(2,5),(0,5)]
    for row in cert['pairs']:
        q={}
        for sq in row['squares']:
            w=F(sq['weight']); assert w>0
            p=parse(sq['terms']); q=add(q,scale(mul(p,p),w))
        a,b=row['a'],row['b']
        assert mul({(1,0):F(1),(0,1):F(-1)},q)=={(b-a,a):F(1),(0,b):F(-1)}

def sos_quadratic(row,theta,L,v):
    total=F(0)
    for sq in row['squares']:
        co=[F(0)]*3
        for i,j,c in sq['terms']: co[j]+=F(c)*theta**i
        # Uses the complete block localizer, including its off-diagonal blocks.
        total+=F(sq['weight'])*sum(co[i]*co[j]*v[r]*v[s]*L[2*i+r][2*j+s]
                    for i in range(3) for j in range(3) for r in range(2) for s in range(2))
    return total

def localizer(K,theta):
    return [[theta*K[2*(i//2+j//2)+i%2][2*(i//2+j//2)+j%2]
             -K[2*(i//2+j//2+1)+i%2][2*(i//2+j//2+1)+j%2]
             for j in range(6)] for i in range(6)]

def block_q(K,k,v):
    return sum(v[i]*v[j]*K[2*k+i][2*k+j] for i in range(2) for j in range(2))

def diag_cov(vals):
    return [[vals[i//2] if i==j else F(0) for j in range(12)] for i in range(12)]

def checks(cert):
    verify_cert(cert)
    # Basis of ALL symmetric 12x12 matrices, including arbitrary time crosses.
    # Signed basis elements need not be covariance matrices: this checks linear identities.
    count=0
    for theta in [F(1,2),F(9,16),F(1)]:
        for i in range(12):
            for j in range(i,12):
                K=[[F(0) for _ in range(12)] for _ in range(12)]
                K[i][j]=K[j][i]=F(1)
                L=localizer(K,theta)
                for v in [(F(1),F(0)),(F(0),F(1)),(F(2),F(-3))]:
                    for row in cert['pairs']:
                        a,b=row['a'],row['b']
                        rhs=theta**(b-a)*block_q(K,a,v)-block_q(K,b,v)
                        assert sos_quadratic(row,theta,L,v)==rhs
                        count+=1
    # Strict algebraic enlargement: old three diagonal constraints all vanish,
    # but a separated-time constraint is negative. This is an alternative, NOT a null.
    theta=F(1,2); vals=[F(1),F(1,2),F(1),F(1,2),F(1),F(1,2)]
    assert all(x>0 for x in vals)
    assert all(theta*vals[2*i]-vals[2*i+1]==0 for i in range(3))
    detected=theta**3*vals[0]-vals[3]; assert detected==F(-3,8)
    # Six pair inequalities are still incomplete: PSD localizer can fail.
    vals2=list(map(F,['1','1/2','1/10','1/20','1/100','1/200']))
    pairs=[(0,1),(2,3),(4,5),(0,3),(2,5),(0,5)]
    means=[theta**(b-a)*vals2[a]-vals2[b] for a,b in pairs]
    assert min(means)>=0 and min(vals2)>0
    L=localizer(diag_cov(vals2),theta)
    w=[F(1),F(0),F(-1),F(0),F(0),F(0)]
    bad=sum(w[i]*L[i][j]*w[j] for i in range(6) for j in range(6))
    assert bad==F(-3,10)
    bad_certs=[]
    c=copy.deepcopy(cert); c['pairs'][0]['squares'][1]['weight']='1/2'; bad_certs.append(c)
    c=copy.deepcopy(cert); c['pairs'][2]['squares'][2]['weight']='-2/3'; bad_certs.append(c)
    c=copy.deepcopy(cert); c['pairs'][1]['a']=1; bad_certs.append(c)
    rejected=0
    for c in bad_certs:
        try: verify_cert(c)
        except AssertionError: rejected+=1
        else: raise AssertionError('Altered certificate accepted')
    return dict(symbolic_identities=3,positive_rational_squares=7,
       full_symmetric_covariance_basis_checks=count,altered_certificates_rejected=rejected,
       strict_enlargement_mean=str(detected),incomplete_family_quadratic_form=str(bad),
       incomplete_family_means=list(map(str,means)),
       new_data=False,solver=False,formalization=False,
       certificate_sha256=hashlib.sha256((ROOT/'separated_time_certificate.json').read_bytes()).hexdigest())

if __name__=='__main__':
    result=checks(json.loads((ROOT/'separated_time_certificate.json').read_text()))
    print(json.dumps(result,indent=2))
