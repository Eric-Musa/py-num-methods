import math


TESTS = [
    (lambda x: math.sin(x) + 2*math.cos(0.5*x) + 2, 0, 5, 1e-3, 4.29124525),
    (lambda x: 4*x**3 - 10*x**2 + 2*x + 3, 0, 2, 1e-3, 0.83918141),
]



def bisection1d(f, a, b, eps=1e-3, v=False):
    assert b > a, f'b ({b}) must be greater than a ({a}) to have a valid search interval'
    fa, fb = f(a), f(b)
    assert (check := fa*fb < 0), f'f({a})*f({b}) = {check} !< 0, cannot prove root on C[a, b] by IVT'
    k = 0
    while abs(b-a) > eps:
        m = (a + b) / 2
        fm = f(m)
        if v: print(f'k={k}: a={a}, b={b}, m={m}, fm={fm}')
        xstar = fm * fa
        if xstar < 0:
            b = m
            fb = fm
        else:
            a = m
            fa = fm
        k += 1
    m = (a + b) / 2
    fm = f(m)
    print(f'Root approximated after {k} steps: x*={fm}... [a, b]=[{a}, {b}] --> [f(a), f(b)]=[{fa}, {fb}]')
    return m, (a, b, fa, fb, k)


def k_min_bisection(a, b, eps):
    return int(math.log((b-a)/eps)/math.log(2))




if __name__ == '__main__':

    def test_bisection1d():
        for (f, a, b, eps, soln) in TESTS:
            try:
                result = bisection1d(f, a, b, eps)
                assert result[0] < soln < result[1]
                assert abs(result[1] - result[0]) < eps
                assert result[4] >= k_min_bisection(a, b, eps)
            except:
                print('errored out')
    

    test_bisection1d()