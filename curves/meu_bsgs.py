import math
from edwards_crypto.curves import edwards_add, Double_and_Add, scalar_multiplication


def ponto_inverso(P, p):
    x, y = P
    return ((-x) % p, y)


def meu_baby_step_giant_step(P, Q, n, a, d, p):
    m = math.ceil(math.sqrt(n))

    baby_steps = {}
    cur = (0, 1)

    for j in range(m):
        if cur not in baby_steps:
            baby_steps[cur] = j
        cur = edwards_add(cur, P, a, d, p)

    Pm = Double_and_Add(P, m, a, d, p)
    inv_Pm = ponto_inverso(Pm, p)

    Qi = Q
    for i in range(m):
        if Qi in baby_steps:
            j = baby_steps[Qi]
            k_candidate = i * m + j

            if scalar_multiplication(P, k_candidate, a, d, p) == Q:
                return k_candidate % n

            for t in range(n):
                if scalar_multiplication(P, t, a, d, p) == Q:
                    return t

            return k_candidate % n

        Qi = edwards_add(Qi, inv_Pm, a, d, p)

    return None
