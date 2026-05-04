import sympy
import edrixs


def test_CT_imp_bath():
    """CT_imp_bath solves for impurity (E_d) and bath (E_L) energies with no
    core hole, using two conditions:
      1. E(d^{n+1} L^9) - E(d^n L^10) = Delta  (charge-transfer energy)
      2. E(d^n L^10) = 0                         (reference state)
    The expected closed-form values are verified against the function output.
    """
    n = 8
    Delta = 3
    U_dd = 2

    E_d_val = (10*Delta - n*(19 + n)*U_dd/2)/(10 + n)
    E_L_val = n*((1 + n)*U_dd/2-Delta)/(10 + n)

    E_d_cal, E_L_cal = edrixs.CT_imp_bath(U_dd, Delta, n)

    assert E_d_val == E_d_cal
    assert E_L_val == E_L_cal


def test_CT_imp_bath_core_hole():
    """CT_imp_bath_core_hole solves for E_dc, E_Lc, E_p with a core hole
    present, using three conditions:
      1. E(d^n  L^10 p^6) = 0          (reference: full core, no ligand holes)
      2. E(d^{n+1} L^10 p^5) = 0       (reference: one core hole)
      3. E(d^{n+1} L^9  p^6) - E(d^n L^10 p^6) = Delta  (same Delta as no-core)

    The test constructs six independent transition-energy equations
    symbolically (all d^n/d^{n+1}/d^{n+2} transitions for np=5 and np=6),
    solves them as an overdetermined system, evaluates numerically, and
    compares against the function.

    With np=6 (full core):
      Eq1: E(d^n  L^10) = 0
      Eq2: E(d^{n+1} L^9)  = Delta
      Eq3: E(d^{n+2} L^8)  = 2*Delta + U_dd

    With np=5 (one core hole):
      Eq4: E(d^{n+1} L^10) = 0
      Eq5: E(d^{n+2} L^9)  = Delta + U_dd - U_pd
      Eq6: E(d^{n+3} L^8)  = 2*Delta + 3*U_dd - 2*U_pd
    """
    n, Delta, E_Lc, E_dc, E_p, U_dd, U_pd = (
        sympy.symbols('n \\Delta E_{Lc} E_{dc} E_p U_{dd} U_{pd}'))

    Eq1 = sympy.Eq(
        6*E_p + 10*E_Lc + n*E_dc + n*(n - 1)*U_dd/2 + 6*n*U_pd,
        rhs=0)
    Eq2 = sympy.Eq(
        6*E_p + 9*E_Lc + (n + 1)*E_dc + (n + 1)*n*U_dd/2 + 6*(n + 1)*U_pd,
        rhs=Delta)
    Eq3 = sympy.Eq(
        6*E_p + 8*E_Lc + (n + 2)*E_dc + (n + 1)*(n+2)*U_dd/2 + 6*(n+2)*U_pd,
        rhs=2*Delta+U_dd)
    Eq4 = sympy.Eq(
        5*E_p + 10*E_Lc + (n + 1)*E_dc + (n + 1)*n*U_dd/2 + 5*(n + 1)*U_pd,
        rhs=0)
    Eq5 = sympy.Eq(
        5*E_p + 9*E_Lc + (n+2)*E_dc + (n + 2)*(n + 1)*U_dd/2 + 5*(n + 2)*U_pd,
        rhs=Delta + U_dd - U_pd)
    Eq6 = sympy.Eq(
        5*E_p + 8*E_Lc + (n + 3)*E_dc + (n + 3)*(n + 2)*U_dd/2 + 5*(n + 3)*U_pd,
        rhs=2*Delta + 3*U_dd - 2*U_pd)

    answer = sympy.solve([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6], [E_dc, E_Lc, E_p])

    E_dc_eq = sympy.Eq(E_dc, rhs=answer[E_dc])
    E_Lc_eq = sympy.Eq(E_dc, rhs=answer[E_Lc])
    E_p_eq = sympy.Eq(E_p, rhs=answer[E_p])

    n_val = 8
    Delta_val = 3
    U_dd_val = 2
    U_pd_val = 1

    subs = {n: n_val,
            Delta: Delta_val,
            U_dd: U_dd_val,
            U_pd: U_pd_val}

    E_dc_val = E_dc_eq.rhs.evalf(subs=subs)
    E_Lc_val = E_Lc_eq.rhs.evalf(subs=subs)
    E_p_val = E_p_eq.rhs.evalf(subs=subs)

    E_dc_cal, E_Lc_cal, E_p_cal = edrixs.CT_imp_bath_core_hole(U_dd_val, U_pd_val, Delta_val, n_val)

    assert E_dc_val == E_dc_cal
    assert E_Lc_val == E_Lc_cal
    assert E_p_val == E_p_cal


def test_get_atom_data():
    """get_atom_data returns Hartree-Fock Slater integrals and SOC parameters
    from the edrixs atomic database. Spot-checks F0, F2, F4 for Ni d^8.
    """
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8)
    slater_i = res['slater_i']
    assert slater_i[0][0] == 'F0_11'
    assert slater_i[0][1] == 0.0
    assert slater_i[1][0] == 'F2_11'
    assert slater_i[1][1] == 12.234
    assert slater_i[2][0] == 'F4_11'
    assert slater_i[2][1] == 7.598
