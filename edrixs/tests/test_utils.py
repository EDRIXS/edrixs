import sympy
import edrixs


def test_CT_imp_bath():
    n, Delta, E_L, E_d, U_dd = sympy.symbols('n \\Delta, E_L, E_d, U_{dd}')

    Eq1 = sympy.Eq(10*E_L + n*E_d + n*(n - 1)*U_dd/2, rhs=0)
    Eq2 = sympy.Eq(9*E_L + (n+1)*E_d + (n + 1)*n*U_dd/2, rhs=Delta)
    Eq3 = sympy.Eq(8*E_L + (n+2)*E_d + (n + 1)*(n + 2)*U_dd/2, rhs=2*Delta + U_dd)

    # from IPython.display import display
    # display(Eq1, Eq2, Eq3)

    answers = sympy.solve([Eq1, Eq2, Eq3], [E_d, E_L])

    E_d_eq = sympy.Eq(E_d, rhs=answers[E_d])
    E_L_eq = sympy.Eq(E_L, rhs=answers[E_L])

    # from IPython.display import display
    # display(E_d_eq, E_L_eq)

    n_val = 8
    Delta_val = 3
    U_dd_val = 2

    E_d_val = E_d_eq.rhs.evalf(subs={n: n_val,
                                     Delta: Delta_val,
                                     U_dd: U_dd_val})

    E_L_val = E_L_eq.rhs.evalf(subs={n: n_val,
                                     Delta: Delta_val,
                                     U_dd: U_dd_val})

    E_d_cal, E_L_cal = edrixs.CT_imp_bath(U_dd_val, Delta_val, n_val)

    assert E_d_val == E_d_cal
    assert E_L_val == E_L_cal


def test_CT_imp_bath_core_hole():
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

    # from IPython.display import display
    # display(Eq1, Eq2, Eq3, Eq4, Eq5, Eq6)

    answer = sympy.solve([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6], [E_dc, E_Lc, E_p])

    E_dc_eq = sympy.Eq(E_dc, rhs=answer[E_dc])
    E_Lc_eq = sympy.Eq(E_dc, rhs=answer[E_Lc])
    E_p_eq = sympy.Eq(E_p, rhs=answer[E_p])

    # from IPython.display import display
    # display(E_dc_eq, E_Lc_eq, E_p_eq)

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
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8)
    slater_i = res['slater_i']
    assert slater_i[0][0] == 'F0_11'
    assert slater_i[0][1] == 0.0
    assert slater_i[1][0] == 'F2_11'
    assert slater_i[1][1] == 12.234
    assert slater_i[2][0] == 'F4_11'
    assert slater_i[2][1] == 7.598
