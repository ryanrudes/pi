#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct tuple {
    mpz_t p;
    mpz_t q;
    mpz_t r;
} tuple_t;

tuple_t tuple_init(mpz_t p, mpz_t q, mpz_t r) {
    tuple_t t;
    mpz_init_set(t.p, p);
    mpz_init_set(t.q, q);
    mpz_init_set(t.r, r);
    return t;
}

// Temporary variables for performing long computations
mpz_t temp, temp2;

// Constants
mpz_t ONE;

// Number of iterations
mpz_t n;

tuple_t binary_split(mpz_t a, mpz_t b) {
    gmp_printf("lo: %Zd, hi: %Zd\n", a, b);
    mpz_t Pab, Qab, Rab;

    mpz_add_ui(temp, a, 1);

    if (mpz_cmp(b, temp) == 0) // if (b == a + 1)
    {
        // Pab = -(6 * a - 5) * (2 * a - 1) * (6 * a - 1)
        //     = a((108 - 72a)a - 46) + 5
        // 
        // Multiplications: 3
        // Additions: 3
        mpz_mul_ui(temp, a, 72); // 72a
        mpz_ui_sub(temp, 108, temp); // 108 - 72a
        mpz_mul(temp, temp, a); // (108 - 72a)a
        mpz_sub_ui(temp, temp, 46); // (108 - 72a)a - 46
        mpz_mul(temp, a, temp); // a((108 - 72a)a - 46)
        mpz_add_ui(temp, temp, 5); // a((108 - 72a)a - 46) + 5
        mpz_init_set(Pab, temp);

        // Qab = c1 * a^3
        mpz_pow_ui(temp, a, 3); // a^3
        mpz_mul_ui(temp, temp, 10939058860032000); // c1 * a^3
        mpz_init_set(Qab, temp);

        // Rab = Pab * (c2 * a + c3)
        mpz_mul_ui(temp, a, 545140134); // c2 * a
        mpz_add_ui(temp, temp, 13591409); // c2 * a + c3
        mpz_mul(temp, Pab, temp); // Pab * (c2 * a + c3)
        mpz_init_set(Rab, temp);
    }
    else {
        mpz_t m;
        mpz_add(temp, a, b); // a + b
        mpz_fdiv_q_ui(temp, temp, 2); // (a + b) / 2
        mpz_init_set(m, temp);

        tuple_t PQRam = binary_split(a, m);
        tuple_t PQRmb = binary_split(m, b);

        mpz_mul(temp, PQRam.p, PQRmb.p);
        mpz_init_set(Pab, temp);

        mpz_mul(temp, PQRam.q, PQRmb.q);
        mpz_init_set(Qab, temp);

        mpz_mul(temp, PQRmb.q, PQRam.r);
        mpz_mul(temp2, PQRam.p, PQRmb.r);
        mpz_add(temp, temp, temp2);
        mpz_init_set(Rab, temp);
    }

    return tuple_init(Pab, Qab, Rab);
}

void chudnovsky(mpz_t n) {
    tuple_t PQR1n = binary_split(ONE, n);

    mpf_t tempf1, tempf2, Q, R;

    mpf_init(Q);
    mpf_init(R);

    mpf_set_z(Q, PQR1n.q);
    mpf_set_z(R, PQR1n.r);

    mpf_init(tempf1);
    mpf_init(tempf2);

    mpf_sqrt_ui(tempf1, 10005); // sqrt(c5)
    mpf_mul_ui(tempf1, tempf1, 426880); // c4 * sqrt(c5)
    mpf_mul(tempf1, tempf1, Q); // c4 * sqrt(c5) * Q

    // c6 * Q + R
    mpf_mul_ui(tempf2, Q, 13591409); // c6 * Q
    mpf_add(tempf2, tempf2, R); // c6 * Q + R

    mpf_div(tempf1, tempf1, tempf2); // (c4 * sqrt(c5) * Q) / (c6 * Q + R

    gmp_printf("π: %.10000Ff\n", tempf1);
}

int main() {
    mpf_set_default_prec(100000000);

    // Initialize constants
    mpz_init_set_si(ONE, 1);

    // Initialize temporary variables
    mpz_init(temp);
    mpz_init(temp2);

    // Initialize number of iterations
    mpz_init_set_si(n, 1000000);

    chudnovsky(n);

    return 0;
}