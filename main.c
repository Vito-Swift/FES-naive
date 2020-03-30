/* implementation of fast exhaustive search over GF(2) */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

/* Functional prototypes */
void read_challenge(const char *cha_file, uint8_t ***mqs, uint64_t *var_num, uint64_t *eq_num, uint64_t *xvar_num);
void fast_exhaustive_search(uint8_t **mqs, uint64_t var_num, uint64_t eq_num, uint64_t xvar_num, uint8_t *solution);
bool verify_sols(uint8_t **mqs, uint8_t *solution, uint64_t var_num, uint64_t eq_num, uint64_t xvar_num);
void alloc_sys(uint8_t ***mqs, uint64_t eq_num, uint64_t xvar_num);
void free_sys(uint8_t **mqs, uint64_t eq_num);

int main() {
    uint8_t **mqs = NULL;
    uint64_t var_num, eq_num, xvar_num;
    const char *cha_file = "cha.txt";
    read_challenge(cha_file, &mqs, &var_num, &eq_num, &xvar_num);

    uint8_t solution[var_num];
    fast_exhaustive_search(mqs, var_num, eq_num, xvar_num, solution);
    if (verify_sols(mqs, solution, var_num, eq_num, xvar_num)) {
        printf("\t\tsolution valid\n");
        printf("\t\tsolution: [");
        for (uint64_t i = 0; i < var_num - 1; i++)
            printf("%d, ", solution[i]);
        printf("%d]\n", solution[var_num - 1]);
    } else {
        printf("\t\tsolution invalid\n");
    }
    free_sys(mqs, eq_num);
}

/* Implementation */
void reduce_sys(uint8_t ***mqs, uint64_t eq_num, uint64_t var_num, uint64_t xvar_num) {
    uint64_t eq_idx, var_idx, i, sqr_term_idx;
    for (eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
        for (var_idx = 0; var_idx < var_num; ++var_idx) {
            for (i = 0, sqr_term_idx = 0; i < var_idx; i++) {
                sqr_term_idx += i + 2;
            }
            // add the coefficient of x_{var_idx}^2 to x_{var_idx}
            (*mqs)[eq_idx][xvar_num - 2 - (var_num - 1 - var_idx)] ^= (*mqs)[eq_idx][sqr_term_idx];
            (*mqs)[eq_idx][sqr_term_idx] = 0;
        }
    }
}

void alloc_sys(uint8_t ***mqs, uint64_t eq_num, uint64_t xvar_num) {
    *mqs = (uint8_t **) malloc(eq_num * sizeof(uint8_t *));

    for (uint64_t i = 0; i < eq_num; i++) {
        (*mqs)[i] = (uint8_t *) malloc(xvar_num * sizeof(uint8_t));
    }
}

void free_sys(uint8_t **mqs, uint64_t eq_num) {
    for (uint64_t i = 0; i < eq_num; i++) {
        free(mqs[i]);
    }
    free(mqs);
}

const size_t MAX_PRE_LEN = 256;
const char *CHA_GF_LINE = "Galois Field";
const char *CHA_VAR_LINE = "Number of variables";
const char *CHA_EQ_LINE = "Number of polynomials";
const char *CHA_SEED_LINE = "Seed";
const char *CHA_EQ_START = "*********";

#define cbinom2(n) \
    ( ((n) * ((n)+1)) / 2)

#define deg2midx2(var1_idx, var2_idx) \
    ((var2_idx) > (var1_idx) ? (cbinom2(var2_idx) + (var1_idx)) : (cbinom2(var1_idx) + (var2_idx)))

#define deg2midx1(vnum, var1_idx) \
    (cbinom2(vnum) + (var1_idx))

bool check_prefix(const char *pre, const char *str) {
    return !strncmp(pre, str, strnlen(pre, MAX_PRE_LEN));
}

bool parse_cha_header(const char *str, uint64_t *var_num, uint64_t *eq_num) {
    if (check_prefix(CHA_EQ_START, str)) {
        printf("\t\treading equations...\n");
        return false;
    }
    uint64_t seed;

    if (check_prefix(CHA_VAR_LINE, str)) {
        if (1 != sscanf(str, "%*s %*s %*s %*s : %" PRIu64, var_num)) {
            printf("[!] cannot parse number of unknowns: %s\n", str);
            exit(1);
        }

        printf("\t\tnumber of variables: %d\n", *var_num);

    } else if (check_prefix(CHA_EQ_LINE, str)) {
        if (1 != sscanf(str, "%*s %*s %*s %*s : %" PRIu64, eq_num)) {
            printf("[!] cannot parse number of equations: %s\n", str);
            exit(1);
        }

        printf("\t\tnumber of equations: %" PRIu64 "\n", *eq_num);

    } else if (check_prefix(CHA_SEED_LINE, str)) {
        if (1 != sscanf(str, "%*s : %d", &seed)) {
            printf("[!] unable to seed: %s\n", str);
            exit(1);
        }

        printf("\t\tseed: %d\n", seed);

    } else if (check_prefix(CHA_GF_LINE, str)) {
        int prime = 0;
        if ((1 != sscanf(str, "%*s %*s : GF(%d)", &prime)) || prime != 2) {
            printf("[!] unable to process GF(%d)\n", prime);
            exit(1);
        }

        printf("\t\tfield: GF(%d)\n", prime);
    }

    return true;
}

void parse_cha_eqs(uint8_t ***mqs, char *str, const uint64_t eq_idx, const uint64_t xvar_num) {
    char *ptr = NULL;

    uint64_t i = 0;
    ptr = strtok(str, " ;");
    while (NULL != ptr) {
        (*mqs)[eq_idx][i++] = atoi(ptr);
        ptr = strtok(NULL, " ;\n");
    }

    assert(i == xvar_num);
}

void read_challenge(const char *cha_file, uint8_t ***mqs, uint64_t *var_num, uint64_t *eq_num, uint64_t *xvar_num) {
    FILE *fp = fopen(cha_file, "r");

    const size_t buf_size = 0x1 << 20; // 1MB per line
    char *buf = (char *) malloc(buf_size * sizeof(char));
    bool before_eq = true;
    uint64_t eq_idx = 0;
    while (NULL != fgets(buf, buf_size, fp)) {
        if (before_eq) {
            before_eq = parse_cha_header(buf, var_num, eq_num);

            if (!before_eq) {
                *xvar_num = (*var_num) * ((*var_num) + 1) / 2 + (*var_num) + 1;
                alloc_sys(mqs, *eq_num, *xvar_num);
            }
        } else {
            parse_cha_eqs(mqs, buf, eq_idx++, *xvar_num);
        }
    }

    reduce_sys(mqs, *eq_num, *var_num, *xvar_num);
}

void diff_eq(uint8_t *func, uint64_t term_num, uint64_t var_num, uint64_t idx,
             bool *result) {
    assert(idx < var_num);

    // set all terms to zero
    memset(result, 0x0, sizeof(bool) * (var_num + 1));

    // ignore constant term in func, it's gone

    // diff x_{i+1}, which becomes the constant term of the result
    result[var_num] = func[term_num - 1 - (var_num - idx)];

    // diff x1x2, x2x3, ...
    uint64_t cursor = term_num - 1 - var_num - 1;  // start with x_{var_num}^2
    uint64_t i, j;
    uint64_t bound = (0 == idx) ? 1 : idx;
    for (i = var_num - 1; i >= bound; --i) {  // check x_?x_{i+1}, and skip x1^2
        if (i == idx) {
            for (j = 1; j <= i; ++j) {
                result[i - j] ^= func[cursor - j];
            }
        } else {
            result[i] ^= func[cursor - (i - idx)];
        }
        cursor -= i + 1;  // move to x_{i-1}^2
    }
//    uint64_t i;
//    for (i = 0; i < idx; i++) {
//        result[i] = func[deg2midx2(i, idx)];
//    }
}

void find_partial_derivs(uint8_t **mqs, bool *derivs, uint64_t eq_num, uint64_t term_num, uint64_t var_num) {
    uint64_t eq_idx, var_idx;
    for (eq_idx = 0; eq_idx < eq_num; eq_idx++) {
        for (var_idx = 0; var_idx < var_num; var_idx++) {
            diff_eq(mqs[eq_idx],
                    term_num,
                    var_num,
                    var_idx,
                    derivs + eq_idx * var_num * (var_num + 1) + var_idx * (var_num + 1));
        }
    }
}

void fast_exhaustive_search(uint8_t **mqs, uint64_t var_num, uint64_t eq_num, uint64_t xvar_num, uint8_t *solution) {
    uint64_t eq_idx, var_idx, i, term;
    bool derivs[eq_num][var_num][var_num + 1];
    find_partial_derivs(mqs, (bool *) derivs, eq_num, xvar_num, var_num);

    // the second order partial dervatives of the system
    // pdiff2[j][i]'s eq_idx-th bit, starting from LSB, holds the partial
    // dervatives of equation eq_idx against first x_i then x_j,
    // which is a constant.
    uint64_t pdiff2[var_num][var_num];
    memset(pdiff2, 0x0, sizeof(uint64_t) * var_num * var_num);
    for (var_idx = 0; var_idx < var_num; ++var_idx) {
        for (i = 0; i < var_num; ++i) {
            for (eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
                term = derivs[eq_idx][var_idx][i];
                assert(0x1UL == term || 0x0UL == term);
                pdiff2[i][var_idx] |= term << eq_idx;
            }
        }
    }

    uint64_t pdiff_eval[var_num];
    memset(pdiff_eval, 0x0, sizeof(uint64_t) * var_num);
    for (var_idx = 0; var_idx < var_num; ++var_idx) {
        for (eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
            if (0 == var_idx) {
                term = mqs[eq_idx][deg2midx1(var_num, 0)];
            } else {
                term = mqs[eq_idx][deg2midx1(var_num, var_idx)] ^ derivs[eq_idx][var_idx][var_idx - 1];
            }
            assert(0x1UL == term || 0x0UL == term);
            pdiff_eval[var_idx] |= term << eq_idx;
        }
    }

    printf("\t\tbrute forcing...\n");

    // func_eval's eq_idx-th bit (from LSB) holds the result of evaluating
    // equation eq_idx at the current solution
    uint64_t func_eval = 0x0UL;
    // evaluate all functions at initial solution, zero vector
    for (eq_idx = 0; eq_idx < eq_num; ++eq_idx) {
        term = mqs[eq_idx][xvar_num - 1];
        assert(0x1UL == term || 0x0UL == term);
        func_eval |= term << eq_idx;
    }

    uint64_t count = 0;
    const uint64_t bound = (0x1UL << var_num) - 1;
    while (func_eval && count < bound) {
        ++count;
        uint64_t fp_idx = __builtin_ctzll(count);

        if (count & (count - 1)) {
            uint64_t pre_fp_idx = __builtin_ctzll(count ^ (0x1UL << fp_idx));
            pdiff_eval[fp_idx] ^= pdiff2[fp_idx][pre_fp_idx];
        }

        func_eval ^= pdiff_eval[fp_idx];
    }

    // fill in the solution
    if (!func_eval) {
        printf("\t\tfound valid solution: %d\n", count);
        for (var_idx = 0; var_idx < var_num; ++var_idx) {
            solution[var_idx] = ((count ^ (count >> 1)) >> var_idx) & 0x1UL;
        }
    }
}

bool verify_sols(uint8_t **sys, uint8_t *solution, uint64_t var_num, uint64_t eq_num, uint64_t xvar_num) {
    for (uint64_t i = 0; i < eq_num; ++i) { // for each equation
        uint64_t res = 0;

        for (uint64_t mul_1 = 0; mul_1 < var_num; mul_1++) {
            for (uint64_t mul_2 = mul_1; mul_2 < var_num; mul_2++) {
                if (sys[i][deg2midx2(mul_1, mul_2)] == 1) {
                    res ^= (solution[mul_1] & solution[mul_2]);
                }
            }
        }

        for (uint64_t var_idx = 0; var_idx < var_num; var_idx++) {
            res ^= (solution[var_idx] & sys[i][deg2midx1(var_num, var_idx)]);
        }

        res ^= sys[i][xvar_num - 1];
        if (res == 1) { // the equation is evaluated to 1
            return false;
        }
    }

    return true;
}
