#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Streaming summary of the supplied CIZ daily CSV. No raw rows are copied. */
#define MAX_DATE_KEY 300000
#define BUF_SIZE (8 * 1024 * 1024)

typedef struct { double n, sx, sy, sxx, sxy; } Cov;
static double market[MAX_DATE_KEY], lag_market[MAX_DATE_KEY], lead_market[MAX_DATE_KEY];
static Cov covs[3];
static double last5[5];
static int ring_count, ring_pos, current_permno = -1, current_month = -1, current_year = -1;
static FILE *tail_out, *beta_out;
static unsigned long long rows, market_mismatch, duplicate_stock_days;
static int previous_permno = -1, previous_key = -1;

static void capture(char *dst, size_t cap, size_t *used, char c) {
    if (*used + 1 < cap) { dst[(*used)++] = c; dst[*used] = 0; }
}
static void reset(char *a, char *b, char *c, char *d,
                  size_t *na, size_t *nb, size_t *nc, size_t *nd) {
    a[0] = b[0] = c[0] = d[0] = 0;
    *na = *nb = *nc = *nd = 0;
}
static void finish_month(void) {
    if (current_permno < 0 || current_month < 0) return;
    if (ring_count >= 5) {
        double gross = 1.0;
        for (int i = 0; i < 5; i++) gross *= last5[i];
        if (isfinite(gross)) fprintf(tail_out, "%d,%d,%.12g\n", current_permno, current_month, gross);
    }
    ring_count = ring_pos = 0;
}
static void finish_year(void) {
    if (current_permno < 0 || current_year < 0) return;
    fprintf(beta_out, "%d,%d", current_permno, current_year);
    for (int j = 0; j < 3; j++)
        fprintf(beta_out, ",%.0f,%.12g,%.12g,%.12g,%.12g",
                covs[j].n, covs[j].sx, covs[j].sy, covs[j].sxx, covs[j].sxy);
    fputc('\n', beta_out);
    memset(covs, 0, sizeof(covs));
}
static int date_key(const char *s) {
    if (strlen(s) < 10) return -1;
    int y = atoi(s), m = atoi(s + 5), d = atoi(s + 8);
    int k = (y - 2000) * 10000 + m * 100 + d;
    return k > 0 && k < MAX_DATE_KEY ? k : -1;
}
static void record(int pass, const char *pno, const char *date, const char *ret, const char *mkt) {
    int k = date_key(date);
    if (k < 0 || !*pno) return;
    rows++;
    if (pass == 1) {
        if (*mkt) {
            double x = strtod(mkt, NULL);
            if (isfinite(x)) {
                if (isfinite(market[k]) && fabs(market[k] - x) > 0.000002) market_mismatch++;
                else market[k] = x;
            }
        }
        return;
    }
    int permno = atoi(pno), year = atoi(date), month = year * 100 + atoi(date + 5);
    if (permno == previous_permno && k == previous_key) {
        duplicate_stock_days++;
        return;
    }
    previous_permno = permno; previous_key = k;
    if (permno != current_permno || month != current_month) finish_month();
    if (permno != current_permno || year != current_year) finish_year();
    current_permno = permno; current_month = month; current_year = year;
    double y = *ret ? strtod(ret, NULL) : NAN;
    double gross = isfinite(y) && y >= -1.0 ? 1.0 + y : NAN;
    last5[ring_pos] = gross; ring_pos = (ring_pos + 1) % 5; ring_count++;
    double xs[3] = {lag_market[k], market[k], lead_market[k]};
    if (isfinite(y)) for (int j = 0; j < 3; j++) if (isfinite(xs[j])) {
        Cov *c = &covs[j]; double x = xs[j];
        c->n++; c->sx += x; c->sy += y; c->sxx += x*x; c->sxy += x*y;
    }
}
static int scan(const char *path, int pass) {
    FILE *f = fopen(path, "rb"); if (!f) { perror(path); return 1; }
    char *buf = malloc(BUF_SIZE); if (!buf) return 1;
    char pno[24], date[24], ret[40], mkt[40];
    size_t na, nb, nc, nd, nr; int field = 1, quoted = 0, header = 1;
    reset(pno, date, ret, mkt, &na, &nb, &nc, &nd);
    while ((nr = fread(buf, 1, BUF_SIZE, f)) > 0) for (size_t i = 0; i < nr; i++) {
        char c = buf[i];
        if (c == '"') quoted = !quoted;
        else if (c == ',' && !quoted) field++;
        else if (c == '\n' && !quoted) {
            if (header) header = 0; else record(pass, pno, date, ret, mkt);
            field = 1; reset(pno, date, ret, mkt, &na, &nb, &nc, &nd);
        } else if (c != '\r') {
            if (field == 1) capture(pno, sizeof(pno), &na, c);
            else if (field == 37) capture(date, sizeof(date), &nb, c);
            else if (field == 48) capture(ret, sizeof(ret), &nc, c);
            else if (field == 90 && pass == 1) capture(mkt, sizeof(mkt), &nd, c);
        }
    }
    if (na && !header) record(pass, pno, date, ret, mkt);
    if (pass == 2) { finish_month(); finish_year(); }
    free(buf); fclose(f); return 0;
}
int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: daily_summary daily.csv tail.csv beta_sums.csv market.csv\n"); return 2; }
    for (int k = 0; k < MAX_DATE_KEY; k++) market[k] = lag_market[k] = lead_market[k] = NAN;
    if (scan(argv[1], 1)) return 1;
    fprintf(stderr, "market pass: %llu rows, %llu inconsistent market values\n", rows, market_mismatch);
    FILE *market_out = fopen(argv[4], "w"); if (!market_out) return 1;
    fprintf(market_out, "date_key,vwretd\n");
    double prev = NAN; int prev_year = -1;
    for (int k = 0; k < MAX_DATE_KEY; k++) if (isfinite(market[k])) {
        int year = 2000 + k / 10000;
        if (year != prev_year) prev = NAN;
        lag_market[k] = prev; prev = market[k]; prev_year = year;
        fprintf(market_out, "%d,%.12g\n", 20000000 + k, market[k]);
    }
    fclose(market_out);
    double next = NAN; int next_year = -1;
    for (int k = MAX_DATE_KEY - 1; k >= 0; k--) if (isfinite(market[k])) {
        int year = 2000 + k / 10000;
        if (year != next_year) next = NAN;
        lead_market[k] = next; next = market[k]; next_year = year;
    }
    tail_out = fopen(argv[2], "w"); beta_out = fopen(argv[3], "w");
    if (!tail_out || !beta_out) return 1;
    fprintf(tail_out, "PERMNO,YYYYMM,last5_gross\n");
    fprintf(beta_out, "PERMNO,year");
    for (int j = 0; j < 3; j++) fprintf(beta_out, ",n%d,sx%d,sy%d,sxx%d,sxy%d", j,j,j,j,j);
    fputc('\n', beta_out);
    rows = 0;
    if (scan(argv[1], 2)) return 1;
    fprintf(stderr, "stock pass: %llu rows, %llu duplicate security-days skipped\n", rows, duplicate_stock_days);
    fclose(tail_out); fclose(beta_out);
    return 0;
}
