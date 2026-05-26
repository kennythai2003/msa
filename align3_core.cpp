#include <algorithm>
#include <string>
#include <cstring>
using std::max;
using std::string;

// the following idea is from claude for speedup
// SP[a][b][c] = s(a,b) + s(a,c) + s(b,c)
// we precompute for all 125 combinations so each cell does a lookup instead of recomputing
// A = 0 C = 1 G = 2 T = 3 gap = 4
static int SP[5][5][5];

static int s(int a, int b)
{
    if (a == 4 && b == 4)
    {
        return 0; // gap-gap
    }
    if (a == 4 || b == 4)
    {
        return -8; // indel
    }
    if (a == b)
    {
        return 5; // match
    }
    return -4; // mismatch
}

// called once at startup to fill SP
static void build_sp()
{
    for (int a = 0; a < 5; a++)
    {
        for (int b = 0; b < 5; b++)
        {
            for (int c = 0; c < 5; c++)
            {
                SP[a][b][c] = s(a, b) + s(a, c) + s(b, c);
            }
        }
    }
}

// used in the traceback to convert an integer code back to its nucleotide character
static const char DECODE[] = "ACGT-";

// encoding characters to integers so they index directly into SP
static inline int encode(char c)
{
    switch (c)
    {
    case 'A':
        return 0;
    case 'C':
        return 1;
    case 'G':
        return 2;
    case 'T':
        return 3;
    default:
        return 4;
    }
}

// computes the dp score plane after processing w[0..w_len-1]
// based on the observation in the textbook (p.232) that only the scores in the preceding column are needed to compute the next column, so earlier planes can be discarded

// we flatten the 2D plane into a 1D array for cache efficiency (idea from claude for speedup)
// cell (i,j) is at index i*stride + j
static int *forward_plane(const int *seqU, int m,
                          const int *seqV, int n,
                          const int *seqW, int w_len)
{
    int stride = n + 1;
    int plane_size = (m + 1) * (n + 1);

    int *prev = new int[plane_size];
    int *curr = new int[plane_size];

    // initialize the l=0 plane with w treated as all gaps
    prev[0] = 0;

    // left edge: u advances, v and w are gaps
    for (int i = 1; i <= m; i++)
    {
        prev[i * stride] = prev[(i - 1) * stride] + SP[seqU[i - 1]][4][4];
    }

    // top edge: v advances, u and w are gaps
    for (int j = 1; j <= n; j++)
    {
        prev[j] = prev[j - 1] + SP[4][seqV[j - 1]][4];
    }

    // interior: standard 2-sequence DP with w treated as all gaps (l=0 init plane)
    for (int i = 1; i <= m; i++)
    {
        int u_char = seqU[i - 1];
        for (int j = 1; j <= n; j++)
        {
            int v_char = seqV[j - 1];
            int cell = i * stride + j;
            int from_top = (i - 1) * stride + j;
            int from_left = i * stride + (j - 1);
            int from_diag = (i - 1) * stride + (j - 1);
            prev[cell] = max({prev[from_top] + SP[u_char][4][4], prev[from_left] + SP[4][v_char][4], prev[from_diag] + SP[u_char][v_char][4]});
        }
    }

    // advance one w character at a time, filling curr from prev
    for (int l = 0; l < w_len; l++)
    {
        int w_char = seqW[l];

        // corner: all three gapped except w
        curr[0] = prev[0] + SP[4][4][w_char];

        // left edge: u and w advance, v is gapped
        for (int i = 1; i <= m; i++)
        {
            int u_char = seqU[i - 1];
            int cell = i * stride;
            int from_above = (i - 1) * stride;
            curr[cell] = max({prev[cell] + SP[4][4][w_char], curr[from_above] + SP[u_char][4][4], prev[from_above] + SP[u_char][4][w_char]});
        }

        // top edge: v and w advance, u is gapped
        for (int j = 1; j <= n; j++)
        {
            int v_char = seqV[j - 1];
            int from_left = j - 1;
            curr[j] = max({prev[j] + SP[4][4][w_char],
                           curr[from_left] + SP[4][v_char][4],
                           prev[from_left] + SP[4][v_char][w_char]});
        }

        // interior: all 7 transitions from the 3-sequence recurrence (textbook p.179, Ch.6.10)
        for (int i = 1; i <= m; i++)
        {
            int u_char = seqU[i - 1];
            for (int j = 1; j <= n; j++)
            {
                int v_char = seqV[j - 1];
                int cell = i * stride + j;
                int above = (i - 1) * stride + j;
                int left = i * stride + (j - 1);
                int diag = (i - 1) * stride + (j - 1);
                // aligned vertically to be read easier
                curr[cell] = max({prev[cell] + SP[4][4][w_char],
                                  curr[above] + SP[u_char][4][4],
                                  curr[left] + SP[4][v_char][4],
                                  prev[above] + SP[u_char][4][w_char],
                                  prev[left] + SP[4][v_char][w_char],
                                  curr[diag] + SP[u_char][v_char][4],
                                  prev[diag] + SP[u_char][v_char][w_char]});
            }
        }

        // discard the older plane, curr becomes the new prev (textbook p.232, Fig 7.2)
        int *tmp = prev;
        prev = curr;
        curr = tmp;
    }

    delete[] curr;
    return prev;
}

// based on textbook p.233: suffix(i) is computed as the longest path in the reversed edit graph
// runs forward_plane on reversed sequences then flips the result so that
// flipped[i][j] = score of aligning u[i..m-1], v[j..n-1], w[w_mid..k-1]
// this gives the suffix scores needed for the middle vertex search
static int *backward_plane(const int *seqU, int m,
                           const int *seqV, int n,
                           const int *seqW, int w_len)
{
    int *revU = new int[m];
    int *revV = new int[n];
    int *revW = new int[w_len];

    for (int i = 0; i < m; i++)
    {
        revU[i] = seqU[m - 1 - i];
    }
    for (int j = 0; j < n; j++)
    {
        revV[j] = seqV[n - 1 - j];
    }
    for (int l = 0; l < w_len; l++)
    {
        revW[l] = seqW[w_len - 1 - l];
    }

    int *plane = forward_plane(revU, m, revV, n, revW, w_len);
    delete[] revU;
    delete[] revV;
    delete[] revW;

    // flip so flipped[i][j] = score of the suffix starting at (i,j)
    int stride = n + 1;
    int plane_size = (m + 1) * (n + 1);
    int *flipped = new int[plane_size];
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            flipped[i * stride + j] = plane[(m - i) * stride + (n - j)];
        }
    }

    delete[] plane;
    return flipped;
}

// textbook p.234: "if source and sink are in consecutive columns, output longest path from source to sink" & claude help
// full 3D DP with traceback.
// used as the base case in hirschberg() when any sequence has length 0 or 1.
// i had to use this because when i tried implementing, i would get stuck in an infinite loop or produce incorrect results (addressed again later in the code)
static void small_align(const int *seqU, int m,
                        const int *seqV, int n,
                        const int *seqW, int k,
                        string &aln_u, string &aln_v, string &aln_w)
{
    int N1 = n + 1, K1 = k + 1;
    long long sz = (long long)(m + 1) * (n + 1) * (k + 1);
    int *dp = new int[sz];
    char *tb = new char[sz]; // traceback direction bits: di|dj|dl

// did this so we don't have to write the same inline expression a bunch of times
#define IDX(i, j, l) ((long long)(i) * N1 * K1 + (j) * K1 + (l))

    dp[IDX(0, 0, 0)] = 0;
    tb[IDX(0, 0, 0)] = 0;
    // u advances, v and w are gaps
    for (int i = 1; i <= m; i++)
    {
        dp[IDX(i, 0, 0)] = dp[IDX(i - 1, 0, 0)] + SP[seqU[i - 1]][4][4];
        tb[IDX(i, 0, 0)] = 4;
    }
    // v advances, u and w are gaps
    for (int j = 1; j <= n; j++)
    {
        dp[IDX(0, j, 0)] = dp[IDX(0, j - 1, 0)] + SP[4][seqV[j - 1]][4];
        tb[IDX(0, j, 0)] = 2;
    }
    // w advances, u and v are gaps
    for (int l = 1; l <= k; l++)
    {
        dp[IDX(0, 0, l)] = dp[IDX(0, 0, l - 1)] + SP[4][4][seqW[l - 1]];
        tb[IDX(0, 0, l)] = 1;
    }
    // i-j face: 2-sequence DP between u and v, w is all gaps
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            int best = dp[IDX(i - 1, j, 0)] + SP[seqU[i - 1]][4][4];
            char dir = 4;
            int v;
            v = dp[IDX(i, j - 1, 0)] + SP[4][seqV[j - 1]][4];
            if (v > best)
            {
                best = v;
                dir = 2;
            }
            v = dp[IDX(i - 1, j - 1, 0)] + SP[seqU[i - 1]][seqV[j - 1]][4];
            if (v > best)
            {
                best = v;
                dir = 6;
            }
            dp[IDX(i, j, 0)] = best;
            tb[IDX(i, j, 0)] = dir;
        }
    }
    // i-l face: 2-sequence DP between u and w, v is all gaps
    for (int i = 1; i <= m; i++)
    {
        for (int l = 1; l <= k; l++)
        {
            int best = dp[IDX(i - 1, 0, l - 1)] + SP[seqU[i - 1]][4][seqW[l - 1]];
            char dir = 5;
            int v;
            v = dp[IDX(i - 1, 0, l)] + SP[seqU[i - 1]][4][4];
            if (v > best)
            {
                best = v;
                dir = 4;
            }
            v = dp[IDX(i, 0, l - 1)] + SP[4][4][seqW[l - 1]];
            if (v > best)
            {
                best = v;
                dir = 1;
            }
            dp[IDX(i, 0, l)] = best;
            tb[IDX(i, 0, l)] = dir;
        }
    }
    // j-l face: 2-sequence DP between v and w, u is all gaps
    for (int j = 1; j <= n; j++)
    {
        for (int l = 1; l <= k; l++)
        {
            int best = dp[IDX(0, j - 1, l - 1)] + SP[4][seqV[j - 1]][seqW[l - 1]];
            char dir = 3;
            int v;
            v = dp[IDX(0, j - 1, l)] + SP[4][seqV[j - 1]][4];
            if (v > best)
            {
                best = v;
                dir = 2;
            }
            v = dp[IDX(0, j, l - 1)] + SP[4][4][seqW[l - 1]];
            if (v > best)
            {
                best = v;
                dir = 1;
            }
            dp[IDX(0, j, l)] = best;
            tb[IDX(0, j, l)] = dir;
        }
    }
    // 3D dp all three sequences advance
    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            for (int l = 1; l <= k; l++)
            {
                int best = dp[IDX(i - 1, j - 1, l - 1)] + SP[seqU[i - 1]][seqV[j - 1]][seqW[l - 1]];
                char dir = 7;
                int v;
                v = dp[IDX(i - 1, j - 1, l)] + SP[seqU[i - 1]][seqV[j - 1]][4];
                if (v > best)
                {
                    best = v;
                    dir = 6;
                }
                v = dp[IDX(i - 1, j, l - 1)] + SP[seqU[i - 1]][4][seqW[l - 1]];
                if (v > best)
                {
                    best = v;
                    dir = 5;
                }
                v = dp[IDX(i, j - 1, l - 1)] + SP[4][seqV[j - 1]][seqW[l - 1]];
                if (v > best)
                {
                    best = v;
                    dir = 3;
                }
                v = dp[IDX(i - 1, j, l)] + SP[seqU[i - 1]][4][4];
                if (v > best)
                {
                    best = v;
                    dir = 4;
                }
                v = dp[IDX(i, j - 1, l)] + SP[4][seqV[j - 1]][4];
                if (v > best)
                {
                    best = v;
                    dir = 2;
                }
                v = dp[IDX(i, j, l - 1)] + SP[4][4][seqW[l - 1]];
                if (v > best)
                {
                    best = v;
                    dir = 1;
                }
                dp[IDX(i, j, l)] = best;
                tb[IDX(i, j, l)] = dir;
            }
        }
    }

    // traceback builds the alignment in reverse order then reverses at the end
    string ru, rv, rw;
    int i = m, j = n, l = k;
    while (i > 0 || j > 0 || l > 0)
    {
        char d = tb[IDX(i, j, l)];
        int di = (d >> 2) & 1;
        int dj = (d >> 1) & 1;
        int dl = d & 1;
        ru += di ? DECODE[seqU[i - 1]] : '-';
        rv += dj ? DECODE[seqV[j - 1]] : '-';
        rw += dl ? DECODE[seqW[l - 1]] : '-';
        i -= di;
        j -= dj;
        l -= dl;
    }
    for (int x = ru.size() - 1; x >= 0; x--)
    {
        aln_u += ru[x];
        aln_v += rv[x];
        aln_w += rw[x];
    }

    delete[] dp;
    delete[] tb;
#undef IDX
}

// hirschberg divide and conquer follows PATH(source, sink) from the textbook (p.234)
// extended to three sequences by dividing along w instead of a single column.

// textbook p.233: length(i) = prefix(i) + suffix(i), find mid = argmax length(i)
// textbook p.233: length(i,j) = fwd[i][j] + bwd[i][j], find (i*,j*) = argmax length(i,j)
// textbook p.234: PATH(source, mid) then PATH(mid, sink)
// hirschberg on left half then hirschberg on right half
static void hirschberg(const int *seqU, int m,
                       const int *seqV, int n,
                       const int *seqW, int k,
                       string &aln_u, string &aln_v, string &aln_w)
{
    // base case: empty sequence, align remaining two sequences against each other with gaps for the empty one
    // textbook p.234: source and sink in consecutive columns
    if (m == 0)
    {
        // u is empty: optimally align v and w against each other, padding u with gaps
        small_align(seqU, 0, seqV, n, seqW, k, aln_u, aln_v, aln_w);
        return;
    }
    if (n == 0)
    {
        // v is empty: optimally align u and w against each other, padding v with gaps
        small_align(seqU, m, seqV, 0, seqW, k, aln_u, aln_v, aln_w);
        return;
    }
    if (k == 0)
    {
        // w is empty: optimally align u and v against each other, padding w with gaps
        small_align(seqU, m, seqV, n, seqW, 0, aln_u, aln_v, aln_w);
        return;
    }

    // divide w at its midpoint (textbook p.233: middle column m/2)
    // when k=1, w_mid=0 so the first recursive call would receive the same k=1 again,
    // causing infinite recursion. same applies when m=1 or n=1.
    // small_align solves these tiny subproblems with a full 3D DP and traceback.
    // textbook p.234: "if source and sink are in consecutive columns, output longest path"
    if (k == 1 || m == 1 || n == 1)
    {
        small_align(seqU, m, seqV, n, seqW, k, aln_u, aln_v, aln_w);
        return;
    }
    int w_mid = k / 2;

    // prefix scores: fwd[i][j] = best score aligning u[0..i-1], v[0..j-1], w[0..w_mid-1]
    // textbook p.233: prefix(i) = s_{i, m/2}
    int *fwd = forward_plane(seqU, m, seqV, n, seqW, w_mid);

    // suffix scores: bwd[i][j] = best score aligning u[i..m-1], v[j..n-1], w[w_mid..k-1]
    // textbook p.233: suffix(i) computed via reversed edit graph
    int *bwd = backward_plane(seqU, m, seqV, n, seqW + w_mid, k - w_mid);

    // find the middle vertex (i*, j*) that maximises fwd[i][j] + bwd[i][j]
    // textbook p.233: mid = argmax_{0<=i<=n} length(i) = prefix(i) + suffix(i)
    int stride = n + 1;
    int best_score = fwd[0] + bwd[0];
    int i_star = 0;
    int j_star = 0;
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            int total = fwd[i * stride + j] + bwd[i * stride + j];
            if (total > best_score)
            {
                best_score = total;
                i_star = i;
                j_star = j;
            }
        }
    }

    delete[] fwd;
    delete[] bwd;

    // recurse on both halves (textbook p.234: PATH(source, mid) then PATH(mid, sink))
    hirschberg(seqU, i_star, seqV, j_star, seqW, w_mid, aln_u, aln_v, aln_w);
    hirschberg(seqU + i_star, m - i_star, seqV + j_star, n - j_star, seqW + w_mid, k - w_mid, aln_u, aln_v, aln_w);
}

extern "C"
{
    int msa_align3(const char *u, int m,
                   const char *v, int n,
                   const char *w, int k,
                   int *aln_len, int *perfect,
                   char *out_u, char *out_v, char *out_w)
    {
        build_sp();

        int seq_len_u = (m > 0) ? m : 1;
        int seq_len_v = (n > 0) ? n : 1;
        int seq_len_w = (k > 0) ? k : 1;

        int *seqU = new int[seq_len_u];
        int *seqV = new int[seq_len_v];
        int *seqW = new int[seq_len_w];

        for (int i = 0; i < m; i++)
        {
            seqU[i] = encode(u[i]);
        }
        for (int j = 0; j < n; j++)
        {
            seqV[j] = encode(v[j]);
        }
        for (int l = 0; l < k; l++)
        {
            seqW[l] = encode(w[l]);
        }

        string aln_u, aln_v, aln_w;
        aln_u.reserve(m + n + k);
        aln_v.reserve(m + n + k);
        aln_w.reserve(m + n + k);
        hirschberg(seqU, m, seqV, n, seqW, k, aln_u, aln_v, aln_w);

        int score = 0;
        int length = 0;
        int perf = 0;
        for (int col = 0; col < (int)aln_u.size(); col++)
        {
            int a = encode(aln_u[col]);
            int b = encode(aln_v[col]);
            int c = encode(aln_w[col]);
            score += SP[a][b][c];
            length++;
            if (a == b && b == c && a != 4)
            {
                perf++;
            }
        }

        *aln_len = length;
        *perfect = perf;
        memcpy(out_u, aln_u.c_str(), aln_u.size() + 1);
        memcpy(out_v, aln_v.c_str(), aln_v.size() + 1);
        memcpy(out_w, aln_w.c_str(), aln_w.size() + 1);

        delete[] seqU;
        delete[] seqV;
        delete[] seqW;
        return score;
    }
}