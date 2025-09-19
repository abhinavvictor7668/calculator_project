#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define f(x) (cos(x) - (x * exp(x)))
#define df(x) (-sin(x) - (x) * exp(x) - exp(x))
#define ddf(x) (-cos(x) - exp(x) - (x) * exp(x))
#define g(x) 1.0 / (1 + x)
#define MAX 500
#define N 3
#define EPSILON 0.0000001

int compare(const void *a, const void *b)
{
    long double diff = *(long double *)a - *(long double *)b;
    return (diff > 0) - (diff < 0);
}

void AdditionSubtraction1(void)
{
    int i = 0, n;
    long double arr[MAX];
    long double sum = 0;

    printf("Enter real numbers one by one\n");
    printf("To stop entering values, enter any an alphabet (like 'q')\n");

    for (i = 0; i < MAX; i++)
    {
        if (scanf("%Lf", &arr[i]) != 1)
        {
            while (getchar() != '\n')
                ;
            break;
        }
        sum = sum + arr[i];
    }
    n = i;
    printf("\nAnswer = %.2Lf\n", sum);
}

void Multiplication2(void)
{
    int i = 0, n;
    long double arr[MAX];
    long double prod = 1;

    printf("Enter real numbers one by one\n");
    printf("To stop entering values, enter any an alphabet (like 'q')\n");

    for (i = 0; i < MAX; i++)
    {
        if (scanf("%Lf", &arr[i]) != 1)
        {
            while (getchar() != '\n')
                ;
            break;
        }
        prod = prod * arr[i];
    }
    printf("\nAnswer = %.2Lf\n", prod);
}

void Divide2Nos3(void)
{
    long double a, b, c;

    printf("\nEnter real numbers\n");
    scanf("%Lf %Lf", &a, &b);

    if (b == 0)
    {
        printf("\nDivision by 0 not defined\n");
    }

    else
    {
        c = a / b;
        printf("\nAnswer = %.2Lf\n", c);
    }
}

void Exponent4(void)
{
    long double num, multi = 1;
    int pow;

    printf("\nEnter any real base and an integer exponent\n");
    scanf("%Lf %d", &num, &pow);

    if (num == 0)
    {
        if (pow > 0)
        {
            printf("\n0\n");
        }

        else
        {
            printf("\nNot defined\n");
        }
    }

    else
    {
        if (pow > 0)
        {
            while (pow > 0)
            {
                multi = (multi * num);
                pow = pow - 1;
            }
            printf("\n%.2Lf\n", multi);
        }

        else if (pow < 0)
        {
            while (pow < 0)
            {
                multi = (multi * num);
                pow = pow + 1;
            }
            printf("\n%.2Lf\n", (1 / multi));
        }

        else
        {
            printf("\n1\n");
        }
    }
}

void Factorial5(void)
{
    long double count, fact, multi = 1;

    printf("\nEnter the natural number for which factorial is to be computed\n");
    scanf("%Lf", &fact);

    if (fact < 0)
    {
        printf("\nFactorial of negative numbers does not exist\n");
    }

    else
    {
        for (count = 1; count <= fact; count++)
        {
            multi = multi * count;
        }
        printf("\nAnswer = %.2Lf\n", multi);
    }
}

void Permutation6(void)
{
    long double count, n, r, n_r, multi1 = 1,
                                  multi2 = 1, multi3 = 1, per, com;

    printf("\nEnter natural N and R for computation of nPr\n");
    scanf("%Lf %Lf", &n, &r);

    n_r = (n - r);

    if (n < r)
    {
        printf("\nCalculation not possible\n");
    }

    else
    {
        for (count = 1; count <= n; count++)
        {
            multi1 = multi1 * count;
        }

        for (count = 1; count <= r; count++)
        {
            multi2 = multi2 * count;
        }

        for (count = 1; count <= n_r; count++)
        {
            multi3 = multi3 * count;
        }
    }

    per = (multi1 / multi3);

    printf("\nAnswer = %.2Lf\n", per);
}

void Combination7(void)
{
    long double count, n, r, n_r, multi1 = 1,
                                  multi2 = 1, multi3 = 1, per, com;

    printf("\nEnter natural N and R for computation of nCr\n");
    scanf("%Lf %Lf", &n, &r);

    n_r = (n - r);

    if (n < r)
    {
        printf("\nCalculation not possible\n");
    }

    else
    {
        for (count = 1; count <= n; count++)
        {
            multi1 = multi1 * count;
        }

        for (count = 1; count <= r; count++)
        {
            multi2 = multi2 * count;
        }

        for (count = 1; count <= n_r; count++)
        {
            multi3 = multi3 * count;
        }
    }

    com = (multi1 / (multi2 * multi3));

    printf("\nAnswer = %.2Lf\n", com);
}

void AP8(void)
{
    long double a, d, term, sum = 0;
    int i, n;

    printf("\nEnter first term, common difference, number of terms to be generated\n");
    scanf("%Lf %Lf %d", &a, &d, &n);

    for (i = 1; i <= n; i++)
    {
        term = (a + ((i - 1) * d));
        printf("%Lf ", term);

        sum = sum + term;
    }
    printf("\nAnswer(sum) = %.3Lf\n", sum);
}

void GP9(void)
{
    long double a, r, term, sum = 0;
    int i, n;

    printf("\nEnter first term, common ratio, number of terms to be generated\n");
    scanf("%Lf %Lf %d", &a, &r, &n);

    term = a;
    for (i = 1; i <= n; i++)
    {
        printf("%Lf ", term);
        sum = sum + term;
        term = term * r;
    }
    printf("\nAnswer(sum) = %.2Lf\n", sum);
}

void HP10(void)
{
    long double a, d, term, sum = 0;
    long long int i, n;

    printf("\nEnter first term, common difference, number of terms to be generated\n");
    scanf("%Lf %Lf %lld", &a, &d, &n);

    for (i = 1; i <= n; i++)
    {
        term = 1.0 / (a + ((i - 1) * d));
        printf("%Lf ", term);

        sum = sum + term;
    }
    printf("\nAnswer(sum) = %.2Lf\n", sum);
}

void SimpleInterest11(void)
{
    long double Rate, SI, Principal, Total;
    int Years;

    printf("\nEnter Principal amount, Rate and Years\n");
    scanf("%Lf %Lf %d", &Principal, &Rate, &Years);

    SI = (Principal * Rate * Years) / 100;
    Total = Principal + SI;

    printf("Principal depostited = Rs.%.2Lf\n", Principal);
    printf("Rate of interest = %.2LF%%\n", Rate);
    printf("Deposite period is %d years\n", Years);
    printf("Simple Interest earned = Rs.%.2Lf\n", SI);
    printf("Total amount = Rs.%.2Lf\n", Total);
}

void CompoundInterest12(void)
{
    long double base, p, r, multi = 1;
    int pow, n, t;

    // P * (1 + r/n)^nt

    printf("\nEnter Principal, Frequency, Rate and Time\n");
    scanf("%Lf %d %Lf %d", &p, &n, &r, &t);

    r = r * 0.01;
    base = (1 + (r / n));
    pow = n * t;

    while (pow > 0)
    {
        multi = (multi * base);
        pow = pow - 1;
    }
    printf("Principal depostited = Rs.%.2Lf\n", p);
    printf("Rate of interest = %.2Lf%%\n", (r * 100));
    printf("Compouning occurs %d time(s)\n", n);
    printf("Deposite period is %d years\n", t);
    printf("Compound Interest earned = Rs.%.2Lf\n", ((p * multi) - p));
    printf("Total amount = Rs.%.2Lf\n", (p * multi));
}

void QuadraticRoots13(void)
{
    long double a, b, c, dis, root1, root2, rp, ip;

    printf("\nEnter the quadratic equation coefficients\n");
    scanf("%Lf %Lf %Lf", &a, &b, &c);

    dis = (b * b) - (4 * a * c);

    if (dis >= 0)
    {
        root1 = ((-1 * b) + sqrt(dis)) / (2 * a);
        root2 = ((-1 * b) - sqrt(dis)) / (2 * a);

        printf("\nFirst real root = %.3Lf\n", root1);
        printf("\nSecond real root = %.3Lf\n", root2);
    }

    else
    {
        rp = ((-1 * b) / (2 * a));
        ip = ((sqrt(-1 * dis)) / (2 * a));

        printf("\nFirst imaginary root = %.3Lf + %.3Lfi\n", rp, ip);
        printf("\nSecond imaginary root = %.3Lf - %.3Lfi\n", rp, ip);
    }
}

void HCF14(void)
{
    unsigned long long int i, num1, num2, hcf;

    printf("\nEnter 2 natural numbers\n");
    scanf("%llu %llu", &num1, &num2);

    if (num1 < num2)
    {
        for (i = num1; i >= 1; i--)
        {
            if ((num1 % i == 0) && (num2 % i == 0))
            {
                hcf = i;
                break;
            }
        }
        printf("\nAnswer = %llu\n", hcf);
    }

    if (num2 < num1)
    {
        for (i = num2; i >= 1; i--)
        {
            if ((num1 % i == 0) && (num2 % i == 0))
            {
                hcf = i;
                break;
            }
        }
        printf("\nAnswer = %llu\n", hcf);
    }

    if (num1 == num2)
    {
        printf("\nAnswer = %llu\n", num1);
    }
}

void LCM15(void)
{
    unsigned long long int i, num1, num2, lcm;

    printf("\nEnter 2 natural numbers\n");
    scanf("%llu %llu", &num1, &num2);

    if (num1 < num2)
    {
        for (i = num2; i >= num2; i++)
        {
            if ((i % num1 == 0) && (i % num2 == 0))
            {
                lcm = i;
                break;
            }
        }
        printf("\nAnswer = %llu\n", lcm);
    }

    if (num2 < num1)
    {
        for (i = num1; i >= num1; i++)
        {
            if ((i % num1 == 0) && (i % num2 == 0))
            {
                lcm = i;
                break;
            }
        }
        printf("\nAnswer = %llu\n", lcm);
    }

    if (num1 == num2)
    {
        printf("\nAnswer = %llu\n", num1);
    }
}

void MatrixAdd16(void)
{
    long double a[MAX][MAX], b[MAX][MAX], c[MAX][MAX];
    int nr = 0, nc = 0;
    int row = 0, col = 0;

    printf("\nEnter number of rows and columns\n");
    scanf("%d %d", &nr, &nc);

    if (nr > MAX || nc > MAX)
    {
        printf("\nThe row or column exceeds the dimension limits (%dx%d)\n", MAX, MAX);
    }

    printf("\nEnter the elements of first matrix\n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            scanf("%Lf", &a[row][col]);
        }
    }

    printf("\nEnter elements of second matrix:\n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            scanf("%Lf", &b[row][col]);
        }
    }

    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            c[row][col] = a[row][col] + b[row][col];
        }
    }

    printf("\nAnswer = \n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            printf("  %.2Lf  ", c[row][col]);
        }
        printf("\n");
    }
}

void MatrixSubtract17(void)
{
    long double a[MAX][MAX], b[MAX][MAX], c[MAX][MAX];
    int nr = 0, nc = 0;
    int row = 0, col = 0;

    printf("\nEnter number of rows and columns\n");
    scanf("%d %d", &nr, &nc);

    if (nr > MAX || nc > MAX)
    {
        printf("\nThe row or column exceeds the dimension limits (%dx%d)\n", MAX, MAX);
    }

    printf("\nEnter the elements of first matrix\n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            scanf("%Lf", &a[row][col]);
        }
    }

    printf("\nEnter elements of second matrix:\n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            scanf("%Lf", &b[row][col]);
        }
    }

    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            c[row][col] = a[row][col] - b[row][col];
        }
    }

    printf("\nAnswer = \n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            printf("  %.2Lf  ", c[row][col]);
        }
        printf("\n");
    }
}

void MatrixTranspose18(void)
{
    long double a[MAX][MAX], t[MAX][MAX];
    int nr, nc;
    int row, col;

    printf("Enter number of rows and columns: ");
    scanf("%d %d", &nr, &nc);

    if (nr > MAX || nc > MAX)
    {
        printf("Matrix size exceeds allowed dimensions (%dx%d).\n", MAX, MAX);
    }

    printf("\nEnter elements of the matrix:\n");
    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            scanf("%Lf", &a[row][col]);
        }
    }

    for (row = 0; row < nr; row++)
    {
        for (col = 0; col < nc; col++)
        {
            t[col][row] = a[row][col];
        }
    }

    printf("\nAnswer =\n");
    for (row = 0; row < nc; row++)
    {
        for (col = 0; col < nr; col++)
        {
            printf("  %.2Lf  ", t[row][col]);
        }
        printf("\n");
    }
}

void MatrixMulti19(void)
{
    long double a[MAX][MAX], b[MAX][MAX], c[MAX][MAX];
    int m, n, p;
    int row, col, k, sum;

    printf("\nEnter number of rows and columns of first matrix:\n");
    scanf("%d %d", &m, &n);

    printf("\nEnter number of columns of second matrix (rows is %d):\n", n);
    scanf("%d", &p);

    if (m > MAX || n > MAX || p > MAX)
    {
        printf("\nMatrix size exceeds allowed dimensions (%dx%d)\n", MAX, MAX);
    }

    printf("\nEnter elements of first matrix:\n");
    for (row = 0; row < m; row++)
    {
        for (col = 0; col < n; col++)
        {
            scanf("%Lf", &a[row][col]);
        }
    }

    printf("\nEnter elements of second matrix:\n");
    for (row = 0; row < n; row++)
    {
        for (col = 0; col < p; col++)
        {
            scanf("%Lf", &b[row][col]);
        }
    }

    for (row = 0; row < m; row++)
    {
        for (col = 0; col < p; col++)
        {
            sum = 0;
            for (k = 0; k < n; k++)
            {
                sum += a[row][k] * b[k][col];
            }
            c[row][col] = sum;
        }
    }

    printf("\nAnswer = \n");
    for (row = 0; row < m; row++)
    {
        for (col = 0; col < p; col++)
        {
            printf("  %.2Lf  ", c[row][col]);
        }
        printf("\n");
    }
}
////////////////////////////////////////////////////

void Raw(void)
{
    int i = 0, n;
    long double arr[MAX];
    long double modes[MAX];
    long double sum = 0, am, gm, prod = 1.0, hm, rsum = 0, median, q1, q3,
                range, iqr, meandevmean, meandevmedian, var, sd, cv,
                mu3, mu4, beta1, beta2, gamma1, gamma2;

    printf("Enter real numbers one by one\n");
    printf("To stop entering values, enter any an alphabet (like 'q')\n");

    // Data input
    for (i = 0; i < MAX; i++)
    {
        if (scanf("%Lf", &arr[i]) != 1)
        {
            while (getchar() != '\n')
                ;
            break;
        }
        sum = sum + arr[i];
    }
    n = i;
    if (n == 0)
    {
        printf("NO DATA ENTERED");
        return;
    }

    qsort(arr, n, sizeof(long double), compare);

    // AM
    am = sum / n;
    //

    // GM
    for (i = 0; i < n; i++)
    {
        prod = prod * arr[i];

        if (arr[i] > 0)
        {
            gm = powl(prod, 1.0 / n);
        }
        else
        {
            gm = 0;
        }
    }
    //

    // HM
    for (i = 0; i < n; i++)
    {
        if (arr[i] == 0)
        {
            hm = 0;
        }
        else
        {
            rsum = rsum + (1.0 / arr[i]);
            hm = n / rsum;
        }
    }
    //

    // Median
    if (n % 2 == 0)
    {
        median = (arr[n / 2] + arr[n / 2 - 1]) / 2;
    }
    else
    {
        median = arr[n / 2];
    }
    //

    // Mode
    int max_count = 1, count = 1, mode = 0;
    modes[0] = arr[0];

    for (i = 1; i < n; i++)
    {
        if (arr[i] == arr[i - 1])
        {
            count++;
        }
        else
        {
            if (count > max_count)
            {
                max_count = count;
            }
            count = 1;
        }
    }
    if (count > max_count)
    {
        max_count = count;
    }

    count = 1;
    for (i = 1; i < n; i++)
    {
        if (arr[i] == arr[i - 1])
        {
            count++;
        }
        else
        {
            if (count == max_count)
            {
                modes[mode++] = arr[i - 1];
            }
            count = 1;
        }
    }

    if (count == max_count)
    {
        modes[mode++] = arr[n - 1];
    }
    //

    long double pos1 = (n + 1) / 4.0, pos2 = 3 * (n + 1) / 4.0;

    // Q1
    int index1 = (int)pos1;
    long double frac1 = pos1 - index1;

    if (index1 <= 0)
    {
        q1 = arr[0];
    }
    else if (index1 >= n)
    {
        q1 = arr[n - 1];
    }
    else
    {
        q1 = arr[index1 - 1] + frac1 * (arr[index1] - arr[index1 - 1]);
    }
    //

    // Q3
    int index2 = (int)pos2;
    long double frac2 = pos2 - index2;

    if (index2 <= 0)
    {
        q3 = arr[0];
    }
    else if (index2 >= n)
    {
        q3 = arr[n - 1];
    }
    else
    {
        q3 = arr[index2 - 1] + frac2 * (arr[index2] - arr[index2 - 1]);
    }
    //

    // Range
    range = arr[n - 1] - arr[0];
    //

    // IQR
    iqr = q3 - q1;
    //

    // MeanDev mean
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + fabsl(arr[i] - am);
    }
    meandevmean = sum / n;
    //

    // MeanDev median
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + fabsl(arr[i] - median);
    }
    meandevmedian = sum / n;
    //

    // Variance
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + ((arr[i] - am) * (arr[i] - am));
    }
    var = sum / n;
    //

    // Stand. Dev
    sd = sqrtl(var);
    //

    // CV
    cv = sd / am;
    //

    // Skewness and Kurtosis from moments
    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + ((arr[i] - am) * (arr[i] - am) * (arr[i] - am));
    }
    mu3 = sum / n;

    sum = 0;
    for (i = 0; i < n; i++)
    {
        sum = sum + ((arr[i] - am) * (arr[i] - am) * (arr[i] - am) * (arr[i] - am));
    }
    mu4 = sum / n;

    if (mu3 > 0)
    {
        beta1 = ((mu3 * mu3) / (var * var * var));
        gamma1 = sqrt(beta1);
    }

    if (mu3 == 0)
    {
        beta1 = 0;
        gamma1 = 0;
    }

    if (mu3 < 0)
    {
        beta1 = ((mu3 * mu3) / (var * var * var));
        gamma1 = -sqrt(beta1);
        beta1 = -beta1;
    }

    beta2 = (mu4 / (var * var));
    gamma2 = beta2 - 3;
    //

    // RESULT
    printf("----CENTRAL TENDANCIES----\n");
    printf("AM\t\t\t%.3Lf\n", am);
    printf("GM\t\t\t%.3Lf\n", gm);
    printf("HM\t\t\t%.3Lf\n", hm);
    printf("Median\t\t\t%.3Lf\n", median);
    printf("Mode(s)\t\t\t");
    for (i = 0; i < mode; i++)
    {
        printf("%.3Lf   ", modes[i]);
    }
    if (mode == 1)
    {
        printf("  (Unimodal)\n");
    }
    if (mode == 2)
    {
        printf("  (Bimodal)\n");
    }
    if (mode > 2)
    {
        printf("  (Multimodal)\n");
    }
    printf("----PARTITION VALUES----\n");
    printf("Q1\t\t\t%.3Lf\n", q1);
    printf("Q3\t\t\t%.3Lf\n", q3);
    printf("----MEASURES OF DISPERSION----\n");
    printf("Range\t\t\t%.3Lf\n", range);
    printf("IQR\t\t\t%.3Lf\n", iqr);
    printf("MDv/AM\t\t\t%.3Lf\n", meandevmean);
    printf("MDv/Md\t\t\t%.3Lf\n", meandevmedian);
    printf("Var.\t\t\t%.3Lf\n", var);
    printf("St.Dv.\t\t\t%.3Lf\n", sd);
    printf("CV\t\t\t%.3Lf\n", cv);
    printf("CV %%\t\t\t%.3Lf%%\n", (cv * 100));
    printf("----SKEWNESS & KURTOSIS----\n");
    printf("Beta1\t\t\t%.3Lf", beta1);
    if (beta1 > 0)
    {
        printf("  (Positively Skewed)\n");
    }
    if (beta1 < 0)
    {
        printf("  (Negatively Skewed)\n");
    }
    if (beta1 == 0)
    {
        printf("  (Symmetric)\n");
    }
    printf("Gamma1\t\t\t%.3Lf", gamma1);
    if (gamma1 > 0)
    {
        printf("  (Positively Skewed)\n");
    }
    if (gamma1 < 0)
    {
        printf("  (Negatively Skewed)\n");
    }
    if (gamma1 == 0)
    {
        printf("  (Symmetric)\n");
    }
    printf("Beta2\t\t\t%.3Lf", beta2);
    if (beta2 > 3)
    {
        printf("  (Leptokurtic (peaked))\n");
    }
    if (beta2 < 3)
    {
        printf("  (Platykurtic (flat))\n");
    }
    if (beta2 == 3)
    {
        printf("  (Mesokurtic (normal))\n");
    }
    printf("Gamma2\t\t\t%.3Lf", gamma2);
    if (gamma2 > 0)
    {
        printf("  (Leptokurtic (peaked))\n");
    }
    if (gamma2 < 0)
    {
        printf("  (Platykurtic (flat))\n");
    }
    if (gamma2 == 0)
    {
        printf("  (Mesokurtic (normal))\n");
    }
    //
}

void Grouped(void) {}

void Continuous(void) {}

void Bivariate(void) {}

////////////////////////////////////////////////////
void BisectionMethod1(void)
{
    {
        long double a, b, m, e;
        int k = 0;

        while (1)
        {
            printf("\nEnetr 2 initial approximations\n");
            scanf("%Lf %Lf", &a, &b);

            if ((f(a) * f(b)) < 0.0)
            {
                break;
            }
        }

        printf("\nEnter error tolerance\n");
        scanf("%Lf", &e);

        while (1)
        {
            m = ((a + b) / 2);
            k = k + 1;
            printf("");
            printf("\n%d | %.4Lf | %.4Lf | %.4Lf\n", k, m, a, b);

            if (fabs(f(m)) < e)
            {
                printf("\nHence, root = %.4Lf\n", m);
                break;
            }

            else if ((f(a) * f(m)) < 0.0)
            {
                a = a;
                b = m;
            }

            else
            {
                a = m;
                b = b;
            }
        }
    }
}

void SecantMethod2(void)
{
    long double x0, x1, x2, e;
    int k = 0;

    printf("\nEnter 2 initial approximations\n");
    scanf("%Lf %Lf", &x0, &x1);

    printf("\nEnter error tolerance\n");
    scanf("%Lf", &e);

    while (1)
    {
        x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
        k = k + 1;
        printf("\n%d | %.4Lf | %.4Lf | %.4Lf\n", k, x2, x0, x1);

        if (fabs(f(x2)) < e)
        {
            printf("\nHence, root = %.4Lf\n", x2);
            break;
        }

        x0 = x1;
        x1 = x2;
    }
}

void RegulaFalsi3(void)
{
    long double x1, x2, x3, x_prev = 0, e, root;
    int k = 0;

    printf("\nEnter Initial Approximations :");
    scanf("%Lf %Lf", &x1, &x2);

    if ((f(x1)) * (f(x2)) < 0.0)
    {
        printf("\nThe interval [%.4Lf, %.4Lf] contains a root.\n", x1, x2);
    }
    else
    {
        printf("\nThe interval [%.4Lf, %.4Lf] does not contain a root.\n", x1, x2);
    }

    printf("\nEnter error tolerance:");
    scanf("%Lf", &e);

    printf("\n%d | %.4Lf | %.4Lf | %.4Lf | %.4Lf\n", k, x1, x2, x1, f(x1));

    while (1)
    {
        x3 = x1 - ((f(x1) * (x2 - x1)) / (f(x2) - f(x1)));

        k++;
        printf("\n%d | %.4Lf | %.4Lf | %.4Lf | %.4Lf\n", k, x1, x2, x3, f(x3));

        if (fabs(f(x3)) < e || fabs(x3 - x_prev) < e)
        {
            root = x3;
            break;
        }

        if (f(x1) * f(x3) < 0)
        {
            x2 = x3;
        }
        else
        {
            x1 = x3;
        }

        x_prev = x3;
    }

    printf("\nHence, root = %.4Lf\n", root);
}

void MullerMethod4(void)
{
    double x0, x1, x2, x3, h1, h2, d1, d2, a, b, c, discriminant, e;
    int k = 0;

    printf("\nEnter three initial approximations: ");
    scanf("%lf %lf %lf", &x0, &x1, &x2);

    printf("Enter error tolerance: ");
    scanf("%lf", &e);

    printf("\n Iteration\t  x0\t\t  x1\t\t  x2\t\t  x3\t\t f(x3)\n");
    printf("----------------------------------------------------------------------------------\n");

    while (1)
    {
        h1 = x1 - x0;
        h2 = x2 - x1;
        d1 = (f(x1) - f(x0)) / h1;
        d2 = (f(x2) - f(x1)) / h2;
        a = (d2 - d1) / (h2 + h1);
        b = a * h2 + d2;
        c = f(x2);

        discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
        {
            printf("Discriminant became negative. Complex root encountered.\n");
            break;
        }

        if (fabs(b + sqrt(discriminant)) > fabs(b - sqrt(discriminant)))
            x3 = x2 - (2 * c) / (b + sqrt(discriminant));
        else
            x3 = x2 - (2 * c) / (b - sqrt(discriminant));

        printf("k = %d\t\t%.6lf\t%.6lf\t%.6lf\t%.6lf\t%.6lf\n", k, x0, x1, x2, x3, f(x3));

        if (fabs(f(x3)) < e || fabs(x3 - x2) < e)
        {
            printf("\nApproximated Root is %lf\n", x3);
            break;
        }

        x0 = x1;
        x1 = x2;
        x2 = x3;
        k++;
    }
}

void Chebyshev5(void)
{
    double x0, x1, tolerance;
    int k = 0;

    printf("\nEnter initial approximation: ");
    scanf("%lf", &x0);

    printf("Enter error tolerance: ");
    scanf("%lf", &tolerance);

    printf("\n Iteration\t  x0\t\t  x1\t\t f(x1)\n");
    printf("--------------------------------------------------------------\n");

    while (1)
    {
        // Chebyshev's formula
        x1 = x0 - (f(x0) / df(x0)) * (1 + (f(x0) * ddf(x0)) / (2 * pow(df(x0), 2)));

        printf("k = %d\t\t%.6lf\t%.6lf\t%.6lf\n", k, x0, x1, f(x1));

        if (fabs(f(x1)) < tolerance || fabs(x1 - x0) < tolerance)
        {
            printf("\nApproximated Root is %lf\n", x1);
            break;
        }

        x0 = x1;
        k++;
    }
}

void Bairstow6(void)
{
    float a[10], b[10], c[10], p, q, cc, den, delp, delq;
    int i, n, m, j, k, l;

    printf("Input p, q, degree, number of iterations:\n");
    scanf("%f %f %d %d", &p, &q, &n, &m);

    printf("Input coefficients of polynomial in decreasing order of n:\n");
    for (i = 0; i <= n; i++)
    {
        scanf("%f", &a[i]);
    }

    for (j = 0; j < m; j++)
    {
        b[0] = a[0];
        b[1] = a[1] - p * b[0];
        b[2] = a[2] - p * b[1] - q * b[0];

        for (k = 3; k <= n; k++)
            b[k] = a[k] - p * b[k - 1] - q * b[k - 2];

        int last = n;
        printf("\n b[0] = %f b[1] = %f b[2] = %f b[%d] = %f \n", b[0], b[1], b[2], last, b[last]);

        c[0] = b[0];
        c[1] = b[1] - p * c[0];
        c[2] = b[2] - p * c[1] - q * c[0];

        l = n - 1;
        for (k = 3; k <= l; k++)
            c[k] = b[k] - p * c[k - 1] - q * c[k - 2];

        cc = c[n - 1] - b[n - 1];
        den = c[n - 2] * c[n - 2] - cc * c[n - 3];

        if (fabs(den) < 1e-10)
        {
            printf("Denominator too small. Iteration stopped to avoid division by zero.\n");
            break;
        }

        delp = -(b[n] * c[n - 3] - b[n - 1] * c[n - 2]) / den;
        delq = -(b[n - 1] * cc - b[n] * c[n - 2]) / den;

        printf("\n den = %f \n delp = %f \n delq = %f \n", den, delp, delq);

        p = p + delp;
        q = q + delq;

        printf("ITERATION = %d, P = %f, Q = %f \n", j + 1, p, q);
    }
}

void ForwardSub7(void)
{
    float L[MAX][MAX], b[MAX], x[MAX];
    int n;

    printf("Enter the number of equations: ");
    scanf("%d", &n);

    printf("Enter the lower triangular matrix L (%d x %d):\n", n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            scanf("%f", &L[i][j]);
        }
    }

    printf("Enter the right-hand side vector b:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &b[i]);
    }

    // Display Lower Triangular Matrix
    printf("\nLower Triangular Matrix (L):\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2f ", L[i][j]);
        }
        printf("\n");
    }

    // Display System of Equations
    printf("\nSystem of Equations:\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2fx%d", L[i][j], j + 1);
            if (j < n - 1)
                printf(" + ");
        }
        printf(" = %.2f\n", b[i]);
    }

    // Forward Substitution
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            x[i] -= L[i][j] * x[j];
        }
        x[i] /= L[i][i];
    }

    // Display Solution
    printf("\nSolution vector x:\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] = %.4f\n", i + 1, x[i]);
    }
}

void BackwardSub8(void)
{
    float L[MAX][MAX], b[MAX], x[MAX];
    int n;

    printf("Enter the number of equations: ");
    scanf("%d", &n);

    printf("Enter the lower triangular matrix L (%d x %d):\n", n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            scanf("%f", &L[i][j]);
        }
    }

    printf("Enter the right-hand side vector b:\n");
    for (int i = 0; i < n; i++)
    {
        scanf("%f", &b[i]);
    }

    // Display Lower Triangular Matrix
    printf("\nLower Triangular Matrix (L):\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2f ", L[i][j]);
        }
        printf("\n");
    }

    // Display System of Equations
    printf("\nSystem of Equations:\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.2fx%d", L[i][j], j + 1);
            if (j < n - 1)
                printf(" + ");
        }
        printf(" = %.2f\n", b[i]);
    }

    // Forward Substitution
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i];
        for (int j = 0; j < i; j++)
        {
            x[i] -= L[i][j] * x[j];
        }
        x[i] /= L[i][i];
    }

    // Display Solution
    printf("\nSolution vector x:\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] = %.4f\n", i + 1, x[i]);
    }
}

void Lagrange9(void)
{
    int n;
    printf("Enter the number of data points: ");
    scanf("%d", &n);

    double x_vals[n], y_vals[n];

    // Taking input from the user
    printf("Enter the data points (x[i], f(x[i])):\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d]: ", i);
        scanf("%lf", &x_vals[i]);
        printf("f(x[%d]): ", i);
        scanf("%lf", &y_vals[i]);
    }

    double x;
    printf("\nEnter the value of x to interpolate: ");
    scanf("%lf", &x);

    // Perform Lagrange Interpolation
    double result = 0.0;
    for (int i = 0; i < n; i++)
    {
        double term = y_vals[i];
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j]);
            }
        }
        result += term;
    }

    printf("\nLagrange Interpolated value at x = %.4lf is %.6lf\n", x, result);
}

void NewtonBackward10(void)
{
    int n;
    float x[10], y[10][10] = {0}, xp;

    // Input number of data points
    printf("Enter number of data points: ");
    scanf("%d", &n);

    // Input x values
    printf("\nEnter x values:\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
    }

    // Input f(x) values
    printf("\nEnter f(x) values:\n");
    for (int i = 0; i < n; i++)
    {
        printf("f(x[%d]) = ", i);
        scanf("%f", &y[i][0]);
    }

    // Compute backward difference table
    for (int j = 1; j < n; j++)
    {
        for (int i = n - 1; i >= j; i--)
        {
            y[i][j] = y[i][j - 1] - y[i - 1][j - 1];
        }
    }

    // Print backward difference table
    printf("\nBackward Difference Table:\n");
    printf("--------------------------------------\n");
    printf("|   x   |  f(x)  ");
    for (int j = 1; j < n; j++)
    {
        printf("|  V%d  ", j);
    }
    printf("|\n");
    printf("--------------------------------------\n");

    for (int i = 0; i < n; i++)
    {
        printf("| %.2f | %.2f ", x[i], y[i][0]);
        for (int j = 1; j < n - i; j++)
        {
            printf("| %.2f ", y[i][j]);
        }
        printf("|\n");
    }
    printf("--------------------------------------\n");

    // Input the x value to interpolate
    printf("\nEnter x value to interpolate: ");
    scanf("%f", &xp);

    // Interpolation using Gregory-Newton Backward formula
    float h = x[1] - x[0]; // Assuming equal spacing
    float u = (xp - x[n - 1]) / h;
    float result = y[n - 1][0];
    float uTerm = 1.0;

    for (int j = 1; j < n; j++)
    {
        uTerm *= (u + j - 1) / j;
        result += uTerm * y[n - 1][j];
    }

    // Output the result
    printf("\nInterpolated value at x = %.2f is %.5f\n", xp, result);
}

void NewtonForward11(void)
{
    int n;
    float x[10], y[10][10] = {0}, xp;

    // Input number of data points
    printf("Enter number of data points: ");
    scanf("%d", &n);

    // Input x values
    printf("\nEnter x values:\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] = ", i);
        scanf("%f", &x[i]);
    }

    // Input f(x) values
    printf("\nEnter f(x) values:\n");
    for (int i = 0; i < n; i++)
    {
        printf("f(x[%d]) = ", i);
        scanf("%f", &y[i][0]);
    }

    // Compute forward difference table
    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            y[i][j] = y[i + 1][j - 1] - y[i][j - 1];
        }
    }

    // Print forward difference table
    printf("\nForward Difference Table:\n");
    printf("--------------------------------------\n");
    printf("|   x   |  f(x)  ");
    for (int j = 1; j < n; j++)
    {
        printf("|  Δ^%d  ", j);
    }
    printf("|\n");
    printf("--------------------------------------\n");

    for (int i = 0; i < n; i++)
    {
        printf("| %.2f | %.2f ", x[i], y[i][0]);
        for (int j = 1; j < n - i; j++)
        {
            printf("| %.2f ", y[i][j]);
        }
        printf("|\n");
    }
    printf("--------------------------------------\n");

    // Input the x value to interpolate
    printf("\nEnter x value to interpolate: ");
    scanf("%f", &xp);

    // Gregory-Newton Forward Interpolation
    float h = x[1] - x[0]; // Assuming equally spaced points
    float u = (xp - x[0]) / h;
    float result = y[0][0];
    float uTerm = 1.0;

    for (int j = 1; j < n; j++)
    {
        uTerm *= (u - j + 1) / j;
        result += uTerm * y[0][j];
    }

    // Output the result
    printf("\nInterpolated value at x = %.2f is %.5f\n", xp, result);
}

void Hermite12(void)
{
    int n;
    printf("Enter the number of data points: ");
    scanf("%d", &n);

    double x_val[n], y_val[n], dy_val[n];

    // Input values
    printf("Enter the data points (x[i], f(x[i]), f'(x[i])):\n");
    for (int i = 0; i < n; i++)
    {
        printf("x[%d] : ", i);
        scanf("%lf", &x_val[i]);
        printf("f(x[%d]) : ", i);
        scanf("%lf", &y_val[i]);
        printf("f'(x[%d]) : ", i);
        scanf("%lf", &dy_val[i]);
    }

    // Display Hermite Table
    printf("\n----------------------------------\n");
    printf("|       HERMITE INTERPOLATION     |\n");
    printf("----------------------------------\n");
    printf("|  x[i]   |  f(x[i])  |  f'(x[i]) |\n");
    printf("----------------------------------\n");
    for (int i = 0; i < n; i++)
    {
        printf("|  %.2lf   |   %.2lf    |   %.2lf   |\n", x_val[i], y_val[i], dy_val[i]);
    }
    printf("----------------------------------\n");

    // Interpolation point
    double x;
    printf("\nEnter the value of x to interpolate: ");
    scanf("%lf", &x);

    // Hermite Interpolation Logic
    double result = 0.0;

    for (int i = 0; i < n; i++)
    {
        // Compute l_i(x)
        double l_i = 1.0;
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                l_i *= (x - x_val[j]) / (x_val[i] - x_val[j]);
            }
        }

        // Compute l_i'(x_i)
        double l_i_prime = 0.0;
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                double term = 1.0 / (x_val[i] - x_val[j]);
                for (int k = 0; k < n; k++)
                {
                    if (k != i && k != j)
                    {
                        term *= (x_val[i] - x_val[k]) / (x_val[i] - x_val[k]);
                    }
                }
                l_i_prime += term;
            }
        }

        // Compute A_i and B_i
        double A_i = (1 - 2 * (x - x_val[i]) * l_i_prime) * (l_i * l_i);
        double B_i = (x - x_val[i]) * (l_i * l_i);

        // Accumulate result
        result += A_i * y_val[i] + B_i * dy_val[i];
    }

    // Display interpolated result
    printf("\nHermite Interpolated value at x = %.2lf is %.5lf\n", x, result);
}

void Jacobi13(void)
{
    float A[N][N] = {
        {4, 1, 1}, // Equation 1: 4x1 + x2 + x3 = 2
        {1, 5, 2}, // Equation 2: x1 + 5x2 + 2x3 = -6
        {1, 2, 3}  // Equation 3: x1 + 2x2 + 3x3 = -4
    };

    float b[N] = {2, -6, -4};       // Right-hand side constants
    float x[N] = {0.5, -0.5, -0.5}; // Initial guess
    float x_new[N] = {0};           // To store new iteration values

    printf("\nJacobi Iteration Method:\n");
    printf("----------------------------------------------------");
    printf("\n%-10s %-10s %-10s %-10s\n", "Iteration", "x[1]", "x[2]", "x[3]");
    printf("----------------------------------------------------\n");

    for (int iter = 1; iter <= MAX; iter++)
    {
        for (int i = 0; i < N; i++)
        {
            float sum = 0;
            for (int j = 0; j < N; j++)
            {
                if (j != i)
                    sum += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }

        // Update x values
        for (int i = 0; i < N; i++)
            x[i] = x_new[i];

        // Print iteration result
        printf("%-10d %-10.4f %-10.4f %-10.4f\n", iter, x[0], x[1], x[2]);
    }
}

void BirgeVieta14(void)
{
    int degree, i, step = 1;
    float coeff[MAX], b[MAX], c[MAX], x, px, px_deriv, new_x;

    printf("Enter degree of polynomial: ");
    scanf("%d", &degree);

    printf("Enter coefficients (from highest degree to constant):\n");
    for (i = 0; i <= degree; i++)
    {
        printf("Coefficient of x^%d: ", degree - i);
        scanf("%f", &coeff[i]);
    }

    printf("Enter initial guess: ");
    scanf("%f", &x);

    while (1)
    {
        // Compute b[] using synthetic division
        b[0] = coeff[0];
        for (i = 1; i <= degree; i++)
        {
            b[i] = coeff[i] + x * b[i - 1];
        }

        // Compute c[] (derivative)
        c[0] = b[0];
        for (i = 1; i < degree; i++)
        {
            c[i] = b[i] + x * c[i - 1];
        }

        // Evaluate P(x) using Horner's method
        px = coeff[0];
        for (i = 1; i <= degree; i++)
        {
            px = px * x + coeff[i];
        }

        px_deriv = c[degree - 1];

        printf("Step %d: x = %.6f, P(x) = %.6f, P'(x) = %.6f\n", step, x, px, px_deriv);

        if (fabs(px_deriv) < EPSILON)
        {
            printf("Derivative too small. Stopping.\n");
            break;
        }

        new_x = x - px / px_deriv;

        if (fabs(new_x - x) < EPSILON)
        {
            break;
        }

        x = new_x;
        step++;
    }

    printf("\nRoot found: %.6f\n", new_x);

    // Deflate polynomial
    int deflated_degree = degree - 1;
    printf("\nDeflated Polynomial:\n");
    for (i = 0; i <= deflated_degree; i++)
    {
        printf("Coefficient of x^%d = %.6f\n", deflated_degree - i, b[i]);
    }
}

void Trapezoidal15(void)
{
    float a, b, result;

    printf("f(x) = 1 / (1 + x)\n");
    printf("Enter the lower limit (a) and upper limit (b): ");
    scanf("%f%f", &a, &b);

    result = ((g(a) + g(b)) / 2.0) * (b - a);

    printf("\nValue of the integral using Trapezoidal Rule: %f\n", result);
}

void Simpson13rd16(void)
{
    float a, b, result, m;

    printf("f(x) = 1 / (1 + x)\n");
    printf("Enter the lower limit (a) and upper limit (b): ");
    scanf("%f%f", &a, &b);

    m = (a + b) / 2.0; // Midpoint
    result = (g(a) + 4 * g(m) + g(b)) * (b - a) / 6.0;

    printf("\nValue of the integral using Simpson's 1/3 Rule: %f\n", result);
}

void Simpson38th17(void)
{
    float a, b, result, x1, x2;

    printf("f(x) = 1 / (1 + x)\n");
    printf("Enter the lower limit (a) and upper limit (b): ");
    scanf("%f%f", &a, &b);

    x1 = (2 * a + b) / 3.0;
    x2 = (a + 2 * b) / 3.0;

    result = (g(a) + 3 * g(x1) + 3 * g(x2) + g(b)) * (b - a) / 8.0;

    printf("\nValue of the integral using Simpson's 3/8 Rule: %f\n", result);
}

//
int Mathematics(void)
{
    int choose = 0;

    while (1)
    {
        printf("\n-----Mathematics menu-----\n");
        printf("\n1> Addition & Subtraction");
        printf("\n2> Multiplication");
        printf("\n3> Division");
        printf("\n4> Exponent");
        printf("\n5> Factorial");
        printf("\n6> Permutation");
        printf("\n7> Combination");
        printf("\n8> Create an AP");
        printf("\n9> Create a GP");
        printf("\n10> Create an HP");
        printf("\n11> Simple interest");
        printf("\n12> Compound interest");
        printf("\n13> Quadratic equation roots");
        printf("\n14> HCF");
        printf("\n15> LCM");
        printf("\n16> Matrix Addition");
        printf("\n17> Matrix Subtraction");
        printf("\n18> Matrix Transpose");
        printf("\n19> Matrix Multiplication");
        printf("\n20> Exit");
        printf("\nEnter your choice\n");
        scanf("%d", &choose);

        switch (choose)
        {
        case 1:
            AdditionSubtraction1();
            break;
        case 2:
            Multiplication2();
            break;
        case 3:
            Divide2Nos3();
            break;
        case 4:
            Exponent4();
            break;
        case 5:
            Factorial5();
            break;
        case 6:
            Permutation6();
            break;
        case 7:
            Combination7();
            break;
        case 8:
            AP8();
            break;
        case 9:
            GP9();
            break;
        case 10:
            HP10();
            break;
        case 11:
            SimpleInterest11();
            break;
        case 12:
            CompoundInterest12();
            break;
        case 13:
            QuadraticRoots13();
            break;
        case 14:
            HCF14();
            break;
        case 15:
            LCM15();
            break;
        case 16:
            MatrixAdd16();
            break;
        case 17:
            MatrixSubtract17();
            break;
        case 18:
            MatrixTranspose18();
            break;
        case 19:
            MatrixMulti19();
            break;
        case 20:
            printf("\nExited Mathematics section\n");
            return 0;
        default:
            printf("\nInvalid choice, please choose a valid option\n");
            break;
        }
    }
}
//

//
int DescStats(void)
{
    int choose = 0;

    while (1)
    {
        printf("\n-----Descriptive Statistics menu-----\n");
        printf("\n1] Raw Data");
        printf("\n2] Grouped Data");
        printf("\n3] Continuous Data");
        printf("\n4] Bivariate Data");
        printf("\n5] Exit");
        printf("\nEnter your choice\n");
        scanf("%d", &choose);

        switch (choose)
        {
        case 1:
            Raw();
            break;
        case 2:
            Grouped();
            break;
        case 3:
            Continuous();
            break;
        case 4:
            Bivariate();
            break;
        case 5:
            printf("\nExited Desc.Stats. section\n");
            return 0;
        default:
            printf("\nInvalid choice, please choose a valid option!!\n");
            break;
        }
    }
}
//

//
int Numerical(void)
{
    int choose = 0;

    while (1)
    {
        printf("\n-----Numerical menu-----\n");
        printf("\n1} Bisection method");
        printf("\n2} Secant method");
        printf("\n3} Regula Falsi method");
        printf("\n4} Muller method");
        printf("\n5} Chebyshev method");
        printf("\n6} Bairstow method");
        printf("\n7} Forward Subs. method");
        printf("\n8} Backwad Subs. method");
        printf("\n9} Lagrange Interp. method");
        printf("\n10} Newton Backward Interp. method");
        printf("\n11} Newton Forward Interp. method");
        printf("\n12} Hermite Interp. method");
        printf("\n13} Jacobi Iter. method");
        printf("\n14} Birge-Vieta method");
        printf("\n15} Trapezoidal int. method");
        printf("\n16} Simpson's 1/3rd int. method");
        printf("\n17} Simpson's 3/8th int. method");
        printf("\n18} Exit");
        printf("\nEnter your choice\n");
        scanf("%d", &choose);

        switch (choose)
        {
        case 1:
            BisectionMethod1();
            break;
        case 2:
            SecantMethod2();
            break;
        case 3:
            RegulaFalsi3();
            break;
        case 4:
            MullerMethod4();
            break;
        case 5:
            Chebyshev5();
            break;
        case 6:
            Bairstow6();
            break;
        case 7:
            ForwardSub7();
            break;
        case 8:
            BackwardSub8();
            break;
        case 9:
            Lagrange9();
            break;
        case 10:
            NewtonBackward10();
            break;
        case 11:
            NewtonForward11();
            break;
        case 12:
            Hermite12();
            break;
        case 13:
            Jacobi13();
            break;
        case 14:
            BirgeVieta14();
            break;
        case 15:
            Trapezoidal15();
            break;
        case 16:
            Simpson13rd16();
            break;
        case 17:
            Simpson38th17();
            break;
        case 18:
            printf("\nExited Numerical section\n");
            return 0;
        default:
            printf("\nInvalid choice, please choose a valid option!!\n");
            break;
        }
    }
}
//

//
int main(void)
{
    int choose = 0;

    while (1)
    {
        printf("\n-----Calculator menu-----\n");
        printf("\n1) Mathematics");
        printf("\n2) Desc. Stats.");
        printf("\n3) Numerical");
        printf("\n4) Exit");
        printf("\nEnter your choice\n");
        scanf("%d", &choose);

        switch (choose)
        {
        case 1:
            Mathematics();
            break; //
        case 2:
            DescStats();
            break; //
        case 3:
            Numerical();
            break; //
        case 4:
            printf("\nExited Calculator\n");
            return 0;
        default:
            printf("\nInvalid choice, please choose a valid option!!\n");
            break;
        }
    }

    return 0;
}
//