using System.Text;

int n = 4;
double E = 10e-3;
bool helper = true;
void PrintMatr<T>(T[,] matr, T[] b)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            Console.Write(matr[i, j] + " ");
            if (j == 3 && i != 1 && helper)
                Console.Write(" ");
        }
        Console.WriteLine("║" + b[i]);
    }
}

void PrintRes<T>(T[] res)
{
    Console.WriteLine("Result:");
    for(int i = 0; i < res.Length; ++i)
    {
        string mess = "x" + (i + 1).ToString();
        Console.WriteLine($"{mess} = {res[i]}");
    }
}

double diffArr(double[] left, double[] right)
{
    double answ = 0;
    for(int i = 0; i < left.Length; ++i)
    {
        answ = Math.Abs(left[i] - right[i]) > answ ? Math.Abs(left[i] - right[i]) : answ;
    }
    return answ;
}

double[] GaussMeth(double[,] matr, double[] b)
{
    double d, s;
    var res = new double[n];
    for (int k = 0; k < n; ++k)
    {
        double max_el = double.MinValue;
        int index = 0;
        for (int i = k; i < n; ++i)
        {
            if (max_el < Math.Abs(matr[i, k]))
            {
                max_el = Math.Abs(matr[i, k]);
                index = i;
            }
        }
        if(index != k)
        {
            for(int i = 0; i < n; ++i)
            {
                var temp = matr[index, i];
                matr[index, i] = matr[k, i];
                matr[k, i] = temp;
            }
            var t = b[index];
            b[index] = b[k];
            b[k] = t;
        }
        for (int j = k + 1; j < n; ++j)
        {
            
            d = matr[j, k] / matr[k, k];
            for (int i = k; i < n; ++i)
            {
                matr[j, i] = matr[j, i] - d * matr[k, i];
            }
            b[j] = b[j] - d * b[k];
        }
    }
    for (int k = n - 1; k >= 0; --k)
    {
        d = 0;
        for (int j = k; j < n; ++j)
        {
            s = matr[k, j] * (k == n - 1 ? 0 : res[j]);
            d += s;
        }

        res[k] = Math.Round(((b[k] - d) / matr[k, k]), 3);
    }


    return res;
}

void Gauss()
{
    var matr = new double[,] { { 4, 3, 1, 0 }, { -2, 2, 6, 1 }, { 0, 5, 2, 3 }, { 0, 1, 2, 7 } };
    var b = new double[] { 14, 31, 33, 45 };
    Console.WriteLine("Метод Гауса для матриці:");
    PrintMatr(matr, b);

    var res = GaussMeth(matr, b);

    PrintRes(res);
    Det();
    Ober();
    Console.WriteLine();
    helper = false;
    n--;
}

void Det()
{
    var matr = new double[,] { { 4, 3, 1, 0 }, { -2, 2, 6, 1 }, { 0, 5, 2, 3 }, { 0, 1, 2, 7 } };
    double det = 1;
    for (int i = 0; i < n; ++i)
    {
        int k = i;
        for (int j = i + 1; j < n; ++j)
        {
            if (Math.Abs(matr[j, i]) > Math.Abs(matr[k, i]))
                k = j;
        }
        for(int j = 0; j < n; ++j)
        {
            var temp = matr[i, j];
            matr[i, j] = matr[k, j];
            matr[k, j] = temp;
        }
        if (i != k)
            det *= -1;
        det *= matr[i, i];

        for(int j = i+1; j < n; ++j)
        {
            matr[i, j] /= matr[i, i];
        }
        for (int j = 0; j < n; ++j)
        {
            if(j != i && (Math.Abs(matr[j, i] ) > E))
            {
                for(k = i+1; k < n; ++k)
                    matr[j, k] -= matr[i, k] * matr[j, i];
            }
        }
    }
    Console.WriteLine($"Визначник = {Math.Round(det, 3)}");
}
void Ober()
{
    var answ = new double[n][];
    for (int i = 0; i < n; i++)
    {
        var matr = new double[,] { { 4, 3, 1, 0 }, { -2, 2, 6, 1 }, { 0, 5, 2, 3 }, { 0, 1, 2, 7 } };
        var b = new double[] { 0, 0, 0, 0 };
        b[i] = 1;
        answ[i] = GaussMeth(matr, b);
    }
    Console.WriteLine("Обернена матриця:");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; ++j)
        {
            Console.Write(answ[j][i] + " ");
        }
        Console.WriteLine();  
     }
            

}

void TridiagonalMatr()
{
    var matr = new int[,] { { 1, 2, 0 }, { 2, 2, 4 }, { 0, 4, 3 } };
    var b = new int[] { 11, 34, 31 };
    var res = new double[3];
    Console.WriteLine("Метод прогонки для матриці:");
    PrintMatr(matr, b);
    double a1 = matr[0, 1] / matr[0, 0] * (-1);
    double b1 = b[0] / matr[0, 0];
    double a2 = matr[1, 2] / (matr[1, 1] * (-1) - a1 * matr[1, 0]);
    double b2 = (b1 * matr[1, 0] - b[1]) / (matr[1, 1] * (-1) - a1 * matr[1, 0]);
    double b3 = (b2 * matr[2, 1] - b[2]) / (matr[2, 2] * (-1) - a2 * matr[2, 1]);
    res[2] = b3;
    res[1] = a2 * res[2] + b2;
    res[0] = a1 * res[1] + b1;
    PrintRes(res);
    Console.WriteLine();
    n++;
}

void Zeidel()
{
    var matr = new double[,] { { 4, 0, 1, 1 }, { 0, 3, 0, 1 }, { 1, 0, 2, 0 }, { 1, 1, 0, 5 } };
    var b = new double[] { 11, 10, 7, 23 };
    Console.WriteLine("Метод Зейделя для матриці:");
    PrintMatr(matr, b);
    var x_curr = new double[4];
    var x_prev = new double[4];
    uint count = 0;
    do
    {
        Array.Copy(x_curr, x_prev, x_curr.Length);
        for(int i = 0; i < n; ++i)
        {
            x_curr[i] = b[i];
            for (int j = 0; j < n; ++j)
                if (j != i)
                    x_curr[i] -= x_prev[j] * matr[i, j];
            x_curr[i] /= matr[i, i];
        }

        ++count;
    } while (diffArr(x_curr, x_prev) > E);
    PrintRes(x_curr);
    Console.Write($"Кількість ітерацій: {count}");
}



Console.OutputEncoding = Encoding.UTF8;

Gauss();
TridiagonalMatr();
Zeidel();