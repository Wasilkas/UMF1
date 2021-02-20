#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <locale.h>

using namespace std;

struct matrix
{
	int n;   // размерность матрицы
	int offset;   // смещение
	vector<double> d_0;   // главна€ диагональ
	vector<double> du_1;   // верхн€€ диагональ со смещением 1
	vector<double> du_offset;   // верхн€€ диагональ со смещением offset
	vector<double> dl_1;	// нижн€€ диагональ со смещением 1
	vector<double> dl_offset;	// нижн€€ диагональ со смещением offset

	matrix(int n, int offset)
	{
		this->n = n;
		this->offset = offset;
		d_0.resize(n);
		du_1.resize(n-1);
		du_offset.resize(n-offset);
		dl_1.resize(n-1);
		dl_offset.resize(n-offset);
	}
};

double CountFirstStep(double begin, double end, int num_steps, double q);
void ReadGrid(string name_file, vector<double>* grid);
void FillingKnotsIN(double lambda, double gamma, vector<double> x, vector<double> y, matrix* A, double* b);
vector <int> SearchBoundKnots(int x_num_all_knots, int y_num_all_knots, int num_all_knots);
void AccountBound(vector<double> x, vector<double> y, matrix* A, double* b, vector <int> bound_knots);
double Skal(int i, matrix* mat, double* xk); // —кал€рное произведение i-ой строки матрицы с вектором xk
double Err(matrix* mat, double* f, double* xk, double f_norm);
int Gausse_Zeidel(matrix* mat, double* f, double* xk, double eps, int maxiter);

double count_bound_func(double X, double Y)
{
	return X * X * X * X;
}

double count_func(double X, double Y)
{
	return 2. * X * X * X * X - 12. * X * X;
}

void write(matrix* A)
{
	
	for (int i = 0; i < A->n; i++)
	{
		for (int j = 0; j < A->n; j++)
		{
			if (i == j)
				cout << A->d_0[i] << "  ";
			else
				if (j == i + 1)
					cout << A->du_1[i] << "  ";
				else
					if (j == i + A->offset)
						cout << A->du_offset[i] << "  ";
					else
						if (j == i - 1)
							cout << A->dl_1[j] << "  ";
						else
							if (j == i - A->offset)
								cout << A->dl_offset[j] << "  ";
							else
								cout << "0  ";
		}

		cout << endl;
	}
}

void main()
{
	setlocale(LC_ALL, "rus");

	vector<double> x_grid;
	ReadGrid("x_input.txt", &x_grid);

	vector<double> y_grid;
	ReadGrid("y_input.txt", &y_grid);

	ifstream coeff_input("coeff_input.txt");
	double lambda, gamma;
	coeff_input >> lambda >> gamma;
	coeff_input.close();

	matrix* A = new matrix(x_grid.size() * y_grid.size(), x_grid.size());
	double* b = new double[A->n];
	FillingKnotsIN(lambda, gamma, x_grid, y_grid, A, b);

	vector <int> bound_knots = SearchBoundKnots(x_grid.size(), y_grid.size(), A->n);
	AccountBound(x_grid, y_grid, A, b, bound_knots);

	//write(A);

	double* q = new double[A->n];
	Gausse_Zeidel(A, b, q, 1e-10, 1e+10);

	cout << "¬ектор правой части : ";

	//for (int i = 0; i < A->n;i++)
	//{
	//	cout << i << " : " << b[i] << endl;
	//}

	cout << "–езультат : \n";
	for (int i = 0; i < A->n; i++)
	{
		int X_ind = i % x_grid.size();
		int Y_ind = i / x_grid.size();
		double resh = count_bound_func(x_grid[X_ind], y_grid[Y_ind]);
		cout << i << " : " << q[i] << " | " << resh << " | " << abs(q[i] - resh) << endl;
	}

}


double CountFirstStep(double begin, double end, int num_steps, double q)
{
	double h0;

	h0 = (q != 1) ? (end - begin) * (1. - q) / (1. - pow(q, num_steps)) : (end - begin) / num_steps;

	return h0;
}

void ReadGrid(string name_file, vector<double>* grid)
{
	ifstream in(name_file);

	double begin, end, num_steps, q;

	in >> begin >> end;

	in >> num_steps;

	in >> q;

	in.close();

	double h0 = CountFirstStep(begin, end, num_steps, q);

	double coord = begin;
	double h = h0;
	while (coord < end)
	{
		grid->push_back(coord);
		coord += h;
		h *= q;
	}

	grid->push_back(end);
}

void FillingKnotsIN(double lambda, double gamma, vector<double> x, vector<double> y, matrix* A, double* b)
{

	int x_num_knots = x.size();
	int min_index = x_num_knots;
	int max_index = A->n - min_index;

	// ѕроходим по узлам, не затрагива€ нижнюю и верхнюю грани
	for (int i = min_index; i < max_index; i++)
	{
		int X_ind = i % x_num_knots;	// определение локального узла на оси х по глобальному

		if (X_ind && X_ind != x_num_knots - 1)
		{
			int Y_ind = i / x_num_knots;	// определение локального узла на оси y по глобальному
			double h_x_i = x[X_ind + 1] - x[X_ind];	// шаг справа от узла
			double h_x_i_1 = x[X_ind] - x[X_ind - 1];	// шаг слева от узла

			double h_y_j = y[Y_ind + 1] - y[Y_ind];	// шаг сверху от узла
			double h_y_j_1 = y[Y_ind] - y[Y_ind - 1];	// шаг снизу от узла


			// 1 / (h_x_i_1 * h_x_i) + 1 / (h_y_j_1 * h_y_j) = (h_y_j_1 * h_y_j + h_x_i_1 * h_x_i) / (h_y_j_1 * h_y_j * h_x_i_1 * h_x_i)
			A->d_0[i] = lambda * 2. * ((h_y_j_1 * h_y_j + h_x_i_1 * h_x_i) / (h_y_j_1 * h_y_j * h_x_i_1 * h_x_i)) + gamma;
			A->dl_1[i - 1] = -lambda * 2. / (h_x_i_1 * (h_x_i + h_x_i_1));
			A->du_1[i] = -lambda * 2. / (h_x_i * (h_x_i + h_x_i_1));
			A->dl_offset[i - A->offset] = -lambda * 2. / (h_y_j_1 * (h_y_j + h_y_j_1));
			A->du_offset[i] = -lambda * 2. / (h_y_j * (h_y_j + h_y_j_1));
			b[i] = count_func(x[X_ind], y[Y_ind]);
		}
	}
}

vector <int> SearchBoundKnots(int x_num_all_knots, int y_num_all_knots, int num_all_knots)
{
	vector <int> bound_knots;

	int X_ind_left_bound = 0;
	int X_ind_right_bound = x_num_all_knots - 1;

	int Y_ind_bottom_bound = 0;
	int Y_ind_top_bound = y_num_all_knots - 1;

	for (int i = 0; i < num_all_knots; i++)
	{
		int X_ind = i % x_num_all_knots;
		int Y_ind = i / x_num_all_knots;

		if (Y_ind == Y_ind_bottom_bound || Y_ind == Y_ind_top_bound)
		{
			bound_knots.push_back(i);
		}
		else
		{
			if (X_ind == X_ind_left_bound || X_ind == X_ind_right_bound)
			{
				bound_knots.push_back(i);
			}
		}
	}

	return bound_knots;

}

void AccountBound(vector<double> x, vector<double> y, matrix* A, double* b, vector <int> bound_knots)
{
	int all_num_bound_knots = bound_knots.size();
	for (int i = 0; i < all_num_bound_knots; i++)
	{
		int ind = bound_knots[i];

		int X_ind = ind % x.size();
		int Y_ind = ind / x.size();
		b[ind] = count_bound_func(x[X_ind], y[Y_ind]);

		A->d_0[ind] = 1.;
		if (ind < A->n - 1)
		{
			A->du_1[ind] = 0.;
			if (ind > 0) A->dl_1[ind - 1] = 0.;
		}

		if (ind < A->n - A->offset)
		{
			A->du_offset[ind] = 0.;
			if (ind > A->offset - 1) A->dl_offset[ind - A->offset] = 0.;
		}



	}
}

double Skal(int i, matrix* mat, double* xk) // —кал€рное произведение i-ой строки матрицы с вектором xk
{
	double sum = mat->d_0[i] * xk[i];
	int n = mat->n;
	int m = mat->offset - 1;
	if (i > 0) sum += mat->dl_1[i - 1] * xk[i - 1];
	if (i < n - 1) sum += mat->du_1[i] * xk[i + 1];

	int ind;
	ind = i - m - 1;
	if (ind > -1)
		sum += mat->dl_offset[ind] * xk[ind];

	ind = i + m + 1;
	if (ind < n)
		sum += mat->du_offset[i] * xk[ind];

	return sum;
}

double Err(matrix* mat, double* f, double* xk, double f_norm)
{
	int n = mat->n;
	double sum = 0;
	double buf;
	for (int i = 0; i < n; i++)
	{
		buf = f[i] - Skal(i, mat, xk);
		sum += buf * buf;
	}

	return sqrt(sum) / f_norm;
}

int Gausse_Zeidel(matrix* mat, double* f, double* xk, double eps, int maxiter)
{
	cout << "Gausse Zeidel" << endl;
	int n = mat->n;
	int k; // —чЄтчик итераций
	int i;
	vector<double> d0 = mat->d_0;
	double f_norm = 0;
	for (i = 0; i < n; i++)
		f_norm += f[i] * f[i];
	f_norm = sqrt(f_norm);

	for (i = 0; i < n; i++)
		xk[i] = 0;

	double W = 1.;
	double err_k = Err(mat, f, xk, f_norm);
	for (k = 0; k < maxiter && err_k > eps; k++)
	{
		for (i = 0; i < n; i++)
			xk[i] += W / d0[i] * (f[i] - Skal(i, mat, xk));

		err_k = Err(mat, f, xk, f_norm);
		cout << "Iteration Number " << k << " Error: " << err_k << endl;
	}
	if (k == maxiter)
		cout << "Limit of iterations is overrided";

	return 0;
}