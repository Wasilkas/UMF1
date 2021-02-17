#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <locale.h>

using namespace std;
typedef double real;

struct matrix
{
	int n;   // размерность матрицы
	int offset;   // смещение
	real* d_0;   // главная диагональ
	real* du_1;   // верхняя диагональ со смещением 1
	real* du_offset;   // верхняя диагональ со смещением offset
	real* dl_1;	// нижняя диагональ со смещением 1
	real* dl_offset;	// нижняя диагональ со смещением offset

	matrix(int n, int offset)
	{
		this->n = n;
		this->offset = offset;
		d_0 = new real[n]{};
		du_1 = new real[n - 1]{};
		du_offset = new real[n - offset]{};
		dl_1 = new real[n - 1]{};
		dl_offset = new real[n - offset]{};
	}
};

struct axis
{
	int num_knots;	// количество узлов на оси
	real main_knots[3]{};	// узлы, формирующие сетку
	real* all_knots;	// все узлы сетки

	axis(real* main_knots, int* num_div, real q)
	{
		for (int i = 0; i < 3; i++)
		{
			this->main_knots[i] = main_knots[i];
		}

		num_knots = num_div[0] + num_div[1] + 3;

		all_knots = new real[num_knots];
		all_knots[0] = main_knots[0];
		real step = (main_knots[2] - main_knots[0]) / (num_knots - 1);

		if (q == 1)
		{
			for (int i = 1; i < num_knots - 1; i++)
			{
				all_knots[i] = all_knots[i - 1] + step;
			}
			cout << "Равномерная сетка" << endl;
		}
		else
		{
			for (int i = 1; i < num_knots - 1; i++)
			{
				//all_knots[i] = all_knots[i - 1] + step;

			}
			cout << "НеРавномерная сетка" << endl;
		}

		all_knots[num_knots - 1] = main_knots[2];
	}
};


void write_A(matrix* A)
{
	cout << "/-------------------" << endl;
	cout << "Матрица А" << endl;
	cout << "Главная диагональ" << endl;

	for (int i = 0; i < A->n; i++)
	{
		cout << i << " : " << A->d_0[i] << endl;
	}

	cout << "/-------------------------------" << endl;
	cout << "Диагональ 1 (верхняя)" << endl;

	for (int i = 0; i < A->n-1; i++)
	{
		cout << i << " : " << A->du_1[i] << endl;
	}

	cout << "/-------------------------------" << endl;
	cout << "Диагональ 1 (нижняя)" << endl;

	for (int i = 0; i < A->n-1; i++)
	{
		cout << i << " : " << A->dl_1[i] << endl;
	}

	cout << "/-------------------------------" << endl;
	cout << "Диагональ " << A->offset << " (верхняя)" << endl;

	for (int i = 0; i < A->n-A->offset; i++)
	{
		cout << i << " : " << A->du_offset[i] << endl;
	}

	cout << "/-------------------------------" << endl;
	cout << "Диагональ " << A->offset << " (нижняя)" << endl;

	for (int i = 0; i < A->n - A->offset; i++)
	{
		cout << i << " : " << A->dl_offset[i] << endl;
	}
}

axis* read_axis(ifstream& in)
{
	real x[3];
	for (int i = 0; i < 3; i++)
		in >> x[i];

	int xn[2];
	for (int i = 0; i < 2; i++)
		in >> xn[i];

	real q;
	in >> q;

	axis* out = new axis(x, xn, q);

	return out;
}

double* count_steps(axis* ax)
{
	int num_steps = ax->num_knots - 1;
	double* steps = new double[num_steps];

	for (int i = 0; i < num_steps; i++)
	{
		steps[i] = ax->all_knots[i + 1] - ax->all_knots[i];
	}

	return steps;
}

matrix* build_A(axis* X_axis, axis* Y_axis, real lambda, real gamma)
{
	int size_matrix = X_axis->num_knots * Y_axis->num_knots;

	matrix* A = new matrix(size_matrix, X_axis->num_knots);

	double* X_steps = count_steps(X_axis);
	double* Y_steps = count_steps(Y_axis);

	int min_index = A->offset;
	int max_index = size_matrix - min_index;


	for (int i = min_index; i < max_index; i++)
	{
		int X_ind = i % X_axis->num_knots;
		int Y_ind = i / X_axis->num_knots;
		real h_x_i_1 = X_steps[X_ind - 1];
		real h_x_i = X_steps[X_ind];
		real h_y_j_1 = Y_steps[Y_ind - 1];
		real h_y_j = Y_steps[Y_ind];
		A->d_0[i] = lambda * 2 * (1 / (h_x_i_1 * h_x_i) + 1 / (h_y_j_1 * h_y_j)) + gamma;
		A->dl_1[i - 1] = -lambda * 2 / h_x_i_1 / (h_x_i + h_x_i_1) + gamma;
		A->du_1[i] = -lambda * 2 / h_x_i / (h_x_i + h_x_i_1) + gamma;
		A->dl_offset[i - A->offset] = -lambda * 2 / h_y_j_1 / (h_y_j + h_y_j_1) + gamma;
		A->du_offset[i] = -lambda * 2 / h_y_j / (h_y_j + h_y_j_1) + gamma;
	}

	return A;
}

double count_func(double X, double Y)
{
	return 2;
}

double* build_b(axis* X_axis, axis* Y_axis)
{
	int X_num_knots = X_axis->num_knots;
	int Y_num_knots = Y_axis->num_knots;
	real* b = new real[X_num_knots * Y_num_knots];
	for (int i = 0; i < Y_num_knots; i++)
	{
		for (int j = 0; j < X_num_knots; j++)
		{
			int knot_num = i * X_num_knots + j;
			b[knot_num] = count_func(X_axis->all_knots[j], Y_axis->all_knots[i]);
		}
	}

	return b;
}

vector <int> count_bound_knots(axis* X_axis, axis* Y_axis)
{
	vector <int> bound_knots;
	int all_num_knots = X_axis->num_knots * Y_axis->num_knots;

	real X_left_bound = X_axis->main_knots[0];
	real X_center_bound = X_axis->main_knots[1];
	real X_right_bound = X_axis->main_knots[2];

	real Y_bottom_bound = Y_axis->main_knots[0];
	real Y_center_bound = Y_axis->main_knots[1];
	real Y_top_bound = Y_axis->main_knots[2];

	for (int i = 0; i < all_num_knots; i++)
	{
		int X_ind = i % X_axis->num_knots;
		int Y_ind = i / X_axis->num_knots;

		real Y_coord = Y_axis->all_knots[Y_ind];
		real X_coord = X_axis->all_knots[X_ind];



		if (Y_coord < Y_center_bound)
		{
			if (Y_coord == Y_bottom_bound)
			{
				bound_knots.push_back(i);
			}
			else
				if (X_coord == X_left_bound || X_coord == X_right_bound)
					bound_knots.push_back(i);
		}
		else
			if (Y_coord == Y_center_bound)
			{
				if (X_coord == X_left_bound || X_coord > X_center_bound)
					bound_knots.push_back(i);
			}
			else
			{
				if (Y_coord == Y_top_bound)
				{
					if (X_coord <= X_center_bound)
						bound_knots.push_back(i);
				}
				else
					if (X_coord == X_left_bound || X_coord == X_center_bound)
						bound_knots.push_back(i);
			}


	}

	return bound_knots;

}

real count_bound_func(real X, real Y)
{
	return 1;
}

void account_bound(axis* X_axis, axis* Y_axis, matrix* A, real* b, vector <int> bound_knots)
{
	int all_num_bound_knots = bound_knots.size();

	for (int i = 0; i < all_num_bound_knots; i++)
	{
		int ind = bound_knots[i];
		int X_ind = ind % X_axis->num_knots;
		real X_coord = X_axis->all_knots[X_ind];

		int Y_ind = ind / Y_axis->num_knots;
		real Y_coord = Y_axis->all_knots[Y_ind];

		b[ind] = count_bound_func(X_coord, Y_coord);

		A->dl_1[ind - 1] = 0;
		A->dl_offset[ind - A->offset] = 0;
		A->du_1[ind] = 0;
		A->du_offset[ind] = 0;
		A->d_0[ind] = 1;
	}
}

vector <int> count_fictitious_knots(axis* X_axis, axis* Y_axis)
{
	vector <int> fictitious_knots;

	real Y_center_coord = Y_axis->main_knots[1];
	real X_center_coord = X_axis->main_knots[1];

	real Y_coord = Y_axis->all_knots[Y_axis->num_knots - 1];

	for (int i = X_axis->num_knots * Y_axis->num_knots - 1; Y_coord > Y_center_coord;)
	{
		int X_ind = i % X_axis->num_knots;

		real X_coord = X_axis->all_knots[X_ind];
		if (X_coord > X_center_coord)
			fictitious_knots.push_back(i);

		int Y_ind = --i / X_axis->num_knots;
		Y_coord = Y_axis->all_knots[Y_ind];
	}

	return fictitious_knots;
}

void account_fictitious_knots(matrix* A, real* b, vector <int> fictitious_knots)
{
	int all_num_fictitious_knots = fictitious_knots.size();

	for (int i = 0; i < all_num_fictitious_knots; i++)
	{
		int ind = fictitious_knots[i];
		b[ind] = 0;

		A->dl_1[ind - 1] = 0;
		A->dl_offset[ind - A->offset] = 0;
		A->du_1[ind] = 0;
		A->du_offset[ind] = 0;
		A->d_0[ind] = 1;
	}
}



void write_axis(axis* ax)
{
	cout << "axis:" << endl;

	cout << "Главные узлы:" << endl;

	for (int i = 0; i < 3; i++)
		cout << ax->main_knots[i] << endl;

	cout << "/-------------------------" << endl;
	cout << "Все узлы. Количество =" << ax->num_knots << endl;

	for (int i = 0; i < ax->num_knots; i++)
		cout << ax->all_knots[i] << endl;
}


real skal(int i, matrix* mat, real* xk) // Скалярное произведение i-ой строки матрицы с вектором xk
{
	real sum = mat->d_0[i] * xk[i];
	int n = mat->n;
	int m = mat->offset-1;
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

real err(matrix* mat, real* f, real* xk, real f_norm)
{
	int n = mat->n;
	real sum = 0;
	real buf;
	for (int i = 0; i < n; i++)
	{
		buf = f[i] - skal(i, mat, xk);
		sum += buf * buf;
	}

	return sqrt(sum) / f_norm;
}

int Gausse_Zeidel(matrix* mat, real* f, real* xk, real eps, int maxiter)
{
	write_A(mat);
	cout << "Gausse Zeidel" << endl;
	int n = mat->n;
	int k; // Счётчик итераций
	int i;
	real* d0 = mat->d_0;
	real f_norm = 0;
	double W = 1;
	for (i = 0; i < n; i++)
		f_norm += f[i] * f[i];
	f_norm = sqrt(f_norm);



	real err_k = err(mat, f, xk, f_norm);
	for (k = 0; k < maxiter && err_k > eps; k++)
	{
		for (i = 0; i < n; i++)
			xk[i] += W / d0[i] * (f[i] - skal(i, mat, xk));
		/*
		cout << "xk" << endl;
		for (i = 0; i < n; i++)
		cout << xk[i] << endl;
		cout << endl;
		*/
		err_k = err(mat, f, xk, f_norm);
		cout << "Iteration Number " << k << " Error: " << "  " << err_k << endl;
	}
	if (k == maxiter)
		cout << "Limit of iterations is overrided";

	return 0;
}

int main()
{
	setlocale(LC_ALL, "rus");

	ifstream X_input("x_input.txt");
	axis* X_axis = read_axis(X_input);

	write_axis(X_axis);

	ifstream Y_input("y_input.txt");
	axis* Y_axis = read_axis(Y_input);

	write_axis(Y_axis);

	ifstream coeff_input("coeff_input.txt");
	double lambda, gamma;
	coeff_input >> lambda >> gamma;

	matrix* A = build_A(X_axis, Y_axis, lambda, gamma);
	/*write_A(A);*/

	real* b = build_b(X_axis, Y_axis);

	/*cout << "/-------------------" << endl;
	cout << "Вектор b" << endl;

	for (int i = 0; i < A->n; i++)
		cout << b[i] << " ";*/

	vector <int> bound_knots = count_bound_knots(X_axis, Y_axis);
	account_bound(X_axis, Y_axis, A, b, bound_knots);

	/*write_A(A);


	cout << "/-------------------" << endl;
	cout << "Вектор b" << endl;

	for (int i = 0; i < A->n; i++)
		cout << b[i] << " ";*/

	vector <int> fictitious_knots = count_fictitious_knots(X_axis, Y_axis);
	account_fictitious_knots(A, b, fictitious_knots);


	write_A(A);


	cout << "/-------------------" << endl;
	cout << "Вектор b" << endl;

	for (int i = 0; i < A->n; i++)
		cout << b[i] << " ";

	double* xk = new double[A->n]{1};
	Gausse_Zeidel(A, b, xk, 1e-8, 1000000000);

	cout << "/-------------------" << endl;
	cout << "Вектор q" << endl;

	for (int i = 0; i < A->n; i++)
		cout << xk[i] << endl;

	return 0;
}