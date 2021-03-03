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
	vector<double> d_0;   // главная диагональ
	vector<double> du_1;   // верхняя диагональ со смещением 1
	vector<double> du_offset;   // верхняя диагональ со смещением offset
	vector<double> dl_1;	// нижняя диагональ со смещением 1
	vector<double> dl_offset;	// нижняя диагональ со смещением offset

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

enum sides
{
	s_left,
	s_right,
	s_bottom,
	s_top
};

double CountFirstStep(double begin, double end, int num_steps, double q);
void ReadGrid(string name_file, vector<double>* grid, vector<int>* nums_main_knots);
void ReadRects(string name_file, int num_rects_x, int num_rects_y, vector<vector<bool>>* fullness);
void FillingKnotsIN(double lambda, double gamma, vector<double> x, vector<double> y, matrix* A, double* b);
void SearchBoundKnots(int x_num_all_knots, int y_num_all_knots, vector<int>* bound_knots);
void SearchFictKnots(vector<int> x_num_mains, vector<int> y_num_mains, vector<vector<bool>> fullness, vector<int>* fict_knots);
void checkFictKnots(vector<int> x_num_mains, vector<int> y_num_mains, vector<int>* fict_knots, vector<int>* bound_knots);
void AccountBound(vector<double> x, vector<double> y, matrix* A, double* b, vector <int> bound_knots);
void AccountFict(int x_size, matrix* A, double* b, vector <int> fict_knots);
double Skal(int i, matrix* mat, double* xk);
double Err(matrix* mat, double* f, double* xk, double f_norm);
int Gausse_Zeidel(matrix* mat, double* f, double* xk, double eps, int maxiter);

double count_bound_func(double X, double Y)
{
	return X * X * X;
}

double count_func(double X, double Y)
{
	return 2.*X*X*X - 6. * X;
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
	vector<int> x_nums_main_knots;
	ReadGrid("x_input.txt", &x_grid, &x_nums_main_knots);

	vector<double> y_grid;
	vector<int> y_nums_main_knots;
	ReadGrid("y_input.txt", &y_grid, &y_nums_main_knots);

	vector<vector<bool>> fullness;
	ReadRects("rects.txt", x_nums_main_knots.size() - 1, y_nums_main_knots.size() - 1, &fullness);

	ifstream coeff_input("coeff_input.txt");
	double lambda, gamma;
	coeff_input >> lambda >> gamma;
	coeff_input.close();

	matrix* A = new matrix(x_grid.size() * y_grid.size(), x_grid.size());
	double* b = new double[A->n];
	FillingKnotsIN(lambda, gamma, x_grid, y_grid, A, b);
	vector <int> bound_knots;
	SearchBoundKnots(x_grid.size(), y_grid.size(), &bound_knots);
	
	vector<int> fict_knots;
	SearchFictKnots(x_nums_main_knots, y_nums_main_knots, fullness, &fict_knots);

	checkFictKnots(x_nums_main_knots, y_nums_main_knots, &fict_knots, &bound_knots);

	AccountBound(x_grid, y_grid, A, b, bound_knots);
	AccountFict(x_grid.size(),A,b,fict_knots);
	
	double* q = new double[A->n];
	Gausse_Zeidel(A, b, q, 1e-11, 1e+11);

	cout << "Результат : \n";
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

void ReadGrid(string name_file, vector<double>* grid, vector<int>* nums_main_knots )
{
	ifstream in(name_file);

	int num_main_knots;

	in >> num_main_knots;
	int all_num_steps = 0;
	double coord;
	for (int i = 0; i < num_main_knots - 1; i++)
	{
		double begin, end, q;
		int  num_steps;
		in >> begin >> end;
		
		in >> num_steps;
		nums_main_knots->push_back(all_num_steps);
		all_num_steps += num_steps;
		in >> q;

		double h0 = CountFirstStep(begin, end, num_steps, q);

		double h = h0;

		coord = i ? begin + h : begin;
	
		while (coord < end)
		{
			grid->push_back(coord);
			coord += h;
			h *= q;
		}

		grid->push_back(end);
	}
	nums_main_knots->push_back(all_num_steps);
	in.close();
}

void ReadRects(string name_file, int num_rects_x, int num_rects_y, vector<vector<bool>>* fullness)
{
	ifstream in(name_file);

	bool is_full;
	for (int i = 0; i < num_rects_y; i++)
	{
		vector<bool> str;
		for (int j = 0; j < num_rects_x; j++)
		{
			in >> is_full;
			str.push_back(is_full);
		}
		fullness->push_back(str);
	}
}

void FillingKnotsIN(double lambda, double gamma, vector<double> x, vector<double> y, matrix* A, double* b)
{

	int x_max_ind = x.size() - 1;
	int y_max_ind = y.size() - 1;
	
	// Проходим по узлам, не затрагивая нижнюю и верхнюю грани
	for (int i = 1; i < x_max_ind; i++)
	{
		for (int j = 1 ; j < y_max_ind; j++)
		{
			double h_right = x[i + 1] - x[i];	// шаг справа от узла
			double h_left = x[i] - x[i - 1];	// шаг слева от узла

			double h_top = y[j + 1] - y[j];		// шаг сверху от узла
			double h_bottom = y[j] - y[j - 1];	// шаг снизу от узла

			int ind = (x_max_ind + 1) * j + i;
			
			A->d_0[ind] = lambda * 2. * ((h_bottom * h_top + h_left * h_right) / (h_bottom * h_top * h_left * h_right)) + gamma;
			A->dl_1[ind - 1] = -lambda * 2. / (h_left * (h_right + h_left));
			A->du_1[ind] = -lambda * 2. / (h_right * (h_right + h_left));
			A->dl_offset[ind - A->offset] = -lambda * 2. / (h_bottom * (h_top + h_bottom));
			A->du_offset[ind] = -lambda * 2. / (h_top * (h_top + h_bottom));
			b[ind] = count_func(x[i], y[j]);
		}
	}
}

void SearchBoundKnots(int x_num_all_knots, int y_num_all_knots, vector<int>* bound_knots)
{
	int left_border = 0;
	int right_border = x_num_all_knots - 1;

	int bottom_border = 0;
	int top_border = y_num_all_knots - 1;

	for (int i = 0; i < x_num_all_knots; i++)
	{
		for (int j = 0; j < y_num_all_knots; j++)
		{
			int ind = x_num_all_knots * j + i;
			if (j == bottom_border || j == top_border)
			{
				bound_knots->push_back(ind);
			}
			else
			{
				if (i == left_border || i == right_border)
				{
					bound_knots->push_back(ind);
				}
			}
		}
	}
}

bool findKnot(vector<int> knots, int knot)
{
	vector<int>::iterator it = knots.begin();

	while (it != knots.end()) {
		if (*it == knot) {
			return 1;
			break;
		}
		it++;
	}
	if (it == knots.end()) return 0;
}

void SearchFictKnots(vector<int> x_num_mains, vector<int> y_num_mains, vector<vector<bool>> fullness, vector<int>* fict_knots)
{
	int x_size = x_num_mains[x_num_mains.size() - 1] + 1;

	int i_num_rects_y = y_num_mains.size() - 1;
	int j_num_rects_x = x_num_mains.size() - 1;

	for (int i = 0; i < i_num_rects_y; i++)
	{
		int y_begin = y_num_mains[i];
		int y_end = y_num_mains[i + 1];

		for (int j = 0; j < j_num_rects_x; j++)
		{
			int x_begin = x_num_mains[j];
			int x_end = x_num_mains[j + 1];

			bool full = fullness[i][j];

			if (!full)	// Если прямоугольник фиктивный
			{
				for (int y = y_begin; y <= y_end; y++)
				{
					for (int x = x_begin; x <= x_end; x++)
					{
						int ind = x_size * y + x;
						if (!findKnot(*fict_knots, ind))
							fict_knots->push_back(ind);	// Записываем все узлы прямоугольника как фиктивные
					}
				}
			}
		}
	}
}

void findIsEmptyNearKnots(vector<int> knots, int ind, int x_size, bool* empty)
{
	empty[0] = findKnot(knots, ind - 1);	// Пустой ли левый узел
	empty[1] = findKnot(knots, ind + 1);	// Пустой ли правый узел
	empty[2] = findKnot(knots, ind - x_size);	// Пустой ли нижний узел
	empty[3] = findKnot(knots, ind + x_size);	// Пустой ли верхний узел
}

void findNearKnots(int ind, int x_size, int* nears)
{
	int x_ind = ind % x_size;
	int y_ind = ind / x_size;

	nears[s_left] = x_ind - 1;
	nears[s_right] = x_ind + 1;
	nears[s_bottom] = y_ind + 1;
	nears[s_top] = y_ind - 1;
}

void checkCornerKnots(int* borders, vector<int> old_knots, vector<int> fict_knots, vector<int>* new_knots)
{
	int x_size = borders[s_right] + 1;
	// перепроверка. вдруг узловые
	for (int i = 0; i < old_knots.size(); i++)
	{
		int ind = old_knots[i];

		int nears[4];
		findNearKnots(ind, x_size, nears);

		if (nears[s_left] >= borders[s_left] &&
			nears[s_right] <= borders[s_right] &&
			nears[s_bottom] >= borders[s_bottom] &&
			nears[s_top] <= borders[s_top])
		{
			// Рассматриваем узлы не на границе
			bool empty[4];
			findIsEmptyNearKnots(fict_knots, ind, x_size, empty);

			if (empty[s_left] || empty[s_right] || empty[s_bottom] || empty[s_top])	// Если хотя бы один узел пустой
			{
				new_knots->push_back(ind);
			}
		}
		else
			new_knots->push_back(ind);
	}
}

void checkAllKnots(int* borders, vector<int> fict_knots, vector<int>* new_fict_knots, vector<int>* new_bound_knots)
{
	int x_size = borders[s_right] + 1;
	// перепроверка. вдруг узловые
	for (int i = 0; i < fict_knots.size(); i++)
	{
		int ind = fict_knots[i];

		int nears[4];
		findNearKnots(ind, x_size, nears);

		bool empty[4];
		findIsEmptyNearKnots(fict_knots, ind, x_size, empty);

		// Проверка если узел не в сетке
		// ----
		if (nears[s_left] < borders[s_left]) empty[s_left] = 1;

		if (nears[s_right] > borders[s_right])	empty[s_right] = 1;

		if (nears[s_bottom] < borders[s_bottom]) empty[s_bottom] = 1;

		if (nears[s_top] > borders[s_top])	empty[s_top] = 1;
		// ----

		if (empty[s_left] && empty[s_right] && empty[s_bottom] && empty[s_top])	// Если все узлы пустые
		{
			new_fict_knots->push_back(ind);		// Записываем узел как фиктивный
		}
		else
		{
			new_bound_knots->push_back(ind);	// Записываем узел как граничный
		}
	}
}

void checkFictKnots(vector<int> x_num_mains, vector<int> y_num_mains, vector<int>* fict_knots, vector<int>* bound_knots)
{
	int borders[4];
	borders[s_left] = 0;
	borders[s_right] = x_num_mains[x_num_mains.size() - 1];
	borders[s_bottom] = 0;
	borders[s_top] = y_num_mains[y_num_mains.size() - 1];

	vector<int> new_fict_knots, new_bound_knots;
	checkAllKnots(borders, *fict_knots, &new_fict_knots, &new_bound_knots);

	// Перепроверка. вдруг узловые
	vector<int> new_new_fict_knots;
	checkCornerKnots(borders, new_fict_knots, new_fict_knots, &new_new_fict_knots);
	*fict_knots = new_new_fict_knots;

	// перепроверка. вдруг узловые
	checkCornerKnots(borders, new_bound_knots, new_new_fict_knots, bound_knots);

}

void AccountBound(vector<double> x, vector<double> y, matrix* A, double* b, vector <int> bound_knots)
{
	for (int i = 0; i < bound_knots.size(); i++)
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

void AccountFict(int x_size, matrix* A, double* b, vector <int> fict_knots)
{
	for (int i = 0; i < fict_knots.size(); i++)
	{
		int ind = fict_knots[i];

		int X_ind = ind % x_size;
		int Y_ind = ind / x_size;

		b[ind] = 0.;

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

double Skal(int i, matrix* mat, double* xk) // Скалярное произведение i-ой строки матрицы с вектором xk
{
	double sum = mat->d_0[i] * xk[i];
	int n = mat->n;
	int m = mat->offset - 1;

	if (i > 0) sum += mat->dl_1[i - 1] * xk[i - 1];
	if (i < n - 1) sum += mat->du_1[i] * xk[i + 1];

	int ind;
	ind = i - m - 1;
	if (ind > -1) sum += mat->dl_offset[ind] * xk[ind];

	ind = i + m + 1;
	if (ind < n) sum += mat->du_offset[i] * xk[ind];

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
	int k; // Счётчик итераций
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