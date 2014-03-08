#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

static const double c_pi = 3.14159265358979323846;
static const double c_2pi = 2 * c_pi;

double clamp(double v, double low, double high)
{
	return v < low ? low : (v > high ? high : v);
}

int total_points(double eps, int Nphi, int * Nthetas = nullptr)
{
	double s = 0.5 * c_pi / (Nphi - 1);
	int num_points = 0;
	for (int j = 0, k = Nphi - 1; j <= k; ++j, --k)
	{
		int Ntheta = 1;
		if (j > 0)
		{
			double phi_hat = j * c_pi / (Nphi - 1);
			double n = cos(eps) - cos(phi_hat) * cos(phi_hat + s);
			double d = sin(phi_hat) * sin(phi_hat + s);
			double q = n / d;
			q = clamp(q, -1.0, 1.0);
			Ntheta = static_cast<int>(ceil(c_pi / acos(q)));
		}

		if (Nthetas != nullptr) Nthetas[j] = Nthetas[k] = Ntheta;	
		num_points += j < k ? 2 * Ntheta : Ntheta;
	}
	return num_points;
}

// Given eps, find the optimal Nphi
int opt_total_points(double eps, int *optNphi = nullptr)
{
	int Nphi = static_cast<int>(ceil(0.5 * c_pi / eps)) + 1;
	int opt_points = total_points(eps, Nphi);
	int max_search_step = 1000;
	while (max_search_step--)
	{
		int points = total_points(eps, Nphi + 1);
		if (points > opt_points) break;
		Nphi += 1;
		opt_points = points;
	}
	if (optNphi != nullptr) *optNphi = Nphi;
	return opt_points;
}

// find the optimal epsilon by binary search
double search_opt_eps()
{
	const int max_points = (1 << 16);
	const double stop_criterion = 1e-9;
	double low = 0.01;
	double high = 1.0;
	while (opt_total_points(low) <= max_points)
	{
		low *= 0.5;
	}

	while (high - low > stop_criterion)
	{
		double mid = (low + high) / 2;
		if (opt_total_points(mid) > max_points)
		{
			low = mid;
		}
		else
		{
			high = mid;
		}
	}
	return high;
}

struct bin
{
	int16_t i;
	int16_t n;
	bool operator < (bin b) const {return n < b.n || (n == b.n && i < b.i);}
};

struct alias
{
	int16_t n;
	int16_t k;
	int16_t d;
};

struct reference
{
	int16_t i;
	int16_t n;
	int16_t k;
	int16_t m;

	bool operator < (reference r) const
	{
		// can't be both equal, so don't need to compare k or m
		return i < r.i || (i == r.i && n < r.n);
	}
};

void precalc()
{
	double eps = search_opt_eps();
	int Nphi = 0;
	opt_total_points(eps, &Nphi);

	std::vector<int> Nthetas(Nphi);
	int num_points = total_points(eps, Nphi, &Nthetas[0]);
	int Ntheta_max = *std::max_element(Nthetas.begin(), Nthetas.end());
	int Ntheta_avg = (num_points + Nphi - 1) / Nphi;

	std::cout << "total points: " << num_points << std::endl;
	std::cout << "max error (in rad): " << eps << std::endl;
	std::cout << "max error (in degree): " << eps * 180 / c_pi << std::endl;

	std::cout << "Nphi: " << Nphi << std::endl;
	std::cout << "Ntheta max: " << Ntheta_max << std::endl;
	std::cout << "Ntheta avg: " << Ntheta_avg << std::endl;

	std::cout << "alias table size: " << sizeof(alias) * Nphi << " bytes" << std::endl;

	// calculate alias table
	std::set<bin> bins;
	for (int i = 0; i < Nphi; ++i)
	{
		bin b = {i, Nthetas[i]};
		bins.insert(b);
	}

	std::vector<alias> alias_table(Nphi);
	std::vector<reference> reference_table;
	bin last_bin;

	while (!bins.empty())
	{
		if (bins.size() == 1)
		{
			bin b = *bins.begin();
			last_bin = b;
			alias_table[b.i].n = b.n;
			alias_table[b.i].k = -1;
			alias_table[b.i].d = -1;
			bins.clear();

			reference ref = {b.i, b.n, b.i, 0};
			reference_table.push_back(ref);
		}
		else
		{
			bin a = *bins.begin();
			bin b = *bins.rbegin();
			alias_table[a.i].n = a.n;
			alias_table[a.i].k = b.i;
			alias_table[a.i].d = b.n - Ntheta_avg;

			bin c = {b.i, b.n - (Ntheta_avg - a.n)};
			bins.erase(a);
			bins.erase(b);
			bins.insert(c);

			reference ref_a = {a.i, a.n, a.i, 0  };
			reference ref_b = {b.i, b.n, a.i, a.n};
			reference_table.push_back(ref_a);
			reference_table.push_back(ref_b);
		}
	}

	std::sort(reference_table.begin(), reference_table.end());

	std::vector<int16_t> reference_count(Nphi);
	std::vector<int16_t> offset_table(Nphi);
	for (int i = 0; i != reference_table.size(); ++i)
	{
		reference_count[reference_table[i].i] += 1;
		if (i > 0 && reference_table[i].i > reference_table[i-1].i)
			offset_table[reference_table[i].i] = i;
	}

	std::cout << "reference table size: " << sizeof(int16_t) * 3 * reference_table.size() << " bytes" << std::endl;
	std::cout << "offset table size: " << sizeof(int16_t) * Nphi << " bytes" << std::endl;
	std::cout << "max reference:" << *max_element(reference_count.begin(), reference_count.end()) << std::endl;
	std::cout << "avg reference:" << std::accumulate(reference_count.begin(), reference_count.end(), 0) / static_cast<float>(Nphi) << std::endl;

	std::ofstream ofs("out.txt");

	// write Nphi
	ofs << "static const int16_t Nphi = " << Nphi << ";" << std::endl << std::endl;

	// write Nthetas
	ofs << "static const int16_t Nthetas[] = {";
	for (int i = 0; i < Nphi; ++i)
	{
		if (i % 16 == 0) ofs << std::endl << "    ";
		ofs << Nthetas[i] << ", ";
	}
	ofs << std::endl << "};" << std::endl << std::endl;

	// write alias table
	ofs << "static const alias alias_table[] = {";
	for (int i = 0; i < Nphi; ++i)
	{
		if (i % 4 == 0) ofs << std::endl << "    ";
		ofs << "{" << alias_table[i].n << ", " << alias_table[i].k << ", " << alias_table[i].d << "}, ";
	}
	ofs << std::endl << "};" << std::endl << std::endl;

	// write gap info
	ofs << "static const int gap_begin  = " << last_bin.i * Ntheta_avg + last_bin.n << ";" << std::endl;
	ofs << "static const int gap_end    = " << (last_bin.i + 1) * Ntheta_avg << ";" << std::endl;
	ofs << "static const int gap_length = gap_end - gap_begin;" << std::endl;
	ofs << std::endl;

	// write out reference table
	ofs << "static const reference reference_table[] = {";
	for (int i = 0; i != reference_table.size(); ++i)
	{
		if (i % 5 == 0) ofs << std::endl << "    ";
		ofs << "{" << reference_table[i].n << ", " << reference_table[i].k << ", " << reference_table[i].m <<  "}, ";
	}
	ofs << std::endl << "};" << std::endl << std::endl;

	// write out offset table
	ofs << "static const int16_t offset_table[] = {";
	for (int i = 0; i < Nphi; ++i)
	{
		if (i % 16 == 0) ofs << std::endl << "    ";
		ofs << offset_table[i] << ", ";
	}
	ofs << std::endl << "};" << std::endl << std::endl;
}
