#include <algorithm>
#include <map>
#include <numbers>
#include <format>
#include <cmath>
#include <limits>
#include <ratio>
#include <span>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <utility>
#include <complex>
#define lapack_complex_float  std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>

using cd = std::complex<double>;
using vecd = std::vector<double>;

double hbar = 1,
	   m = 1,
	   a = 1;

struct CNSystem
{
	vecd x, V;
	double dx, dt;
	cd r, beta;
	std::vector<cd> Amain, Asub, Bmain;
};

struct EvolutionResult
{
	vecd time, PL, norm;
};

struct EigenResult
{
	vecd E;
	std::vector<vecd> psi;
};

struct ComparisonResult
{
	double dE, dECN, dEm, methodError, measureError;
};

double measureTime(auto&& f)
{
	const auto t0 = std::chrono::steady_clock::now();
	f();
	const auto t1 = std::chrono::steady_clock::now();
	return std::chrono::duration<double, std::milli>(t1-t0).count();
}

CNSystem buildSystem(const double V0, const double L, const int Nx, const double Tfin,
		const int Nt)
{
	using namespace std::complex_literals;
	CNSystem s;

	s.dx = 2.0 * L / (Nx - 1);
	s.dt = Tfin / Nt;

	s.x.resize(Nx);
	s.V.resize(Nx);
	for(int i = 0; i < Nx; i++)
	{
		s.x[i] = -L + i * s.dx;
		double t = 1.0 - s.x[i] * s.x[i] / (a * a);
		s.V[i] = V0 * t * t;
	}

	s.r = 1i * hbar * s.dt / (4 * m * s.dx * s.dx);
	s.beta = 1i * s.dt / (2 * hbar);
	
	const int n = Nx - 2;
	s.Amain.resize(n);
	s.Asub.assign(n-1, -s.r);
	s.Bmain.resize(n);

	for(int i = 0; i < n; i++)
	{
		s.Amain[i] = 1.0 + 2.0 * s.r + s.beta * s.V[i+1];
		s.Bmain[i] = 1.0 - 2.0 * s.r - s.beta * s.V[i+1];
	}

	return s;
}

EigenResult eigenPsi(const vecd& V, const vecd& x, int k = 2)
{
	const lapack_int N = static_cast<lapack_int>(x.size()),
		  n = N - 2;
	const double dx = x[1] - x[0],
		  t = hbar * hbar / (2.0 * m * dx * dx);

	vecd Hmain(n), Hsub(n - 1, -t);
	for (lapack_int i = 0; i < Hmain.size(); i++)
		Hmain[i] = 2.0 * t + V[i+1];
	
	lapack_int found = 0;
	vecd w(n), z(static_cast<size_t>(n) * k);
	std::vector<lapack_int> ifail(n);
	const double abstol = 2.0 * LAPACKE_dlamch('S');

	lapack_int info = LAPACKE_dstevx(
			LAPACK_COL_MAJOR,
			'V',
			'I',
			n, Hmain.data(), Hsub.data(),
			0.0, 0.0,
			1, k,
			abstol, &found,
			w.data(), z.data(), n,
			ifail.data()
			);

	if(info != 0)
		throw std::runtime_error("LAPACKE_dstevx failed, info = " +
				std::to_string(info));

	EigenResult res;
	res.E.assign(w.begin(), w.begin() + found);
	res.psi.assign(found, vecd(N, 0.0));

	const double s = 1.0 / std::sqrt(dx);
	for(lapack_int i = 0; i < found; i++)
		for(lapack_int j = 0; j < n; j++)
			res.psi[i][j+1] = s * z[i * n + j];

	return res;
}

std::vector<cd> localizedState(std::vector<vecd>& eig, int sign = +1)
{
	const auto& p0 = eig[0], p1 = eig[1];
	std::vector<cd> psi(p0.size());
	for(size_t i = 0; i < p0.size(); i++)
		psi[i] = (p0[i] + sign * p1[i]) / std::sqrt(2.0);
	return psi;
}

double totalNorm(const std::span<cd>& psi, const double dx)
{
	double s = 0.0;
	for(const auto& z : psi) s += std::norm(z);
	return s * dx;
}

EvolutionResult evolve(const CNSystem& s, std::vector<cd>& psi, const int Nt, 
		const int n_out)
{
	int save_every = Nt / n_out;
	const int Nx = static_cast<int>(psi.size());
	const lapack_int n = Nx - 2;

	std::vector<cd> dl = s.Asub, d = s.Amain, du = s.Asub, du2(n-2);
	std::vector<lapack_int> ipiv(n);
	lapack_int info = LAPACKE_zgttrf(n, dl.data(), d.data(), du.data(), du2.data(),
			ipiv.data());

	if(info != 0) throw std::runtime_error("zgttrf: info = " + std::to_string(info));

	const size_t k = std::lower_bound(s.x.begin(), s.x.end(), -0.5 * s.dx) - s.x.begin();
	const bool centerOnGrid = std::abs(s.x[k]) < 0.5 * s.dx;

	EvolutionResult res;
	std::vector<cd> rhs(Nx, 0.0);
	for(int i = 0; i < Nt; i++)
	{
		if(i % save_every == 0)
		{
			res.time.push_back(i * s.dt);
			res.PL.push_back(totalNorm(std::span(psi).first(k), s.dx));
			if(centerOnGrid) res.PL.back() += 0.5 * std::norm(psi[k]) * s.dx;
			res.norm.push_back(totalNorm(psi, s.dx));
		}

		for(int j = 1; j < Nx - 1; j++)
			rhs[j] = s.r * psi[j-1] + s.Bmain[j-1] * psi[j] + s.r * psi[j+1];

		info = LAPACKE_zgttrs(LAPACK_COL_MAJOR, 'N', n, 1, 
				dl.data(), d.data(), du.data(), du2.data(), 
				ipiv.data(), rhs.data() + 1, n);

		if(info != 0) throw std::runtime_error("zgttrs: info = " + std::to_string(info));

		std::swap(psi, rhs);
	}

	return res;
}

double measurePeriod(const vecd& time, const vecd& PL)
{    
	if (PL.size() < 2) throw std::runtime_error("measurePeriod: too few points");

	vecd tc;
	for(size_t i = 0; i < PL.size() - 1; i++)
	{
		const double f0 = PL[i] - 0.5, f1 = PL[i + 1] - 0.5;
		if(f0 * f1 < 0.0)
			tc.push_back(time[i] - f0 * (time[i+1] - time[i]) / (f1 - f0));
	}

	if(tc.size() < 2)
		return std::numeric_limits<double>::quiet_NaN();

	return 2.0 * (tc.back() - tc.front()) / (tc.size() - 1);
}

ComparisonResult compareSplitting(const vecd& E, const double dt, const double T)
{
	ComparisonResult res;
	res.dE = E[1] - E[0];
	res.dECN = hbar * 2 * (std::atan(E[1] * dt / (2 * hbar)) - 
			std::atan(E[0] * dt / (2 * hbar))) / dt;
	res.dEm = 2.0 * std::numbers::pi * hbar / T;
	res.methodError = (res.dECN - res.dE) / res.dE;
	res.measureError = (res.dEm - res.dECN) / res.dECN;
	return res;
}

std::map<int, double> evolutionTime;
int iterations = 50;

void Crank_Nicolson(double V0, double L, int Nx, double Tfin, int Nt, int n_out)
{
	auto sys = buildSystem(V0, L, Nx, Tfin, Nt);

	auto [E, psi] = eigenPsi(sys.V, sys.x, 2);
	auto psi0 = localizedState(psi);

	EvolutionResult evo;
	double tEvolve = measureTime([&] { evo = evolve(sys, psi0, Nt, n_out); });
	auto& [time, PL, norm] = evo;
	evolutionTime[Nt] += tEvolve / iterations;
	//std::cout << "Evolution time with Nt = " << Nt << ":\t" << tEvolve << '\n';

	double period = measurePeriod(time, PL);

	ComparisonResult compare = compareSplitting(E, sys.dt, period);

//	std::cout << std::setprecision(15);
//
//	std::cout << Nt << "\n\n";
//	std::cout << std::format(
//  	  "dE (spectrum)  = {:.8e}\n"
//  	  "dE (CN)        = {:.8e}\n"
//  	  "dE_m           = {:.8e}\n"
//  	  "method error   = {:.3e}\n"
//  	  "measure error  = {:.3e}\n",
//  	  compare.dE, compare.dECN, compare.dEm, compare.methodError, compare.measureError);
//	std::cout << "\n\n";
}

int main()
{
	std::vector<int> times {2000, 5000, 10000, 20000};

	for(auto i : times) evolutionTime[i] = 0;

	for(int j = 0; j < iterations; j++)
	{
		for(auto i : times) 
			Crank_Nicolson(8.0, 3.0, 1001, 100.0, i, 2000);
	}

	for(auto i : evolutionTime)
		std::cout << i.first << ":\t" << i.second << '\n';

	return 0;
}
