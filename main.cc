#include <fstream>
#include <vector>
#include <cstdint>
#include <random>
#include <iostream>
#include <algorithm>
#include <utility>
#include <string>

void read_file(const char*, std::vector<std::vector<double>>&);
void initialize_p(std::vector<double>&, double, unsigned);
std::vector<std::pair<unsigned, double>> simulate_network(
  std::vector<double>&, 
  std::vector<std::vector<double>>&, 
  double, 
  double,
  size_t
);
void write_output(double, std::vector<std::pair<unsigned, double>>&);

int main(int argc, char const *argv[])
{
  const char* rij_data = argv[1];   // The filename of the link probability matrix
  size_t m = atoi(argv[2]);         // The matrix dimension
  unsigned lambda = atoi(argv[3]);  // The number of trials
  double beta = atof(argv[4]);      // The transmission rate
  double mu = atof(argv[5]);        // The recovery rate
  double rho_0 = atof(argv[6]);     // The initial infection probability
  size_t max_step = atoi(argv[7]);  // The number of simulation steps
  unsigned seed = int(argv[8]);     // The seed
  std::vector<std::vector<double>> r(m);  // The link probability matrix
  for (auto& row : r)
    row.resize(m);
  std::vector<double> p(m);

  read_file(rij_data, r); // read the r matrix
  initialize_p(p, rho_0, seed); // initialize p
  auto rho_ts = std::move(simulate_network(p, r, beta, mu, max_step));
  write_output(beta, rho_ts);
  return 0;
}

void read_file(const char* filename, std::vector<std::vector<double>>& r) {
  std::ifstream ifs(filename);
  double x;
  for (auto& row : r)
    for (auto& element: row)
      ifs >> element;
  ifs.close();
}

void initiliaze_p(std::vector<double>& p, double rho_0, unsigned seed) {
  std::mt19937 rng(seed);
  std::bernoulli_distribution bernoulli(rho_0);
  for (auto& p_i : p)
    p_i = bernoulli(rng);
}

std::vector<std::pair<unsigned, double>> simulate_network(
  std::vector<double>& p, 
  std::vector<std::vector<double>>& r,
  double beta,
  double mu,
  size_t max_step
) {
  double rho;
  std::vector<std::pair<unsigned, double>> rho_ts;
  std::vector<double> q(p.size()*p.size());
  for (size_t step = 0; step < max_step; step++) {
    for (size_t i = 0; i < p.size(); i++) {
      q[i] = 1;
      for (size_t j = 0; j < p.size(); j++) 
        q[i] = (1 - beta*r[j][i]*p[j]);
      p[i] = (1 - q[i]) * (1 - p[i]) + (1 - mu)*(p[i]) + mu*(1 - q[i])*p[i];
    }
    rho = std::accumulate(p.begin(), p.end(), 0.0) / p.size(); 
    rho_ts.push_back(std::make_pair(step, rho));
  }
  return rho_ts;
}

void write_output(
  double beta,
  std::vector<std::pair<unsigned, double>>& rho_ts
) {
  std::string filename = "beta_" + std::to_string(beta) + ".dat";
  std::ofstream ofs(filename);
  for (size_t i = 0; i < rho_ts.size(); i++) {
    ofs << rho_ts[i].first << ' ' << rho_ts[i].second;
    if (i < rho_ts.size())
      ofs << '\n';
  }
}

// double InfectionObserver::compute_rho() {
//   //(1−q_i(t))(1−p_i(t)) + (1−μ)p_i(t) + μ(1−q_i(t))p_i(t)
//   //q_i (t)= ∏24_(j=1)^N▒(1−βr_ij p_j (t)) 

//   for (size_t i = 0; i < host_num; i++) {
//     (*q)[i] = 1.0; //not infected by a neighbor
//     for (auto&& entry : *(adjacency_matrix->at(i)))
//       (*q)[i] *= (1 - beta * (*p)[entry.node_index]);
//     (*next_p)[i]                         //p(t+1)
//       = (1 - (*q)[i]) * (1 - (*p)[i])    //infected by nbr & not infected
//       + (1 - mu) * (*p)[i]               //not recovery & infected
//       + mu * (1 - (*q)[i]) * (*p)[i];    //infected after recovery
//     EV_INFO << "next_p[" << i << "] = " << (*next_p)[i] << '\n';
//   }
//   return std::accumulate(next_p->begin(), next_p->end(), 0.0) / host_num; //rho
// }