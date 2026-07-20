#include "dlp.h"
#include <future>
#include <mutex>


namespace dlp {
  static std::mutex crt_solve_mtx;
  
  DLP::DLP(u128 p) : params(p) {}

  u128 DLP::solve(u128 g, u128 b) {

    size_t k = params.factor_base.size();

    RelationMatrix Mg, Mb;
    U128Vector Xg, Xb;

    for (size_t blks = 1; blks <= 16; blks *= 2) {
      
      auto fut_g = std::async(std::launch::async, [&]() -> std::vector<std::pair<uint32_t, u128>>
      {
        infra::addRelations(g, blks * k, params, Mg, Xg);
	
	U128Vector Lg;
	{
	  std::lock_guard<std::mutex> lock(crt_solve_mtx);
	  Lg = infra::crtSolve(params, Mg, Xg);
	}
	
	if (Lg.empty()) return {};
	return infra::filterDetermined(g, params, Lg);
      });
      
      auto fut_b = std::async(std::launch::async, [&]() -> std::vector<std::pair<uint32_t, u128>>
      {
	infra::addRelations(b, blks * k, params, Mb, Xb);
	U128Vector Lb;
	{
	  std::lock_guard<std::mutex> lock(crt_solve_mtx);
	  Lb = infra::crtSolve(params, Mb, Xb);
	}
	
	if (Lb.empty()) return {};
	return infra::filterDetermined(b, params, Lb);
      });

      auto og = fut_g.get();
      auto ob = fut_b.get();

      if (og.empty() || ob.empty()) continue;

      size_t ib = 0;

      for (size_t ig = 0; ig < og.size(); ig++) {
	while (ib < ob.size() && ob[ib].first < og[ig].first) ib++;
	if (ib >= ob.size()) break;

	auto [p1, alpha] = og[ig];
	auto [p2, beta] = ob[ib];

	if (p1 == p2) {
	  u128 beta_inv = infra::modInv(beta, params.p_factorization);
	  if (beta_inv == 0) continue;

	  return mont::mulmod_any(alpha, beta_inv, params.p - 1);
	}
      }
    }

    return 0;
  }
}
