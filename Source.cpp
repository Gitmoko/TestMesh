
#include <iostream>
#include <cstdlib>
#include <map>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <random>
//-0.5,-0.5,0.5の位置をNo.1ノードとする


struct Pos {
	double x = 0, y = 0, z = 0;
	std::string show() {
		std::string ret;
		ret += std::to_string(x);
		ret += ",";
		ret += std::to_string(y);
		ret += ",";
		ret += std::to_string(z);
		return ret;
	}

	Pos operator+(const Pos& pos) {
		auto ret = Pos{ x + pos.x,y + pos.y,z + pos.z };
		return ret;
	}

	Pos operator*(float c) {
		auto ret = Pos{ x*c,y*c,z*c };
		return ret;
	}
};

#define P (1)
#define Z (0)
#define M (-1)

std::array<Pos, 28> init() {
	std::array<Pos, 28> RelPos;
	RelPos[0] = { Z,Z,Z };
	RelPos[1] = { M,M,M };
	RelPos[2] = { P,M,M };
	RelPos[3] = { P,P,M };
	RelPos[4] = { M,P,M };
	RelPos[5] = { M,M,P };
	RelPos[6] = { P,M,P };
	RelPos[7] = { P,P,P };
	RelPos[8] = { M,P,P };
	RelPos[9] = { Z,M,M };
	RelPos[10] = { P,Z,M };
	RelPos[11] = { Z,P,M };
	RelPos[12] = { M,Z,M };
	RelPos[13] = { Z,M,P };
	RelPos[14] = { P,Z,P };
	RelPos[15] = { Z,P,P };
	RelPos[16] = { M,Z,P };
	RelPos[17] = { M,M,Z };
	RelPos[18] = { P,M,Z };
	RelPos[19] = { P,P,Z };
	RelPos[20] = { M,P,Z };
	RelPos[21] = { Z,Z,M };
	RelPos[22] = { Z,Z,P };
	RelPos[23] = { Z,M,Z };
	RelPos[24] = { P,Z,Z };
	RelPos[25] = { Z,P,Z };
	RelPos[26] = { M,Z,Z };
	RelPos[27] = { Z,Z,Z };
	return RelPos;

}
void move_offset(std::array<Pos, 28>& RelPos,Pos pos) {
	for (auto& e : RelPos) {
		e = e + pos;
	}
}

struct Cell {
	double rho = 0;//密度
	std::tuple<double, double, double> velocity;
	std::tuple<double, double, double> grad_rho = { 0.0,0.0,0.0 };
	double div_velocity = 0;
	std::array<int, 27> twn;
};
struct Generate {

	std::map<int, Pos> Nodes;
	std::map <std::tuple<int,int,int>,Cell> cells;
	double celllength;
	Pos offset = Pos{0,0,0};
	int rho_randseed;
	int velo_randseed;
	int xnum = 3, ynum = 3, znum = 3;
	

	int getindex(int x, int y, int z) {
		auto ret = z + y*(znum*2+1) + x*(znum*2+1)*(ynum*2+1);
		return ret;
	}

	//(x,y,z)の位置にある頂点から双対メッシュをつくる
	std::array<std::tuple<int, int, int>, 27> gettwn(int x, int y, int z) {
		std::array<std::tuple<int, int, int>, 27> ret;
		auto RelPos = init();
	
		for (int i = 1, e = 27; i <= e;i++) {
			std::tuple<int, int, int> n;
			std::get<0>(n) = x + RelPos[i].x;
			std::get<1>(n) = y + RelPos[i].y;
			std::get<2>(n) = z + RelPos[i].z;

			ret[i-1] = n;
		}
		return ret;
	}

	Pos getpos(int x, int y, int z) {
		Pos ret = offset + Pos{ celllength*x/2.0 , celllength*y/2.0 , celllength*z/2.0 };
		return ret;
	}

	Generate(double celllength_) :celllength(celllength_) {
		std::random_device randdev;
		rho_randseed = randdev();
		velo_randseed = randdev();
		std::mt19937 rhoengine(rho_randseed);
		std::mt19937 veloengine(velo_randseed);
		std::uniform_real_distribution<> dist(-5.0, 5.0);

		auto RelPos = init();
		move_offset(RelPos, Pos{ 1,1,1 });

		int index = 0;
		for (auto x = 0; x < xnum; x++) {
			for (auto y = 0; y < ynum; y++) {
				for (auto z = 0; z < znum; z++) {
					
					std::array<int, 27> twn;
					for (auto i = 1, e = 27; i <= e;i++) {
						auto& relpos = RelPos[i];
						auto node_relpos = Pos{x*2.0,y*2.0,z*2.0}+relpos;
						auto node_realpos = getpos(node_relpos.x, node_relpos.y, node_relpos.z);
						auto index = getindex(node_relpos.x,node_relpos.y,node_relpos.z);
						Nodes[index] = node_realpos;
						twn[i - 1] = index;
					}

					cells[{x, y, z}].twn = twn;
					cells[{x, y, z}].rho = dist(rhoengine);
					cells[{x, y, z}].velocity = { dist(veloengine),dist(veloengine),dist(veloengine) };

				}
			}
		}

		calc_grad();
		calc_div();
	}

	template<class Func>
	void cell_foreach(Func f) {
		for (int x = 0; x < xnum; x++) {
			for (int y = 0; y < ynum; y++) {
				for (int z = 0; z < znum; z++) {
					f(cells[{x, y, z}]);
				}
			}
		}
	}

	template<class Func>
	void cell_foreach_withend(Func f) {
		for (int x = 0; x < xnum; x++) {
			for (int y = 0; y < ynum; y++) {
				for (int z = 0; z < znum; z++) {
					f(cells[{x, y, z}], x == xnum - 1 && y == ynum - 1 && z == znum - 1);
				}
			}
		}
	}

	void calc_grad() {
		std::map<std::tuple<int, int, int>, double> ret;
		for (int x = 0; x < xnum; x++) {
			for (int y = 0; y < ynum; y++) {
				for (int z = 0; z < znum; z++) {
					cells[{x, y, z}].grad_rho = { (cells[{x, y, z}].rho - cells[{x - 1, y, z}].rho)/celllength,
						(cells[{x, y, z}].rho - cells[{x, y - 1, z}].rho)/celllength,
							(cells[{x, y, z}].rho - cells[{x, y, z - 1}].rho)/celllength };
				}
			}
		}
	}

	void calc_div() {
		std::map<std::tuple<int, int, int>, double> ret;
		for (int x = 0; x < xnum; x++) {
			for (int y = 0; y < ynum; y++) {
				for (int z = 0; z < znum; z++) {

					cells[{x, y, z}].div_velocity = (std::get<0>(cells[{x, y, z}].velocity) - std::get<0>(cells[{x - 1, y, z}].velocity))
						+ (std::get<1>(cells[{x, y, z}].velocity) - std::get<1>(cells[{x, y - 1, z}].velocity)
							+ std::get<2>(cells[{x, y, z}].velocity) - std::get<2>(cells[{x, y, z - 1}].velocity)) / celllength;
				}
			}
		}
	}

};





int main()
{
	Generate gen{2};


	//node data
	std::ofstream nodefs{"nodes.csv"};
	for (auto i = gen.Nodes.begin(), e = gen.Nodes.end(); i != e; i++) {
		nodefs << i->second.show();
		if (std::next(i) != e) {
			nodefs << std::endl;
		}
	}


	//twenty vertices data
	std::ofstream twnfs{ "twn.csv" };
	auto twnoutput = [&](Cell& cell,bool last) {
		for (auto twnb = 0, size = 20; twnb < size; twnb++) {
			twnfs << cell.twn[twnb] + 1;
			if (twnb + 1 != size) {
				twnfs << ",";
			}
		}

		if (!last) {
			twnfs << std::endl;
		}
	};
	gen.cell_foreach_withend(twnoutput);

	//rho data
	std::string rhoname = "rho.csv";
	std::ofstream rhofs{ rhoname };
	rhofs << "rho_randseed:" <<gen.rho_randseed << std::endl;
	auto rhooutput = [&](Cell& cell, bool last) {
		rhofs << cell.rho;

		if (!last) {
			rhofs << std::endl;
		}
	};
	gen.cell_foreach_withend(rhooutput);


	//velocity data
	std::string veloname = "velocity.csv";
	std::ofstream velofs{ veloname };
	velofs << "velo_randseed" << gen.velo_randseed << std::endl;
	auto velooutput = [&](Cell& cell, bool last) {
		velofs << std::get<0>(cell.velocity) << "," << std::get<1>(cell.velocity) << "," << std::get<0>(cell.velocity);

		if (!last) {
			velofs << std::endl;
		}
	};
	gen.cell_foreach_withend(velooutput);



	//dual-mesh twenty vertices datas
	std::array<std::array<std::tuple<int, int, int>,27>, 8> duals;
	duals[0] = gen.gettwn(2, 2, 2);
	duals[1] = gen.gettwn(4, 2, 2);
	duals[2] = gen.gettwn(4, 4, 2);
	duals[3] = gen.gettwn(2, 4, 2);
	duals[4] = gen.gettwn(2, 2, 4);
	duals[5] = gen.gettwn(2, 4, 4);
	duals[6] = gen.gettwn(4, 4, 4);
	duals[7] = gen.gettwn(4, 2, 4);

	std::ofstream dualtwnfs{ "dualstwn.csv" };
	for (auto i = duals.begin(), e = duals.end(); i != e; i++) {
		for (auto twnb = 0, size = 20; twnb < size; twnb++) {
			auto pos = (*i)[twnb];
			dualtwnfs << gen.getindex(std::get<0>(pos),std::get<1>(pos),std::get<2>(pos)) + 1;
			if (twnb + 1 != size) {
				dualtwnfs << ",";
			}
		}

		if (i + 1 != e) {
			dualtwnfs << std::endl;
		}
	}

	std::string divrhoname = "divrho.csv";
	std::ofstream divrhofs{ divrhoname };
	auto divrhooutput = [&](Cell& cell, bool last) {
		divrhofs << std::get<0>(cell.grad_rho) << "," << std::get<1>(cell.grad_rho) << "," << std::get<2>(cell.grad_rho);

		if (!last) {
			divrhofs << std::endl;
		}
	};
	gen.cell_foreach_withend(divrhooutput);
	

	//rho data
	std::string divveloname = "divvelocity.csv";
	std::ofstream divvelofs{ divveloname };
	divvelofs << "divvelo_randseed:" << gen.velo_randseed << std::endl;
	auto divvelooutput = [&](Cell& cell, bool last) {
		divvelofs << cell.div_velocity;

		if (!last) {
			divvelofs << std::endl;
		}
	};
	gen.cell_foreach_withend(divvelooutput);


}


