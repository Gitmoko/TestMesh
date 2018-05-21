
#include <iostream>
#include <cstdlib>
#include <map>
#include <array>
#include <string>
#include <vector>
#include <fstream>
#include <random>
//-0.5,-0.5,0.5の位置をNo.1ノードとする


struct Vec {
	enum class Type :int {Corner,Surface,Edge,Center};
	double x = 0, y = 0, z = 0;
	Type type;

	std::string show() {
		std::string ret;
		ret += std::to_string(x);
		ret += ",";
		ret += std::to_string(y);
		ret += ",";
		ret += std::to_string(z);
		return ret;
	}

	Vec operator+(const Vec& pos) {
		auto ret = Vec{ x + pos.x,y + pos.y,z + pos.z };
		return ret;
	}

	Vec operator*(float c) {
		auto ret = Vec{ x*c,y*c,z*c };
		return ret;
	}
};

#define P (1)
#define Z (0)
#define M (-1)

std::array<Vec, 28>& getRelPos() {
	static std::array<Vec, 28> RelPos;
	RelPos[0] = { Z,Z,Z,Vec::Type::Center};
	RelPos[1] = { M,M,M,Vec::Type::Corner };
	RelPos[2] = { P,M,M,Vec::Type::Corner };
	RelPos[3] = { P,P,M,Vec::Type::Corner };
	RelPos[4] = { M,P,M,Vec::Type::Corner };
	RelPos[5] = { M,M,P,Vec::Type::Corner };
	RelPos[6] = { P,M,P,Vec::Type::Corner };
	RelPos[7] = { P,P,P ,Vec::Type::Corner };
	RelPos[8] = { M,P,P,Vec::Type::Corner };
	RelPos[9] = { Z,M,M,Vec::Type::Edge };
	RelPos[10] = { P,Z,M,Vec::Type::Edge };
	RelPos[11] = { Z,P,M,Vec::Type::Edge };
	RelPos[12] = { M,Z,M,Vec::Type::Edge };
	RelPos[13] = { Z,M,P,Vec::Type::Edge };
	RelPos[14] = { P,Z,P,Vec::Type::Edge };
	RelPos[15] = { Z,P,P,Vec::Type::Edge };
	RelPos[16] = { M,Z,P,Vec::Type::Edge };
	RelPos[17] = { M,M,Z ,Vec::Type::Edge };
	RelPos[18] = { P,M,Z ,Vec::Type::Edge };
	RelPos[19] = { P,P,Z ,Vec::Type::Edge };
	RelPos[20] = { M,P,Z ,Vec::Type::Edge };
	RelPos[21] = { Z,Z,M ,Vec::Type::Surface};
	RelPos[22] = { Z,Z,P ,Vec::Type::Surface };
	RelPos[23] = { Z,M,Z ,Vec::Type::Surface };
	RelPos[24] = { P,Z,Z ,Vec::Type::Surface };
	RelPos[25] = { Z,P,Z ,Vec::Type::Surface };
	RelPos[26] = { M,Z,Z ,Vec::Type::Surface };
	return RelPos;

}
void move_offset(std::array<Vec, 28>& RelPos,Vec pos) {
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



struct Edge {
	int node1;
	int node2;
};

struct Surface {
	std::array<Edge, 4> edges;
};


struct FormBase {
	enum class Type :int { Three, Two, One, Zero };
	int index;
	Type type;
	virtual ~FormBase() {}
};

struct ThreeForm:public FormBase {
	Cell origin_cell;
	int dual_point;
	double rho;
	virtual ~ThreeForm() override {}
};

struct TwoForm:public FormBase {
	Surface origin_surf;
	Edge dual_edge;
	std::tuple<double,double,double> vel;
	virtual ~TwoForm() override {}
};

struct OneForm:public FormBase {
	Edge origin_edge;
	Surface dual_sruface;
	virtual ~OneForm() override {}
};

struct ZeroForm:public FormBase {
	int origin_point;
	Cell dual_cell;
	virtual ~ZeroForm() override {}
};

//node index -> Form
//Func:: Int->Pos
std::map<int, FormBase*> createForms(std::map <std::tuple<int, int, int>, Cell>& cells, std::map<int, Vec>& Nodes) {
	std::map<int, FormBase*> ret;
	auto RelPos = getRelPos();
	for (auto& cell : cells) {
		for (int i = 0; i < cell.second.twn.size(); i++){
			auto type = RelPos[i].type;
			auto index = cell.second.twn[i];

			if (type == Vec::Type::Center) {
				if (ret.count(index) > 0) {
					auto ptr = new ThreeForm{};
					ptr->index = index;
					ptr->type = FormBase::Type::Three;
					ptr->dual_point = index;
					
				}
			}
			else if (type == Vec::Type::Corner) {
				if (ret.count(index) > 0) {
					auto ptr = new ZeroForm{};
					ptr->index = index;
					ptr->type = FormBase::Type::Three;
					ptr->dual_cell;
				}
			}
			else if (type == Vec::Type::Edge) {
				if (ret.count(index) > 0) {

				}
			}
			else if (type == Vec::Type::Surface) {
				if (ret.count(index) > 0) {

				}
			}
		}
		return ret;
	}
}

struct Generate {

	std::map<int, Vec> Nodes;
	std::map <std::tuple<int,int,int>,Cell> cells;
	double celllength;
	Vec offset = Vec{0,0,0};
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
		auto RelPos = getRelPos();
	
		for (int i = 1, e = 27; i <= e;i++) {
			std::tuple<int, int, int> n;
			std::get<0>(n) = x + RelPos[i].x;
			std::get<1>(n) = y + RelPos[i].y;
			std::get<2>(n) = z + RelPos[i].z;

			ret[i-1] = n;
		}
		return ret;
	}

	Vec getpos(int x, int y, int z) {
		Vec ret = offset + Vec{ celllength*x/2.0 , celllength*y/2.0 , celllength*z/2.0 };
		return ret;
	}

	Generate(double celllength_) :celllength(celllength_) {
		std::random_device randdev;
		rho_randseed = randdev();
		velo_randseed = randdev();
		std::mt19937 rhoengine(rho_randseed);
		std::mt19937 veloengine(velo_randseed);
		std::uniform_real_distribution<> dist(-5.0, 5.0);

		auto RelPos = getRelPos();
		move_offset(RelPos, Vec{ 1,1,1 });

		int index = 0;
		for (auto x = 0; x < xnum; x++) {
			for (auto y = 0; y < ynum; y++) {
				for (auto z = 0; z < znum; z++) {
					
					std::array<int, 27> twn;
					for (auto i = 1, e = 27; i <= e;i++) {
						auto& relpos = RelPos[i];
						auto node_relpos = Vec{x*2.0,y*2.0,z*2.0}+relpos;
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




void output(Generate& gen) {


	//node data
	std::ofstream nodefs{ "nodes.csv" };
	for (auto i = gen.Nodes.begin(), e = gen.Nodes.end(); i != e; i++) {
		nodefs << i->second.show();
		if (std::next(i) != e) {
			nodefs << std::endl;
		}
	}


	//twenty vertices data
	std::ofstream twnfs{ "twn.csv" };
	auto twnoutput = [&](Cell& cell, bool last) {
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
	rhofs << "rho_randseed:" << gen.rho_randseed << std::endl;
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
	std::array<std::array<std::tuple<int, int, int>, 27>, 8> duals;
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
			dualtwnfs << gen.getindex(std::get<0>(pos), std::get<1>(pos), std::get<2>(pos)) + 1;
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

int main()
{
	Generate gen{2};
	output(gen);



}


