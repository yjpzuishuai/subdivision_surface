#include<iostream>
#include"meshView.h"
#include"meshOperation.h"
using namespace NUBELib;


/*
得到网格中每个点的最大的度及其个数
*/
std::map<int, int, std::greater<int>> get_max_degree_and_count(M* pMesh)
{
	std::map<int, int, std::greater<int>> degree_count;
	int max_degree = 0;
	for (auto vertex : pMesh->vertices())
	{
		max_degree = std::max(max_degree, vertex->degree());
		degree_count[vertex->degree()]++;
	}
	return degree_count;
}


/*
红绿细分测试
*/
void test250724(int argc, char* argv[])
{
	M trimesh;

	trimesh.read_obj(argv[1]);
	std::map<int, bool> mask;
	std::map<int, bool> mask1;
	int n = 10;
	std::vector<M> trimeshs1(11);
	std::vector<M> trimeshs2(11);
	compute_vertex_degree(&trimesh);
	std::map<int, int, std::greater<int>> max_degree_and_count;
	max_degree_and_count = get_max_degree_and_count(&trimesh);
	std::cout << "细分前图像：" << std::endl;
	for (auto degree_and_count : max_degree_and_count)
		std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
	trimeshs1[0] = trimesh;
	trimeshs2[0] = trimesh;
	for (int i = 0;i < 5;i++)
	{
		std::queue<M::V*> q;
		std::map<int, bool> p_ed;//点是否被扩展过
		q.push(trimeshs1[i].vertices().front());
		int times = 100*i;
		while (!q.empty() && times--)
		{
			auto point = q.front();
			q.pop();
			if (p_ed.count(point->id()))continue;
			p_ed[point->id()] = true;
			for (M::It::VCcwOutHEIterator it(&trimeshs1[i], point);it != it.end();++it)
			{
				auto he = *it;
				q.push(trimeshs1[i].halfedgeTarget(he));
				if (mask.count(he->face()->id()))continue;
				mask[he->face()->id()] = true;
			}
		}
		std::queue<M::V*> q1;
		std::map<int, bool> p_ed1;//点是否被扩展过
		q1.push(trimeshs2[i].vertices().front());
		times = 100*i;
		while (!q1.empty() && times--)
		{
			auto point = q1.front();
			q1.pop();
			if (p_ed1.count(point->id()))continue;
			p_ed1[point->id()] = true;
			for (M::It::VCcwOutHEIterator it(&trimeshs2[i], point);it != it.end();++it)
			{
				auto he = *it;
				q1.push(trimeshs2[i].halfedgeTarget(he));
				if (mask1.count(he->face()->id()))continue;
				mask1[he->face()->id()] = true;
			}
		}
		red_green_division(&trimeshs1[i], &trimeshs1[i + 1], mask);
		red_green_division_pro(&trimeshs2[i], &trimeshs2[i + 1], mask1);
		compute_vertex_degree(&trimeshs1[i + 1]);
		std::map<int, int, std::greater<int>> max_degree_and_count;
		max_degree_and_count = get_max_degree_and_count(&trimeshs1[i + 1]);
		int point_count = 0;
		int degree_count = 0;
		std::cout << "red_green_division第" << i + 1 << "次细分" << std::endl;
		for (auto degree_and_count : max_degree_and_count)
		{
			std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
			point_count += degree_and_count.second;
			degree_count += degree_and_count.first * degree_and_count.second;
		}
		std::cout << "总点数： " << point_count << std::endl;
		std::cout << "点平均度数： " << 1.0 * degree_count / point_count << std::endl;
		std::cout << std::endl;
		compute_vertex_degree(&trimeshs2[i + 1]);
		max_degree_and_count = get_max_degree_and_count(&trimeshs2[i + 1]);
		std::cout << "red_green_division_pro第" << i + 1 << "次细分" << std::endl;
		point_count = 0;
		degree_count = 0;
		for (auto degree_and_count : max_degree_and_count)
		{
			std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
			point_count += degree_and_count.second;
			degree_count += degree_and_count.first * degree_and_count.second;
		}
		std::cout << "总点数： " << point_count << std::endl;
		std::cout << "点平均度数： " << 1.0 * degree_count / point_count << std::endl;
		std::cout << std::endl;
	}


	viewMesh(&trimeshs2[4]);

	
}
void test250725_1(int argc, char* argv[])
{
	M trimesh;
	int n = 10;
	std::vector<M> trimeshs1(11);
	std::vector<M> trimeshs2(11);
	trimesh.read_obj(argv[1]);
	compute_vertex_degree(&trimesh);
	std::map<int, int, std::greater<int>> max_degree_and_count;
	max_degree_and_count = get_max_degree_and_count(&trimesh);
	std::cout << "细分前图像：" << std::endl;
	for(auto degree_and_count: max_degree_and_count)
		std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
	trimeshs1[0] = trimesh;
	trimeshs2[0] = trimesh;
	std::map<int, bool> mask;
	/*for(auto face:trimesh)*/
	mask[1] = true;
	mask[2] = true;
	for (int i = 0;i < 6;i++)
	{
		red_green_division(&trimeshs1[i], &trimeshs1[i+1], mask);
		red_green_division_pro(&trimeshs2[i], &trimeshs2[i + 1], mask);
		compute_vertex_degree(&trimeshs1[i+1]);
		std::map<int, int, std::greater<int>> max_degree_and_count;
		max_degree_and_count = get_max_degree_and_count(&trimeshs1[i + 1]);
		int point_count = 0;
		int degree_count = 0;
		std::cout << "red_green_division第" << i + 1 << "次细分" << std::endl; 
		for (auto degree_and_count : max_degree_and_count)
		{
			std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
			point_count += degree_and_count.second;
			degree_count += degree_and_count.first * degree_and_count.second;
		}
		std::cout << "总点数： " << point_count << std::endl;
		std::cout << "点平均度数： " << 1.0*degree_count/point_count << std::endl;
		std::cout << std::endl;
		compute_vertex_degree(&trimeshs2[i + 1]);
		max_degree_and_count = get_max_degree_and_count(&trimeshs2[i + 1]);
		std::cout << "red_green_division_pro第" << i + 1 << "次细分" << std::endl;
		point_count = 0;
		degree_count = 0;
		for (auto degree_and_count : max_degree_and_count)
		{
			std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
			point_count += degree_and_count.second;
			degree_count += degree_and_count.first * degree_and_count.second;
		}
		std::cout << "总点数： " << point_count << std::endl;
		std::cout << "点平均度数： " << 1.0 * degree_count / point_count << std::endl;
		std::cout << std::endl;
	}


	viewMesh(&trimeshs1[6]);
}

/*
loop细分测试
*/
void test250724_2(int argc, char* argv[])
{
	M trimesh;
	std::vector<M> trimeshs(6);
	trimesh.read_obj(argv[1]);
	trimeshs[0] = trimesh;
	for (int i = 0;i < 5;i++)
	{
		std::cout << "正在进行第" << i + 1 << "次" << std::endl;
		loop_subdivision(&trimeshs[i], &trimeshs[i+1]);
		std::string output_filename = "../data/output/loop_" + std::to_string(i + 1) + ".obj";
		trimeshs[i + 1].write_obj(output_filename.c_str());
	}
	viewMesh(&trimeshs[trimeshs.size()-1]);
}
/*
改进butterfly算法测试
*/
void test250724_3(int argc, char* argv[])
{
	//局部细分
	M trimesh;
	trimesh.read_obj(argv[1]);
	std::map<int, bool> mask;
	std::vector<M> trimeshs1(11);
	compute_vertex_degree(&trimesh);
	std::map<int, int, std::greater<int>> max_degree_and_count;
	max_degree_and_count = get_max_degree_and_count(&trimesh);
	std::cout << "细分前图像：" << std::endl;
	for (auto degree_and_count : max_degree_and_count)
		std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
	trimeshs1[0] = trimesh;
	for (int i = 0;i < 4;i++)
	{
		std::queue<M::V*> q;
		std::map<int, bool> p_ed;//点是否被扩展过
		q.push(trimeshs1[i].vertices().front());
		int times = 100 * i;
		while (!q.empty() && times--)
		{
			auto point = q.front();
			q.pop();
			if (p_ed.count(point->id()))continue;
			p_ed[point->id()] = true;
			for (M::It::VCcwOutHEIterator it(&trimeshs1[i], point);it != it.end();++it)
			{
				auto he = *it;
				q.push(trimeshs1[i].halfedgeTarget(he));
				if (mask.count(he->face()->id()))continue;
				mask[he->face()->id()] = true;
			}
		}
		butterfly_division_pro(&trimeshs1[i], &trimeshs1[i + 1], mask);
		compute_vertex_degree(&trimeshs1[i + 1]);
		std::map<int, int, std::greater<int>> max_degree_and_count;
		max_degree_and_count = get_max_degree_and_count(&trimeshs1[i + 1]);
		int point_count = 0;
		int degree_count = 0;
		std::cout << "butterfly_division_pro第" << i + 1 << "次细分" << std::endl;
		for (auto degree_and_count : max_degree_and_count)
		{
			std::cout << "degree: " << degree_and_count.first << "degree count: " << degree_and_count.second << std::endl;
			point_count += degree_and_count.second;
			degree_count += degree_and_count.first * degree_and_count.second;
		}
		std::cout << "总点数： " << point_count << std::endl;
		std::cout << "点平均度数： " << 1.0 * degree_count / point_count << std::endl;
		std::cout << std::endl;

	}
	viewMesh(&trimeshs1[4]);

	//全局细分
	/*M trimesh;
	M trimesh2;
	M trimesh3;
	std::vector<M> trimeshs(3);
	
	trimesh.read_obj(argv[1]);
	trimeshs[0] = trimesh;
	for (int i = 0;i < 2 ;i++)
	{
		std::map<int, bool> mask;
		for (auto face : trimeshs[i].faces())
		{
			mask[face->id()] = true;
		}
		std::cout << "正在进行第" << i + 1 << "次" << std::endl;
		butterfly_division_pro(&trimeshs[i], &trimeshs[i + 1], mask);
	}
	viewMesh(&trimeshs[trimeshs.size() - 1]);*/
}
void test()
{
	std::map<int, int> degree_count;
	degree_count[0]++;
	std::cout << degree_count[1] << std::endl;
}
int main(int argc, char* argv[])
{
	//test250724(argc, argv);
	//test250724_2(argc, argv);
	test250724_3(argc, argv);
	//test();
	//test250725_1(argc, argv);
}
