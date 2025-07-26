#include "meshOperation.h"
#include"..\3rdParty\MeshLib\core\Parser\parser.h"
#include<queue>
#define M_PI 3.14159265358979323846

void NUBELib::loop_subdivision(M* pMesh,M* pMesh1)
{
	std::map<int, int> edge_vertex;
	int vid = 1;
	int fid = 1;
	for (auto ver : pMesh->vertices())
	{
		std::vector<M::V*> vccwv;
		std::vector<M::V*> prev_vccwv;
		int degree = 0;
		for (M::It::VCcwOutHEIterator it(pMesh, ver);it != it.end();++it)
		{
			degree++;
			auto he = *it;
			prev_vccwv.push_back(pMesh->halfedgeTarget(he));
			CPoint p =  get_point_on_halfedge(he);
			int edge_id = he->edge()->id();
			if (edge_vertex.count(edge_id))
			{
				vccwv.push_back(pMesh1->idVertex(edge_vertex[edge_id]));
			}
			else
			{
				M::V* vertex = pMesh1->createVertex(vid);
				vertex->point() = p;
				edge_vertex[edge_id] = vid;
				vid++;
				vccwv.push_back(vertex);
			}
		}
		double beta = 1.0 * (5.0 / 8 - (3.0 / 8 + 1.0 / 4 * cos(M_PI * 2 / degree)) * (3.0 / 8 + 1.0 / 4 * cos(M_PI * 2 / degree))) / degree;
		CPoint new_center(0,0,0);
		for (int i = 0;i < degree;i++)
		{
			new_center += prev_vccwv[i]->point() * beta;

		}
		new_center += ver->point() * (1.0 - beta * degree);
		M::V* vertex = pMesh1->createVertex(vid);
		vid++;
		vertex->id() = vid;
		vertex->point() = new_center;
		vccwv.push_back(vccwv[0]);
		//建面
		for (int i = 0;i < vccwv.size() - 1;i++)
		{
			std::vector<M::V*> new_face_points;
			new_face_points.push_back(vccwv[i]);
			new_face_points.push_back(vccwv[i+1]);
			new_face_points.push_back(vertex);
			pMesh1->createFace(new_face_points, fid);
			fid++;
		}
	}
	for (auto face : pMesh->faces())
	{
		std::vector<M::V*> newpoints;
		auto he = face->halfedge();
		int edge_id1 = he->edge()->id();
		int edge_id2 = he->he_next()->edge()->id();
		int edge_id3 = he->he_next()->he_next()->edge()->id();
		newpoints.push_back(pMesh1->idVertex(edge_vertex[edge_id1]));
		newpoints.push_back(pMesh1->idVertex(edge_vertex[edge_id2]));
		newpoints.push_back(pMesh1->idVertex(edge_vertex[edge_id3]));
		pMesh1->createFace(newpoints, fid);
		fid++;
	}
	pMesh1->labelBoundary();
}


CPoint NUBELib::get_point_on_halfedge(M::HE* he)
{
	auto he1 = he->he_sym();
	CPoint a = he->target()->point();
	CPoint b = he1->target()->point();
	CPoint c = he->he_next()->target()->point();
	CPoint d = he1->he_next()->target()->point();
	return (a + b) * 3 / 8 + (c + d) * 1 / 8;
}

void NUBELib::red_green_division(M* pMesh, M* pMesh1, std::map<int, bool>& faces_to_refine)
{
	//旧边->pMesh1中点
	std::map<typename M::E*, typename M::V*> edge_to_midpoint_vertex;
	//旧顶点->pMesh1新顶点
	std::map<typename M::V*, typename M::V*> old_v_to_new_v;
	int new_vid = 1;
	int new_fid = 1;
	//复制所有旧顶点到新网格,建立点到点映射
	for (auto v_old : pMesh->vertices())
	{
		M::V* new_v = pMesh1->createVertex(new_vid++);
		new_v->point() = v_old->point();
		new_v->id() = v_old->id();
		old_v_to_new_v[v_old] = new_v;
	}
	//根据mask内容区分需要细分的面，在这些面中创建出每条边的中点，并建立边->顶点映射
	for (auto const face_to_refine : faces_to_refine)
	{
		if (!face_to_refine.second) continue;
		M::F* f_old = pMesh->idFace(face_to_refine.first);
		if (!f_old) continue;
		M::HE* he = pMesh->faceHalfedge(f_old);
		do 
		{
			M::E* edge = pMesh->halfedgeEdge(he);
			if (edge_to_midpoint_vertex.count(edge)) 
			{
				he = pMesh->halfedgeNext(he);
				continue;
			}
			M::V* v1 = pMesh->edgeVertex1(edge);
			M::V* v2 = pMesh->edgeVertex2(edge);
			CPoint new_point = (v1->point() + v2->point()) * 0.5;
			while (pMesh1->idVertex(new_vid))
			{
				new_vid++;
			}
			M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);
			new_midpoint_v->point() = new_point;
			edge_to_midpoint_vertex[edge] = new_midpoint_v;
			he = pMesh->halfedgeNext(he);
		} while (he != pMesh->faceHalfedge(f_old));
	}
	//遍历每个面，根据每个面所有边上的中点数量进行分类
	/*
	1.case0:面的每条边无中点，直接建面
	2.case1：1个中点，二分
	3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
	4.case3：正常1分为4
	*/
	for (auto f_old : pMesh->faces())
	{
		std::vector<M::V*> old_vertices;
		M::HE* he_start = pMesh->faceHalfedge(f_old);
		M::HE* he_iter = he_start;
		do 
		{
			old_vertices.push_back(pMesh->halfedgeTarget(he_iter));
			he_iter = pMesh->halfedgeNext(he_iter);
		} while (he_iter != he_start);
		std::vector<M::E*> old_edges;
		old_edges.push_back(pMesh->vertexEdge(old_vertices[0], old_vertices[1]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[1], old_vertices[2]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[2], old_vertices[0]));
		std::vector<M::V*> midpoints;
		std::vector<int> split_edge_indices;
		for (int i = 0; i < 3; ++i) 
		{
			if (edge_to_midpoint_vertex.count(old_edges[i])) 
			{
				midpoints.push_back(edge_to_midpoint_vertex.at(old_edges[i]));
				split_edge_indices.push_back(i);
			}
		}
		int hanging_node_count = midpoints.size();
		switch (hanging_node_count)
		{
			//1.case0:面的每条边无中点，直接建面
			case 0:
			{
				pMesh1->createFace(std::vector<M::V*>{ old_v_to_new_v.at(old_vertices[0]), old_v_to_new_v.at(old_vertices[1]), old_v_to_new_v.at(old_vertices[2]) }, new_fid++);
				break;
			}
			//2.case1：1个中点，二分
			case 1:
			{
				M::V* m = midpoints[0];
				int edge_idx = split_edge_indices[0];
				M::V* v_start = old_v_to_new_v.at(old_vertices[edge_idx]);
				M::V* v_end = old_v_to_new_v.at(old_vertices[(edge_idx + 1) % 3]);
				M::V* v_opposite = old_v_to_new_v.at(old_vertices[(edge_idx + 2) % 3]);

				pMesh1->createFace(std::vector<M::V*>{ v_start, m, v_opposite }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ m, v_end, v_opposite }, new_fid++);
				break;
			}
			//3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
			case 2:
			{
				bool edge_flip = true;
				M::V* m1 = midpoints[0];
				M::V* m2 = midpoints[1];
				int e_idx1 = split_edge_indices[0];
				int e_idx2 = split_edge_indices[1];
				M::V* v_common_old, * v_A_old, * v_B_old;
				//保持逆时针，正方向与原来一致
				if (e_idx1 == 0 && e_idx2 == 1) 
				{
					v_common_old = old_vertices[1],v_A_old = old_vertices[0], v_B_old = old_vertices[2];
				}
				else if (e_idx1 == 1 && e_idx2 == 2) 
				{
					v_common_old = old_vertices[2], v_A_old = old_vertices[1], v_B_old = old_vertices[0];
				}
				else if(e_idx1 == 0 && e_idx2 == 2)
				{
					//保持逆时针，正方向与原来一致
					v_common_old = old_vertices[0], v_A_old = old_vertices[2], v_B_old = old_vertices[1];
					/*v_common_old = old_vertices[0], v_A_old = old_vertices[1], v_B_old = old_vertices[2];*/
				}
				
				M::V* v_common_new = old_v_to_new_v.at(v_common_old);
				M::V* v_A_new = old_v_to_new_v.at(v_A_old);
				M::V* v_B_new = old_v_to_new_v.at(v_B_old);
				M::V* m_A, * m_B;
				M::E* edge_A_common_old = pMesh->vertexEdge(v_A_old, v_common_old);

				if (edge_to_midpoint_vertex.at(edge_A_common_old) == m1) 
				{
					m_A = m1; m_B = m2;
				}
				else 
				{
					//e_idx1 == 0 && e_idx2 == 2 特殊情况处理
					m_A = m2; m_B = m1;
				}
				//是否加入边翻转
				if (edge_flip)
				{
					double A_r1 = get_circleRadius(m_A->point(), v_A_new->point(), m_B->point());
					double A_r2 = get_circleRadius(m_B->point(), v_A_new->point(), v_B_new->point());
					double max_A_r = std::max(A_r1, A_r2);
					double B_r1 = get_circleRadius(m_A->point(), v_B_new->point(), m_B->point());
					double B_r2 = get_circleRadius(m_A->point(), v_B_new->point(), v_A_new->point());
					double max_B_r = std::max(B_r1, B_r2);
					pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
					if (max_A_r <= max_B_r)
					{
						pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
						pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
					}
					else
					{
						pMesh1->createFace(std::vector<M::V*>{ v_B_new, m_A, m_B }, new_fid++);
						pMesh1->createFace(std::vector<M::V*>{ v_B_new, v_A_new,m_A  }, new_fid++);
					}
				}
				else
				{
					pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
				}
				break;
			}
			//4.case3：正常1分为4
			case 3:
			{
				M::V* v0_new = old_v_to_new_v.at(old_vertices[0]);
				M::V* v1_new = old_v_to_new_v.at(old_vertices[1]);
				M::V* v2_new = old_v_to_new_v.at(old_vertices[2]);
				M::V* m01_new = edge_to_midpoint_vertex.at(old_edges[0]);
				M::V* m12_new = edge_to_midpoint_vertex.at(old_edges[1]);
				M::V* m20_new = edge_to_midpoint_vertex.at(old_edges[2]);

				pMesh1->createFace(std::vector<M::V*>{ v0_new, m01_new, m20_new }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v1_new, m12_new, m01_new }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v2_new, m20_new, m12_new }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ m01_new, m12_new, m20_new }, new_fid++);
				break;
			}
		}
	}
	pMesh1->labelBoundary();
}

void NUBELib::red_green_division_pro(M* pMesh, M* pMesh1, std::map<int, bool>& faces_to_refine)
{
	//旧边->pMesh1中点
	std::map<M::E*, M::V*> edge_to_midpoint_vertex;
	//旧顶点->pMesh1新顶点
	std::map<M::V*, M::V*> old_v_to_new_v;
	//face->最长边
	std::map<M::F*, M::E*> face_to_maxlength_edge;
	int new_vid = 1;
	int new_fid = 1;
	//复制所有旧顶点到新网格,建立点到点映射
	for (auto v_old : pMesh->vertices())
	{
		M::V* new_v = pMesh1->createVertex(new_vid++);
		new_v->point() = v_old->point();
		old_v_to_new_v[v_old] = new_v;
	}
	//根据mask内容区分需要细分的面，在这些面中创建出每条边的中点，并建立边->顶点映射
	for (auto const face_to_refine : faces_to_refine)
	{
		if (!face_to_refine.second) continue;
		M::F* f_old = pMesh->idFace(face_to_refine.first);
		if (!f_old) continue;
		M::HE* he = pMesh->faceHalfedge(f_old);
		do
		{
			M::E* edge = pMesh->halfedgeEdge(he);
			if (edge_to_midpoint_vertex.count(edge))
			{
				he = pMesh->halfedgeNext(he);
				continue;
			}
			M::V* v1 = pMesh->edgeVertex1(edge);
			M::V* v2 = pMesh->edgeVertex2(edge);
			CPoint new_point = (v1->point() + v2->point()) * 0.5;
			while (pMesh1->idVertex(new_vid))
			{
				new_vid++;
			}
			M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);
			new_midpoint_v->point() = new_point;
			edge_to_midpoint_vertex[edge] = new_midpoint_v;
			he = pMesh->halfedgeNext(he);
		} while (he != pMesh->faceHalfedge(f_old));
	}

	for (auto f_old : pMesh->faces())
	{
		auto he = pMesh->faceHalfedge(f_old);
		auto he_start = he;
		auto max_he = he;
		double max_length = 0;
		do
		{
			double length = get_triangle_edge_length(he->edge());
			if (length > max_length)
			{
				max_he = he;
				max_length = length;
			}
			he =pMesh->halfedgeNext(he);
		} while (he != he_start);
		face_to_maxlength_edge[f_old] = pMesh->halfedgeEdge(max_he);
	}
	/*
	最长边细分优化：
	* 这里我们采用两种方法
	* 方法1：在满足三角形扭曲不是很大的情况下，尽量不加更多点和边
	情况1（最小角对应的边 & 次小角对应的边）：对于一条边（边上带中点）对应的角度如果<30°，将最长边的中点加入
	情况2（最大角对应的边）：直接连
	* 方法2：在满足三角形扭曲不是很大的情况下，尽量加更多点和边(两次方法1可以达到方法2效果)
	情况1（最小角对应的边 & 次小角对应的边）：对于一条边（无论边上带不带中点）对应的角度如果<30°，将最长边的中点加入
	*/
	//模式选择
	int subdivision_method = 1;
	for (auto f_old : pMesh->faces())
	{
		std::queue<M::F*> q_f;
		q_f.push(f_old);
		while(!q_f.empty())
		{
			auto face = q_f.front();
			q_f.pop();
			if (edge_to_midpoint_vertex.count(face_to_maxlength_edge[face]))continue;
			std::vector<M::V*> old_vertices;
			M::HE* he_start = pMesh->faceHalfedge(face);
			M::HE* he_iter = he_start;
			do
			{
				old_vertices.push_back(pMesh->halfedgeTarget(he_iter));
				he_iter = pMesh->halfedgeNext(he_iter);
			} while (he_iter != he_start);
			std::vector<M::E*> old_edges;
			old_edges.push_back(pMesh->vertexEdge(old_vertices[0], old_vertices[1]));
			old_edges.push_back(pMesh->vertexEdge(old_vertices[1], old_vertices[2]));
			old_edges.push_back(pMesh->vertexEdge(old_vertices[2], old_vertices[0]));
			std::vector<M::V*> midpoints;
			std::vector<int> split_edge_indices;
			for (int i = 0; i < 3; ++i)
			{
				if (edge_to_midpoint_vertex.count(old_edges[i]))
				{
					midpoints.push_back(edge_to_midpoint_vertex.at(old_edges[i]));
					split_edge_indices.push_back(i);
				}
			}
			int hanging_node_count = midpoints.size();
			if (subdivision_method == 0 || subdivision_method == 1 && hanging_node_count == 2)
			{
				int max_flag = 0, min_flag = 0;
				if (hanging_node_count == 0 || hanging_node_count == 3)continue;
				for (int i = 1;i < old_edges.size();i++)
				{
					if (pMesh->edgeLength(old_edges[i]) > pMesh->edgeLength(old_edges[max_flag])) { max_flag = i; continue; }
					if (pMesh->edgeLength(old_edges[i]) < pMesh->edgeLength(old_edges[min_flag]))min_flag = i;
				}
				auto he_min = pMesh->edgeHalfedge(old_edges[min_flag], 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(old_edges[min_flag], 0) : pMesh->edgeHalfedge(old_edges[min_flag], 1);
				double he_min_angle = get_div_halfedge_angle(he_min);
				//he_angle_max为角度弧度制的最大值（角度最小值），不能更大
				double he_angle_max = cos(M_PI / 6.1);
				if (he_min_angle > he_angle_max)
				{
					if (!edge_to_midpoint_vertex.count(old_edges[max_flag]))
					{
						while (pMesh1->idVertex(new_vid))
						{
							new_vid++;
						}
						M::V* v1 = pMesh->edgeVertex1(old_edges[max_flag]);
						M::V* v2 = pMesh->edgeVertex2(old_edges[max_flag]);
						CPoint new_point = (v1->point() + v2->point()) * 0.5;
						M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);
						new_midpoint_v->point() = new_point;
						edge_to_midpoint_vertex[old_edges[max_flag]] = new_midpoint_v;
						auto he_max = pMesh->edgeHalfedge(old_edges[max_flag], 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(old_edges[max_flag], 0) : pMesh->edgeHalfedge(old_edges[max_flag], 1);
						if (pMesh->halfedgeEdge(he_max)->boundary())continue;
						else q_f.push(pMesh->halfedgeFace(pMesh->halfedgeSym(he_max)));
					}
				}
			}
			else if (subdivision_method == 1)
			{
				switch (hanging_node_count)
				{
				case 1:
				{
					auto old_edge = old_edges[split_edge_indices[0]];
					auto he = pMesh->edgeHalfedge(old_edge, 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(old_edge, 0) : pMesh->edgeHalfedge(old_edge, 1);
					double he_angle = get_div_halfedge_angle(he);
					//he_angle_max为角度弧度制的最大值（角度最小值），不能更大
					double he_angle_max = cos(M_PI / 6);
					if (he_angle > he_angle_max)
					{
						M::E* max_e = face_to_maxlength_edge[face];
						M::HE* max_he = pMesh->edgeHalfedge(max_e, 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(max_e, 0) : pMesh->edgeHalfedge(max_e, 1);
						M::V* v1 = pMesh->halfedgeSource(max_he);
						M::V* v2 = pMesh->halfedgeTarget(max_he);
						CPoint new_point = (v1->point() + v2->point()) * 0.5;
						while (pMesh1->idVertex(new_vid))
						{
							new_vid++;
						}
						M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);
						new_midpoint_v->point() = new_point;
						edge_to_midpoint_vertex[max_e] = new_midpoint_v;
						if (pMesh->halfedgeEdge(max_he)->boundary())continue;
						else q_f.push(pMesh->halfedgeFace(pMesh->halfedgeSym(max_he)));
					}
					break;
				}
				//case 2:
				//{
				//	auto old_edge = old_edges[split_edge_indices[0]];
				//	auto he = pMesh->edgeHalfedge(old_edge, 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(old_edge, 0) : pMesh->edgeHalfedge(old_edge, 1);
				//	double he_angle = get_div_halfedge_angle(he);
				//	//he_angle_max为角度弧度制的最大值（角度最小值），不能更大
				//	double he_angle_max = cos(M_PI / 6);
				//	if (he_angle > he_angle_max)
				//	{
				//		M::E* max_e = face_to_maxlength_edge[face];
				//		M::HE* max_he = pMesh->edgeHalfedge(max_e, 0)->face()->id() == face->id() ? pMesh->edgeHalfedge(max_e, 0) : pMesh->edgeHalfedge(max_e, 1);
				//		M::V* v1 = pMesh->halfedgeSource(max_he);
				//		M::V* v2 = pMesh->halfedgeTarget(max_he);
				//		CPoint new_point = (v1->point() + v2->point()) * 0.5;
				//		while (pMesh1->idVertex(new_vid))
				//		{
				//			new_vid++;
				//		}
				//		M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);
				//		new_midpoint_v->point() = new_point;
				//		edge_to_midpoint_vertex[max_e] = new_midpoint_v;
				//		if (pMesh->halfedgeEdge(max_he)->boundary())continue;
				//		else q_f.push(pMesh->halfedgeFace(pMesh->halfedgeSym(max_he)));
				//	}
				//	break;
				//}
				}
				
			}
		}
	}
	//遍历每个面，根据每个面所有边上的中点数量进行分类
	/*
	1.case0:面的每条边无中点，直接建面
	2.case1：1个中点，二分
	3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
	4.case3：正常1分为4
	*/
	for (auto f_old : pMesh->faces())
	{
		std::vector<M::V*> old_vertices;
		M::HE* he_start = pMesh->faceHalfedge(f_old);
		M::HE* he_iter = he_start;
		do
		{
			old_vertices.push_back(pMesh->halfedgeTarget(he_iter));
			he_iter = pMesh->halfedgeNext(he_iter);
		} while (he_iter != he_start);
		std::vector<M::E*> old_edges;
		old_edges.push_back(pMesh->vertexEdge(old_vertices[0], old_vertices[1]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[1], old_vertices[2]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[2], old_vertices[0]));
		std::vector<M::V*> midpoints;
		std::vector<int> split_edge_indices;
		for (int i = 0; i < 3; ++i)
		{
			if (edge_to_midpoint_vertex.count(old_edges[i]))
			{
				midpoints.push_back(edge_to_midpoint_vertex.at(old_edges[i]));
				split_edge_indices.push_back(i);
			}
		}
		int hanging_node_count = midpoints.size();

		

		switch (hanging_node_count)
		{
		//1.case0:面的每条边无中点，直接建面
		case 0:
		{
			pMesh1->createFace(std::vector<M::V*>{ old_v_to_new_v.at(old_vertices[0]), old_v_to_new_v.at(old_vertices[1]), old_v_to_new_v.at(old_vertices[2]) }, new_fid++);
			break;
		}
		//2.case1：1个中点，二分
		case 1:
		{
			
			M::V* m = midpoints[0];
			int edge_idx = split_edge_indices[0];
			M::V* v_start = old_v_to_new_v.at(old_vertices[edge_idx]);
			M::V* v_end = old_v_to_new_v.at(old_vertices[(edge_idx + 1) % 3]);
			M::V* v_opposite = old_v_to_new_v.at(old_vertices[(edge_idx + 2) % 3]);

			pMesh1->createFace(std::vector<M::V*>{ v_start, m, v_opposite }, new_fid++);
			pMesh1->createFace(std::vector<M::V*>{ m, v_end, v_opposite }, new_fid++);
			break;
		}
		//3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
		case 2:
		{
			bool edge_flip = true;
			M::V* m1 = midpoints[0];
			M::V* m2 = midpoints[1];
			int e_idx1 = split_edge_indices[0];
			int e_idx2 = split_edge_indices[1];
			M::V* v_common_old, * v_A_old, * v_B_old;
			//保持逆时针，正方向与原来一致
			if (e_idx1 == 0 && e_idx2 == 1 || e_idx2 == 0 && e_idx1 == 1)
			{
				v_common_old = old_vertices[1], v_A_old = old_vertices[0], v_B_old = old_vertices[2];
			}
			else if (e_idx1 == 1 && e_idx2 == 2 || e_idx2 == 1 && e_idx1 == 2)
			{
				v_common_old = old_vertices[2], v_A_old = old_vertices[1], v_B_old = old_vertices[0];
			}
			else if (e_idx1 == 0 && e_idx2 == 2 || e_idx2 == 0 && e_idx1 == 2)
			{
				//保持逆时针，正方向与原来一致
				v_common_old = old_vertices[0], v_A_old = old_vertices[2], v_B_old = old_vertices[1];
				/*v_common_old = old_vertices[0], v_A_old = old_vertices[1], v_B_old = old_vertices[2];*/
			}

			M::V* v_common_new = old_v_to_new_v.at(v_common_old);
			M::V* v_A_new = old_v_to_new_v.at(v_A_old);
			M::V* v_B_new = old_v_to_new_v.at(v_B_old);
			M::V* m_A, * m_B;
			M::E* edge_A_common_old = pMesh->vertexEdge(v_A_old, v_common_old);

			if (edge_to_midpoint_vertex.at(edge_A_common_old) == m1)
			{
				m_A = m1; m_B = m2;
			}
			else
			{
				//e_idx1 == 0 && e_idx2 == 2 特殊情况处理
				m_A = m2; m_B = m1;
			}
			//是否加入边翻转
			if (edge_flip)
			{
				double A_r1 = get_circleRadius(m_A->point(), v_A_new->point(), m_B->point());
				double A_r2 = get_circleRadius(m_B->point(), v_A_new->point(), v_B_new->point());
				double max_A_r = std::max(A_r1, A_r2);
				double B_r1 = get_circleRadius(m_A->point(), v_B_new->point(), m_B->point());
				double B_r2 = get_circleRadius(m_A->point(), v_B_new->point(), v_A_new->point());
				double max_B_r = std::max(B_r1, B_r2);
				pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
				if (max_A_r <= max_B_r)
				{
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
				}
				else
				{
					pMesh1->createFace(std::vector<M::V*>{ v_B_new, m_A, m_B }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_B_new, v_A_new, m_A  }, new_fid++);
				}
			}
			else
			{
				pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
			}
			break;
		}
		//4.case3：正常1分为4
		case 3:
		{
			M::V* v0_new = old_v_to_new_v.at(old_vertices[0]);
			M::V* v1_new = old_v_to_new_v.at(old_vertices[1]);
			M::V* v2_new = old_v_to_new_v.at(old_vertices[2]);
			M::V* m01_new = edge_to_midpoint_vertex.at(old_edges[0]);
			M::V* m12_new = edge_to_midpoint_vertex.at(old_edges[1]);
			M::V* m20_new = edge_to_midpoint_vertex.at(old_edges[2]);
			//faces_to_refine[new_fid] = true;
			pMesh1->createFace(std::vector<M::V*>{ v0_new, m01_new, m20_new }, new_fid++);
			//faces_to_refine[new_fid] = true;
			pMesh1->createFace(std::vector<M::V*>{ v1_new, m12_new, m01_new }, new_fid++);
			//faces_to_refine[new_fid] = true;
			pMesh1->createFace(std::vector<M::V*>{ v2_new, m20_new, m12_new }, new_fid++);
			//faces_to_refine[new_fid] = true;
			pMesh1->createFace(std::vector<M::V*>{ m01_new, m12_new, m20_new }, new_fid++);
			break;
		}
		}
	}
	pMesh1->labelBoundary();

}
//r = abc/4/S
double NUBELib::get_circleRadius(CPoint a, CPoint b, CPoint c)
{
	CPoint ab = b - a;
	CPoint ac = c - a;
	CPoint bc = c - b;
	double len_ab = ab.norm();
	double len_ac = ac.norm();
	double len_bc = bc.norm();
	CPoint cross_product = ab ^ ac;
	double area_parallelogram_sq = cross_product.norm() * cross_product.norm();
	//判断是否共线
	const double epsilon = 1e-12;
	if (area_parallelogram_sq < epsilon) {
		return std::numeric_limits<double>::infinity();
	}
	double triangle_area = 0.5 * std::sqrt(area_parallelogram_sq);
	//R = (len_ab * len_ac * len_bc) / (4 * Area)     即abc/4/S
	double radius = (len_ab * len_ac * len_bc) / (4.0 * triangle_area);
	return radius;
}

void NUBELib::butterfly_division_pro(M* pMesh, M* pMesh1, std::map<int, bool>& faces_to_refine)
{
	compute_vertex_degree(pMesh);
	//旧边->pMesh1中点
	std::map<typename M::E*, typename M::V*> edge_to_midpoint_vertex;
	//旧顶点->pMesh1新顶点
	std::map<typename M::V*, typename M::V*> old_v_to_new_v;
	int new_vid = 1;
	int new_fid = 1;
	//复制所有旧顶点到新网格,建立点到点映射
	for (auto v_old : pMesh->vertices())
	{
		M::V* new_v = pMesh1->createVertex(new_vid++);
		new_v->point() = v_old->point();
		new_v->point() = v_old->point();
		old_v_to_new_v[v_old] = new_v;
	}
	//根据mask内容区分需要细分的面，在这些面中创建出每条边的中点，并建立边->顶点映射
	for (auto face_to_refine : faces_to_refine)
	{
		if (!face_to_refine.second) continue;
		M::F* f_old = pMesh->idFace(face_to_refine.first);
		if (!f_old) continue;
		M::HE* he = pMesh->faceHalfedge(f_old);
		do
		{
			M::E* edge = pMesh->halfedgeEdge(he);
			if (edge_to_midpoint_vertex.count(edge))
			{
				he = pMesh->halfedgeNext(he);
				continue;
			}
			M::V* v1 = pMesh->halfedgeSource(he);
			M::V* v2 = pMesh->halfedgeTarget(he);
			CPoint new_point(0, 0, 0);
			/*
			* 已解决
			边界情况不完整，缺少
			1.边不是边界边，点是边界点（目前用不到，先不考虑，遍历会中断）
			使用VCcwVIterator2 VClwVIterator2 第二个构造函数，分别遍历到边界停止，并把点按索引存储 优化非边界情况即可（已解决）
			*/
			if (edge->boundary())
			{
				////自动从边界半边开始遍历
				//M::It::VCcwVIterator it1(pMesh, v1);
				//M::It::VClwVIterator it2(pMesh, v2);

				//new_point += v1->point() * 9 / 16;
				//new_point += v2->point() * 9 / 16;
				//for (int i = 1;i < v1->degree();i++)
				//{
				//	++it1;
				//}
				//for (int i = 1;i < v2->degree();i++)
				//{
				//	++it2;
				//}
				//assert((*it1)->boundary());
				//assert((*it2)->boundary());
				//new_point -= (*it1)->point() / 16;
				//new_point -= (*it2)->point() / 16;
			}
			//情况2：边两端点度均为6，权重分别为
			/*
			* 组成：与共有此边的两面存在公共边的所有面上的点
				1.两端点：1/2，即以此边的半边开始做顺逆时针遍历分别得到的第一个点
				2.以此边的半边开始做顺逆时针遍历分别得到的第二个点（端点是第一个）
				3.以此边的半边开始做顺逆时针遍历分别得到的第三个点
			*/
			else if (v1->degree() == 6 && v2->degree() == 6)
			{
				M::It::VCcwVIterator2 it1(pMesh, v1, he);
				M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
				std::vector<CPoint> verts1(7);//v1周围点
				std::vector<CPoint> verts2(7);//v2周围点
				int begin_index = 0,end_index = 6;
				while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
				{
					if (!(*it1)->boundary())
					{
						verts1[begin_index++] = (*it1)->point();
						++it1;
					}
					if (!(*it2)->boundary())
					{
						verts1[end_index--] = (*it2)->point();
						++it2;
					}
				}
				it1 = M::It::VCcwVIterator2(pMesh, v2, pMesh->halfedgeSym(he));
				it2 = M::It::VClwVIterator2(pMesh, v2, he);
				begin_index = 0, end_index = 6;
				while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
				{
					if (!(*it1)->boundary())
					{
						verts2[begin_index++] = (*it1)->point();
						++it1;
					}
					if (!(*it2)->boundary())
					{
						verts2[end_index--] = (*it2)->point();
						++it2;
					}
				}
				new_point = (verts1[0] + verts2[0]) / 2 + (verts1[1] + verts2[1]) / 8 - (verts1[2] + verts1[4] + verts2[2] + verts2[4]) / 16;
			}
			else if (v1->degree() != 6 && v2->degree() == 6 && v1->degree() >= 3)
			{
				double w_v1 = 3.0 / 4;
				new_point += v1->point() * w_v1;
				switch (v1->degree())
				{
				case 3:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(4);//v1周围点
					int begin_index = 0, end_index = 3;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts1[0] * 5 / 12 - (verts1[1] + verts1[2]) / 12;
					break;
				}
				case 4:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(5);//v1周围点
					int begin_index = 0, end_index = 4;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts1[0] * 3 / 8 - (verts1[2]) / 8;
					break;
				}
				default:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(v1->degree()+1);//v1周围点
					int begin_index = 0, end_index = v1->degree();
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					for (int i = 0;i < v1->degree();i++)
					{
						double w_vi = 1.0 / v1->degree() * (1.0 / 4 + cos(2 * M_PI * i / v1->degree()) + 1.0 / 2 * cos(4 * M_PI * i / v1->degree()));
						new_point += verts1[i] * w_vi;
					}
					break;
				}
				}
			}
			else if (v1->degree() == 6 && v2->degree() != 6 && v2->degree() >= 3)
			{
				double w_v2 = 3.0 / 4;
				new_point += v2->point() * w_v2;
				switch (v2->degree())
				{
				case 3:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(4);//v2周围点
					int begin_index = 0, end_index = 3;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts2[0] * 5 / 12 - (verts2[1] + verts2[2]) / 12;
					break;
				}
				case 4:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(5);//v2周围点
					int begin_index = 0, end_index = 4;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts2[0] * 3 / 8 - (verts2[2]) / 8;
					break;
				}
				default:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(v2->degree()+1);//v2周围点
					int begin_index = 0, end_index = v2->degree();
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					for (int i = 0;i < v2->degree();i++)
					{
						double w_vi = 1.0 / v2->degree() * (1.0 / 4 + cos(2 * M_PI * i / v2->degree()) + 1.0 / 2 * cos(4 * M_PI * i / v2->degree()));
						new_point += verts2[i] * w_vi;
					}
					break;
				}
				}
			}
			else if (v1->degree() != 6 && v2->degree() != 6 && v1->degree() >= 3 && v2->degree() >= 3)
			{
				double w_v1 = 3.0 / 4;
				new_point += v1->point() * w_v1;
				switch (v1->degree())
				{
				case 3:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(4);//v1周围点
					int begin_index = 0, end_index = 3;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts1[0] * 5 / 12 - (verts1[1] + verts1[2]) / 12;
					break;
				}
				case 4:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(5);//v1周围点
					int begin_index = 0, end_index = 4;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts1[0] * 3 / 8 - (verts1[2]) / 8;
					break;
				}
				default:
				{
					M::It::VCcwVIterator2 it1(pMesh, v1, he);
					M::It::VClwVIterator2 it2(pMesh, v1, pMesh->halfedgeSym(he));
					std::vector<CPoint> verts1(v1->degree() + 1);//v1周围点
					int begin_index = 0, end_index = v1->degree();
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts1[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts1[end_index--] = (*it2)->point();
							++it2;
						}
					}
					for (int i = 0;i < v1->degree();i++)
					{
						double w_vi = 1.0 / v1->degree() * (1.0 / 4 + cos(2 * M_PI * i / v1->degree()) + 1.0 / 2 * cos(4 * M_PI * i / v1->degree()));
						new_point += verts1[i] * w_vi;
					}
					break;
				}
				}
				double w_v2 = 3.0 / 4;
				new_point += v2->point() * w_v2;
				switch (v2->degree())
				{
				case 3:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(4);//v2周围点
					int begin_index = 0, end_index = 3;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts2[0] * 5 / 12 - (verts2[1] + verts2[2]) / 12;
					break;
				}
				case 4:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(5);//v2周围点
					int begin_index = 0, end_index = 4;
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					new_point += verts2[0] * 3 / 8 - (verts2[2]) / 8;
					break;
				}
				default:
				{
					M::It::VCcwVIterator2 it1(pMesh, v2, pMesh->halfedgeSym(he));
					M::It::VClwVIterator2 it2(pMesh, v2, he);
					std::vector<CPoint> verts2(v2->degree() + 1);//v2周围点
					int begin_index = 0, end_index = v2->degree();
					while (it1 != it1.end() && it2 != it2.end() && (!(*it1)->boundary() || !(*it2)->boundary()))
					{
						if (!(*it1)->boundary())
						{
							verts2[begin_index++] = (*it1)->point();
							++it1;
						}
						if (!(*it2)->boundary())
						{
							verts2[end_index--] = (*it2)->point();
							++it2;
						}
					}
					for (int i = 0;i < v2->degree();i++)
					{
						double w_vi = 1.0 / v2->degree() * (1.0 / 4 + cos(2 * M_PI * i / v2->degree()) + 1.0 / 2 * cos(4 * M_PI * i / v2->degree()));
						new_point += verts2[i] * w_vi;
					}
					break;
				}
				}
				new_point /= 2;
			}
			while (pMesh1->idVertex(new_vid))
			{
				new_vid++;
			}
			M::V* new_midpoint_v = pMesh1->createVertex(new_vid++);

			new_midpoint_v->point() = new_point;
			edge_to_midpoint_vertex[edge] = new_midpoint_v;
			he = pMesh->halfedgeNext(he);
		} while (he != pMesh->faceHalfedge(f_old));
	}
	//遍历每个面，根据每个面所有边上的中点数量进行分类
	/*
	1.case0:面的每条边无中点，直接建面
	2.case1：1个中点，二分
	3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
	4.case3：正常1分为4
	*/
	for (auto f_old : pMesh->faces())
	{
		std::vector<M::V*> old_vertices;
		M::HE* he_start = pMesh->faceHalfedge(f_old);
		M::HE* he_iter = he_start;
		do
		{
			old_vertices.push_back(pMesh->halfedgeTarget(he_iter));
			he_iter = pMesh->halfedgeNext(he_iter);
		} while (he_iter != he_start);
		std::vector<M::E*> old_edges;
		old_edges.push_back(pMesh->vertexEdge(old_vertices[0], old_vertices[1]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[1], old_vertices[2]));
		old_edges.push_back(pMesh->vertexEdge(old_vertices[2], old_vertices[0]));
		std::vector<M::V*> midpoints;
		std::vector<int> split_edge_indices;
		for (int i = 0; i < 3; ++i)
		{
			if (edge_to_midpoint_vertex.count(old_edges[i]))
			{
				midpoints.push_back(edge_to_midpoint_vertex.at(old_edges[i]));
				split_edge_indices.push_back(i);
			}
		}
		int hanging_node_count = midpoints.size();
		switch (hanging_node_count)
		{
			//1.case0:面的每条边无中点，直接建面
		case 0:
		{
			pMesh1->createFace(std::vector<M::V*>{ old_v_to_new_v.at(old_vertices[0]), old_v_to_new_v.at(old_vertices[1]), old_v_to_new_v.at(old_vertices[2]) }, new_fid++);
			break;
		}
		//2.case1：1个中点，二分
		case 1:
		{
			M::V* m = midpoints[0];
			int edge_idx = split_edge_indices[0];
			M::V* v_start = old_v_to_new_v.at(old_vertices[edge_idx]);
			M::V* v_end = old_v_to_new_v.at(old_vertices[(edge_idx + 1) % 3]);
			M::V* v_opposite = old_v_to_new_v.at(old_vertices[(edge_idx + 2) % 3]);

			pMesh1->createFace(std::vector<M::V*>{ v_start, m, v_opposite }, new_fid++);
			pMesh1->createFace(std::vector<M::V*>{ m, v_end, v_opposite }, new_fid++);
			break;
		}
		//3.case2:两个中点，连接两个中点，剩一个四边形随便连(也可以考虑边翻转，最小角最大)
		case 2:
		{
			bool edge_flip = true;
			M::V* m1 = midpoints[0];
			M::V* m2 = midpoints[1];
			int e_idx1 = split_edge_indices[0];
			int e_idx2 = split_edge_indices[1];
			M::V* v_common_old, * v_A_old, * v_B_old;
			//保持逆时针，正方向与原来一致
			if (e_idx1 == 0 && e_idx2 == 1)
			{
				v_common_old = old_vertices[1], v_A_old = old_vertices[0], v_B_old = old_vertices[2];
			}
			else if (e_idx1 == 1 && e_idx2 == 2)
			{
				v_common_old = old_vertices[2], v_A_old = old_vertices[1], v_B_old = old_vertices[0];
			}
			else if (e_idx1 == 0 && e_idx2 == 2)
			{
				//保持逆时针，正方向与原来一致
				v_common_old = old_vertices[0], v_A_old = old_vertices[2], v_B_old = old_vertices[1];
				/*v_common_old = old_vertices[0], v_A_old = old_vertices[1], v_B_old = old_vertices[2];*/
			}

			M::V* v_common_new = old_v_to_new_v.at(v_common_old);
			M::V* v_A_new = old_v_to_new_v.at(v_A_old);
			M::V* v_B_new = old_v_to_new_v.at(v_B_old);
			M::V* m_A, * m_B;
			M::E* edge_A_common_old = pMesh->vertexEdge(v_A_old, v_common_old);

			if (edge_to_midpoint_vertex.at(edge_A_common_old) == m1)
			{
				m_A = m1; m_B = m2;
			}
			else
			{
				//e_idx1 == 0 && e_idx2 == 2 特殊情况处理
				m_A = m2; m_B = m1;
			}
			//是否加入边翻转
			if (edge_flip)
			{
				double A_r1 = get_circleRadius(m_A->point(), v_A_new->point(), m_B->point());
				double A_r2 = get_circleRadius(m_B->point(), v_A_new->point(), v_B_new->point());
				double max_A_r = std::max(A_r1, A_r2);
				double B_r1 = get_circleRadius(m_A->point(), v_B_new->point(), m_B->point());
				double B_r2 = get_circleRadius(m_A->point(), v_B_new->point(), v_A_new->point());
				double max_B_r = std::max(B_r1, B_r2);
				pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
				if (max_A_r <= max_B_r)
				{
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
				}
				else
				{
					pMesh1->createFace(std::vector<M::V*>{ v_B_new, m_A, m_B }, new_fid++);
					pMesh1->createFace(std::vector<M::V*>{ v_B_new, v_A_new, m_A  }, new_fid++);
				}
			}
			else
			{
				pMesh1->createFace(std::vector<M::V*>{ v_common_new, m_B, m_A }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_A, m_B }, new_fid++);
				pMesh1->createFace(std::vector<M::V*>{ v_A_new, m_B, v_B_new }, new_fid++);
			}
			break;
		}
		//4.case3：正常1分为4
		case 3:
		{
			M::V* v0_new = old_v_to_new_v.at(old_vertices[0]);
			M::V* v1_new = old_v_to_new_v.at(old_vertices[1]);
			M::V* v2_new = old_v_to_new_v.at(old_vertices[2]);
			M::V* m01_new = edge_to_midpoint_vertex.at(old_edges[0]);
			M::V* m12_new = edge_to_midpoint_vertex.at(old_edges[1]);
			M::V* m20_new = edge_to_midpoint_vertex.at(old_edges[2]);
			pMesh1->createFace(std::vector<M::V*>{ m01_new, m12_new, m20_new }, new_fid++);
			pMesh1->createFace(std::vector<M::V*>{ v0_new, m01_new, m20_new }, new_fid++);
			pMesh1->createFace(std::vector<M::V*>{ v1_new, m12_new, m01_new }, new_fid++);
			pMesh1->createFace(std::vector<M::V*>{ v2_new, m20_new, m12_new }, new_fid++);
			break;
		}
		}
	}
	pMesh1->labelBoundary();
}


void NUBELib::compute_vertex_degree(M* pMesh)
{
	for (auto vertex : pMesh->vertices())
	{
		int count = 0;
		for (M::It::VCcwVIterator it(pMesh, vertex);it != it.end();++it)
		{
			count++;
		}

		vertex->degree() = count;
	}
}

double NUBELib::get_div_halfedge_angle(M::HE* he)
{
	CPoint a = (he->source()->point() + he->target()->point()) / 2 - he->he_next()->target()->point();
	CPoint b = he->source()->point()- he->he_next()->target()->point();
	CPoint c = he->target()->point()- he->he_next()->target()->point();
	/*std::cout << "he->source()->point():" << he->source()->point()[0] << " " << he->source()->point()[1] << " " << he->source()->point()[2] << std::endl;
	std::cout << "he->target()->point():" << he->target()->point()[0] << " " << he->target()->point()[1] << " " << he->target()->point()[2] << std::endl;
	std::cout << "he->he_next()->target()->point():" << he->he_next()->target()->point()[0] << " " << he->he_next()->target()->point()[1] << " " << he->he_next()->target()->point()[2] << std::endl;*/
	return std::max((a * b) / a.norm() / b.norm(), (a * c) / a.norm() / c.norm());

}

double NUBELib::get_triangle_edge_length(CEdge* e)
{
	return (e->halfedge(0)->source()->point() - e->halfedge(0)->target()->point()).norm();
}
