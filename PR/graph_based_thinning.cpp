# include "PR\graph_based_thinning.h"

#include <iostream>

using namespace std;

void delete_diag_at_concaves(const Nodes & nodes, Edges & edges)
{
	vector<Edge> to_remove;

	for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
	{
		const Node & curr = *it;
		// Get neighbours from the graph (deleted nodes are represented as (-1,-1))
		Nodes n = nodes.getNeighbours(curr);
		if (!n[1].exist() && n[0].exist() && n[2].exist())
			to_remove.push_back(Edge(n[0], n[2]));
		if (!n[3].exist() && n[2].exist() && n[4].exist())
			to_remove.push_back(Edge(n[4], n[2]));
		if (!n[5].exist() && n[4].exist() && n[6].exist())
			to_remove.push_back(Edge(n[6], n[4]));
		if (!n[7].exist() && n[6].exist() && n[0].exist())
			to_remove.push_back(Edge(n[6], n[0]));
	}

	if (to_remove.size())
	{
		make_vector_set(to_remove);
		// cerr << "delete_diag_at_concaves: deleting " << to_remove << endl; 
		edges.remove(to_remove);
	}
}

void delete_border_nodes(Nodes & nodes, Edges & edges)
{
	Nodes::iterator itEnd = nodes.end();
	for (Nodes::iterator it = nodes.begin(); it != nodes.end(); )
	{
		const Node & curr = *it;
		unsigned connectivity = 0;
		Edges circleClosure = edges.getClosedCircleEdges(curr);

		for (unsigned k = 0; k < circleClosure.size(); ++k)
			if (circleClosure[k].exist())
				++connectivity;

		if (connectivity == circleClosure.size())
			connectivity = 0;
		else
		{
			connectivity = 0;
			for (unsigned k = 0; k < 8; ++k)
				if (circleClosure[k].exist())
					++connectivity;
			for (unsigned k = 0; k < 8; ++k)
				if (circleClosure[k].exist() && 
					circleClosure[k + 1].exist())
					--connectivity;
			for (unsigned k = 0; k < 4; ++k)
				if (circleClosure[2 * k].exist() && 
					circleClosure[2 * k + 1].exist() &&
					circleClosure[2 * k + 2].exist())
					++connectivity;
			for (unsigned k = 0; k < 4; ++k)
				if (circleClosure[2 * k].exist() && 
					circleClosure[2 * k + 2].exist() &&
					circleClosure[k + 9].exist())
					--connectivity;
		}

		if (connectivity == 1)
		{
			vector<Edge> edges_to_remove;
			// Remove 8 edges, which are connected to curr node
			for (unsigned k = 0; k < 8; ++k)
				if (circleClosure[k].exist())
					edges_to_remove.push_back(circleClosure[k]);

			make_vector_set(edges_to_remove);
			//cerr << "delete_border_nodes: deleting edges: " << edges_to_remove << endl;
			edges.remove(edges_to_remove);

			it = nodes.erase(it);
		}
		else ++it;
	}
}

void delete_extra_diag_edges(const Nodes & nodes, Edges & edges)
{
	vector<Edge> to_remove;

	// Delete extra diagonal edges (i.e. edges 1, 3, 5, 7)
	for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
	{
		const Node & curr = *it;
		Nodes n = nodes.getNeighbours(curr);
		Edges e = edges.get16CircleEdges(curr);
				
		if (n[1].exist() && n[2].exist())
			if (e[2].exist() && e[9].exist())
				if (!e[0].exist() || !e[8].exist())
					to_remove.push_back(e[1]);

		if (n[3].exist() && n[4].exist())
			if (e[4].exist() && e[11].exist())
				if (!e[2].exist() || !e[10].exist())
					to_remove.push_back(e[3]);

		if (n[5].exist() && n[6].exist())
			if (e[6].exist() && e[13].exist())
				if (!e[4].exist() || !e[12].exist())
					to_remove.push_back(e[5]);

		if (n[7].exist() && n[0].exist())
			if (e[0].exist() && e[15].exist())
				if (!e[6].exist() || !e[14].exist())
					to_remove.push_back(e[7]);
	}

	if (to_remove.size() > 0)
	{
		make_vector_set(to_remove);
		cerr << "delete_extra_diag_edges: deleting edges: " << to_remove << endl;
		edges.remove(to_remove);
	}
}

void delete_extra_vert_and_hor_edges(const Nodes & nodes, Edges & edges)
{
	vector<Edge> to_remove;

	for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
	{
		const Node & curr = *it;
		Nodes n = nodes.getNeighbours(curr);
		Edges e = edges.get20CircleEdges(curr);
				
		if (n[2].exist() && n[3].exist() && n[4].exist())
			if (e[4].exist() && e[17].exist() && e[3].exist())
				if ((!e[5].exist() && !e[18].exist()) ||
					(!e[5].exist() && !e[6].exist()) ||
					(!e[12].exist() && !e[18].exist()))
					to_remove.push_back(e[4]);

		if (n[4].exist() && n[5].exist() && n[6].exist())
			if (e[6].exist() && e[18].exist() && e[5].exist())
				if ((!e[7].exist() && !e[19].exist()) ||
					(!e[7].exist() && !e[8].exist()) ||
					(!e[14].exist() && !e[19].exist()))
					to_remove.push_back(e[6]);

		if (n[6].exist() && n[7].exist() && n[0].exist())
			if (e[0].exist() && e[19].exist() && e[7].exist())
				if ((!e[9].exist() && !e[16].exist()) ||
					(!e[9].exist() && !e[10].exist()) ||
					(!e[14].exist() && !e[16].exist()))
					to_remove.push_back(e[0]);

		if (n[0].exist() && n[1].exist() && n[2].exist())
			if (e[2].exist() && e[16].exist() && e[1].exist())
				if ((!e[11].exist() && !e[17].exist()) ||
					(!e[11].exist() && !e[12].exist()) ||
					(!e[14].exist() && !e[17].exist()))
					to_remove.push_back(e[2]);
	}

	if (to_remove.size() > 0)
	{
		make_vector_set(to_remove);
		cerr << "delete_extra_vert_and_hor_edges: deleting edges: " << to_remove << endl;
		edges.remove(to_remove);
	}
}

bool delete_intersections(const Nodes & nodes, Edges & edges)
{
	vector<Edge> to_add;
	vector<Edge> to_remove;

	for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
	{
		const Node & curr = *it;
		Nodes n = nodes.getNeighbours(curr);
		Edges e = edges.get20CircleEdges(curr);

		// Intersection (e[1], e[16])
		if (n[0].exist() && n[1].exist() && n[2].exist())
		{
			// Rule S1
			if (e[1].exist() && e[16].exist())
				if (!e[0].exist() && !e[2].exist() && !e[8].exist() && !e[9].exist())
				{
					to_remove.push_back(e[16]);
					to_add.push_back(Edge(n[0], n[1]));
					to_add.push_back(Edge(n[1], n[2]));
				}

			// Rule S2
			if (e[1].exist() && e[2].exist() && e[16].exist())
				if (!e[0].exist() && !e[8].exist() && !e[9].exist())
				{
					to_remove.push_back(e[1]);
					to_add.push_back(Edge(n[1], n[2]));
				}

			// Rule S3
			if (e[1].exist() && e[8].exist() && e[16].exist())
				if (!e[0].exist() && !e[2].exist() && !e[9].exist())
				{
					to_remove.push_back(e[16]);
					to_add.push_back(Edge(n[1], n[2]));
				}

			// Rule S4
			if (e[1].exist() && e[9].exist() && e[16].exist())
				if (!e[0].exist() && !e[2].exist() && !e[8].exist())
				{
					to_remove.push_back(e[16]);
					to_add.push_back(Edge(n[0], n[1]));
				}

			// Rule S5
			if (e[0].exist() && e[1].exist() && e[16].exist())
				if (!e[2].exist() && !e[8].exist() && !e[9].exist())
				{
					to_remove.push_back(e[1]);
					to_add.push_back(Edge(n[0], n[1]));
				}
		}

		// Intersection (e[3], e[17])
		if (n[2].exist() && n[3].exist() && n[4].exist())
		{
			// Rule S1
			if (e[3].exist() && e[17].exist())
				if (!e[2].exist() && !e[4].exist() && !e[10].exist() && !e[11].exist())
				{
					to_remove.push_back(e[17]);
					to_add.push_back(Edge(n[2], n[3]));
					to_add.push_back(Edge(n[4], n[3]));
				}

			// Rule S2
			if (e[3].exist() && e[4].exist() && e[17].exist())
				if (!e[2].exist() && !e[10].exist() && !e[11].exist())
				{
					to_remove.push_back(e[3]);
					to_add.push_back(Edge(n[4], n[3]));
				}

			// Rule S3
			if (e[3].exist() && e[10].exist() && e[17].exist())
				if (!e[2].exist() && !e[4].exist() && !e[11].exist())
				{
					to_remove.push_back(e[17]);
					to_add.push_back(Edge(n[4], n[3]));
				}

			// Rule S4
			if (e[3].exist() && e[11].exist() && e[17].exist())
				if (!e[2].exist() && !e[4].exist() && !e[10].exist())
				{
					to_remove.push_back(e[17]);
					to_add.push_back(Edge(n[2], n[3]));
				}

			// Rule S5
			if (e[2].exist() && e[3].exist() && e[17].exist())
				if (!e[4].exist() && !e[10].exist() && !e[11].exist())
				{
					to_remove.push_back(e[3]);
					to_add.push_back(Edge(n[2], n[3]));
				}
		}

		// Intersection (e[5], e[18])
		if (n[4].exist() && n[5].exist() && n[6].exist())
		{
			// Rule S1
			if (e[5].exist() && e[18].exist())
				if (!e[4].exist() && !e[6].exist() && !e[12].exist() && !e[13].exist())
				{
					to_remove.push_back(e[18]);
					to_add.push_back(Edge(n[5], n[4]));
					to_add.push_back(Edge(n[6], n[5]));
				}

			// Rule S2
			if (e[5].exist() && e[6].exist() && e[18].exist())
				if (!e[4].exist() && !e[12].exist() && !e[13].exist())
				{
					to_remove.push_back(e[5]);
					to_add.push_back(Edge(n[6], n[5]));
				}

			// Rule S3
			if (e[5].exist() && e[12].exist() && e[18].exist())
				if (!e[4].exist() && !e[6].exist() && !e[13].exist())
				{
					to_remove.push_back(e[18]);
					to_add.push_back(Edge(n[6], n[5]));
				}

			// Rule S4
			if (e[5].exist() && e[14].exist() && e[18].exist())
				if (!e[4].exist() && !e[6].exist() && !e[12].exist())
				{
					to_remove.push_back(e[18]);
					to_add.push_back(Edge(n[5], n[4]));
				}

			// Rule S5
			if (e[4].exist() && e[5].exist() && e[18].exist())
				if (!e[6].exist() && !e[12].exist() && !e[13].exist())
				{
					to_remove.push_back(e[5]);
					to_add.push_back(Edge(n[5], n[4]));
				}
		}

		// Intersection (e[7], e[19])
		if (n[6].exist() && n[7].exist() && n[0].exist())
		{
			// Rule S1
			if (e[7].exist() && e[19].exist())
				if (!e[6].exist() && !e[0].exist() && !e[14].exist() && !e[15].exist())
				{
					to_remove.push_back(e[19]);
					to_add.push_back(Edge(n[7], n[6]));
					to_add.push_back(Edge(n[7], n[0]));
				}

			// Rule S2
			if (e[7].exist() && e[0].exist() && e[19].exist())
				if (!e[6].exist() && !e[14].exist() && !e[15].exist())
				{
					to_remove.push_back(e[7]);
					to_add.push_back(Edge(n[7], n[0]));
				}

			// Rule S3
			if (e[7].exist() && e[14].exist() && e[19].exist())
				if (!e[6].exist() && !e[0].exist() && !e[15].exist())
				{
					to_remove.push_back(e[18]);
					to_add.push_back(Edge(n[7], n[0]));
				}

			// Rule S4
			if (e[7].exist() && e[8].exist() && e[19].exist())
				if (!e[6].exist() && !e[0].exist() && !e[14].exist())
				{
					to_remove.push_back(e[19]);
					to_add.push_back(Edge(n[7], n[6]));
				}

			// Rule S5
			if (e[6].exist() && e[7].exist() && e[19].exist())
				if (!e[8].exist() && !e[14].exist() && !e[15].exist())
				{
					to_remove.push_back(e[7]);
					to_add.push_back(Edge(n[7], n[6]));
				}
		}
	}

	make_vector_set(to_remove);
	edges.remove(to_remove);
	edges.add(to_add);

	return to_remove.size() > 0;
}