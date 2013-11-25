#ifndef _GRAPH_BASED_THINNING_H_
#define _GRAPH_BASED_THINNING_H_

#include "PR/binary_graph.h"
#include <algorithm>

void delete_diag_at_concaves(const Nodes & nodes, Edges & edges);
void delete_border_nodes(Nodes & nodes, Edges & edges);
void delete_extra_diag_edges(const Nodes & nodes, Edges & edges);
void delete_extra_vert_and_hor_edges(const Nodes & nodes, Edges & edges);

// Returns true if at least one intersection was removed
bool delete_intersections(const Nodes & nodes, Edges & edges);

template <class T>
Matrix<T> graph_based_thinning(const Matrix<T> & img)
{
	cerr << "Processing " << endl << img << endl;
	
	// Build the graph
	Nodes nodes(img);
	Edges edges(img);
	size_t graph_size = nodes.size() + edges.size();

	while (true)
	{
		// Delete diagonal edges at concave corners
		delete_diag_at_concaves(nodes, edges);
		
		// Delete border nodes
		delete_border_nodes(nodes, edges);

#if 1
		// Build a new image from the nodes left after thinning
		Matrix<T> t(img.row, img.col);
		for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
		{
			const Node & node = *it;
			t[node.i()][node.j()] = img[node.i()][node.j()]; 
		}

		cerr << "After deleting diag at concaves and deleting border nodes:" << endl << t << endl;
#endif
		size_t edges_before;
		do 
		{
			do 
			{
				edges_before = edges.size();

				// Delete extra diagonal edges (i.e. edges 1, 3, 5, 7)
				delete_extra_diag_edges(nodes, edges);

				// Delete extra horizontal and vertical edges (i.e. edges 0, 2, 4, 6)
				delete_extra_vert_and_hor_edges(nodes, edges);
			
			}  while (edges.size() < edges_before);

		// Delete intersections
		} while (delete_intersections(nodes, edges));

		// Check exit condition
		unsigned next_graph_size = nodes.size() + edges.size();
		if (next_graph_size == graph_size)
			break;
		
		graph_size = next_graph_size;
	}

	// Build a new image from the nodes left after thinning
	Matrix<T> thinned(img.row, img.col);
	for (Nodes::const_iterator it = nodes.begin(), itEnd = nodes.end(); it != itEnd; ++it)
	{
		const Node & node = *it;
		thinned[node.i()][node.j()] = img[node.i()][node.j()]; 
	}

	cerr << "Thinned:" << endl << thinned << endl;
	_getch();

	return thinned;
}

#endif