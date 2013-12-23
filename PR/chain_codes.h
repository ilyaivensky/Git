#ifndef _CHAIN_CODES_H_
#define _CHAIN_CODES_H_

#include "PR/binary_graph.h"
#include <algorithm>
#include <conio.h>

void delete_diag_at_concaves(const Nodes & nodes, Edges & edges);
vector<Node> get_border_nodes(const Nodes & nodes, const Edges & edges);

template <class Image, class Feature>
vector<Feature> chain_codes(const Image & img)
{
	cerr << "Processing " << endl << img << endl;

	vector<Feature> result(8);
	
	// Build the graph
	Nodes nodes(img);
	Edges edges(img);
	
	// Delete diagonal edges at concave corners
	delete_diag_at_concaves(nodes, edges);
	vector<Node> border = get_border_nodes(nodes, edges);

	Image b(img.row, img.col);
	for (vector<Node>::const_iterator it = border.begin(), itEnd = border.end(); it != itEnd; ++it)
	{
		const Node & node = *it;
		b[node.i()][node.j()] = img[node.i()][node.j()];

	}

	cerr << "Border nodes:" << endl << b << endl;

	return result;
}

#endif