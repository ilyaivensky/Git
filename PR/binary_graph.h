#ifndef _BINARY_GRAPH_H_
#define _BINARY_GRAPH_H_

#include <ostream>

#include "LA/matrix.h"

struct Node : public std::pair<signed, signed>
{
	typedef pair<signed, signed> Super;

	Node() : Super(-1, -1) {}

	Node(signed i, signed j) : Super(i, j) {}

	bool deleted() const { return first == -1 || second == -1; }
	bool exist() const { return !deleted(); }

	signed & i() { return first; }
	const signed & i() const { return first; }

	signed & j() { return second; }
	const signed & j() const { return second; }

	Node adjacent(unsigned code) const;

	// Returns 8 neighbours (they are not neccessarily in the graph)
	vector<Node> getNeighbours() const;

	// Returns 0-7 if adjacent, -1 otherwise
	signed getAdjacentCode(const Node & other) const;
};

ostream & operator<<(ostream & os, const Node & n);

struct Edge : public std::pair<Node, Node>
{
	typedef pair<Node, Node> Super;

	Edge () {}
	
	Edge(const Node & src, const Node & dest) :
		Super(src, dest) {}

	bool deleted() const { return first.deleted() || second.deleted(); }
	bool exist() const { return !deleted(); }

	Node & src() { return first; }
	const Node & src() const { return first; }

	Node & dest() { return second; }
	const Node & dest() const { return second; }
};

ostream & operator<<(ostream & os, const Edge & e);

struct Nodes : public std::vector<Node>
{ 
	Nodes(unsigned n = 0) : vector<Node>(n) {}

	template <class T>
	Nodes(const Matrix<T> & img)
	{
		for (signed i = 0; i < (signed)img.row; ++i)
			for (signed j = 0; j < (signed)img.col; ++j)
				if (img[i][j])
					push_back(Node(i, j));

		std::sort(this->begin(), this->end());
	}

	// Returns 8 neighbours of given node.
	// If some of them do not exist in the graph
	// then such nodes will be returned as deleted) 
	Nodes getNeighbours(const Node & node) const;

	// Batch remove
	Nodes & remove(const vector<Node> &ns);

	bool isRemoved(const Node & n);
};

struct Edges : public std::vector<Edge>
{
	Edges(unsigned n = 0) : vector<Edge>(n) {}

	template <class T>
	Edges(const Matrix<T> & img);
	
	Edges getClosedCircleEdges(const Node & node) const;

	Edges get16CircleEdges(const Node & node) const;
	Edges get20CircleEdges(const Node & node) const;

	// Batch remove
	Edges & remove(const vector<Edge> & es);
	bool isRemoved(const Edge & e);

	void add(const Edge & e);
	// Batch add
	void add(const vector<Edge> & es);
};

////////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////////

template <class T>
Edges::Edges(const Matrix<T> & img) 
{
	// Here we define direction of each edge in the graph.
	// Each two adjacent nodes are connected only by a single one-directional edge

	// Add all horizontal edges (direction left->right)
	for (signed i = 0; i < (signed)img.row; ++i)
		for (signed j = 0; j < (signed)img.col - 1; ++j)
			if (img[i][j] && img[i][j + 1])
				push_back(Edge(Node(i, j), Node(i, j + 1)));

	// Add all vertical edges (direction up->down)
	for (signed j = 0; j < (signed)img.col; ++j)
		for (signed i = 0; i < (signed)img.row - 1; ++i)
			if (img[i][j] && img[i + 1][j])
				push_back(Edge(Node(i, j), Node(i + 1, j)));

	// Add all left diagonal edges (direction [up, left]->[down,right])
	for (signed i = 0; i < (signed)img.row - 1; ++i)
		for (signed j = 0; j < (signed)img.col - 1; ++j)
			if (img[i][j] && img[i + 1][j + 1])
				push_back(Edge(Node(i, j), Node(i + 1, j + 1)));

	// Add all right diagonal edges (direction [up, right]->[down,left])
	for (signed i = 1; i < (signed)img.row; ++i)
		for (signed j = 1; j < (signed)img.col; ++j)
			if (img[i][j] && img[i + 1][j - 1])
				push_back(Edge(Node(i, j), Node(i + 1, j - 1)));

	// ... and always keep edges sorted
	std::sort(this->begin(), this->end());
}

#endif