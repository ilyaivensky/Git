#include "PR/binary_graph.h"

#include <exception>
#include <algorithm>
#include <iostream>

using namespace std;

/*****************************************************************
*
* Node
*
******************************************************************/

Node Node::adjacent(unsigned code) const
{
	switch (code)
	{
	case 0:
		return Node(i(), j() - 1);
	case 1:
		return Node(i() + 1, j() - 1);
	case 2:
		return Node(i() + 1, j());
	case 3:
		return Node(i() + 1, j() + 1);
	case 4:
		return Node(i(), j() + 1);
	case 5:
		return Node(i() - 1, j() + 1);
	case 6:
		return Node(i() - 1, j());
	case 7:
		return Node(i() - 1, j() - 1);
	default:
		throw exception("Node::adjacent(): invalid code");
	}
}

vector<Node> Node::getNeighbours() const
{
	vector<Node> n(8);
	n[0] = Node(i(), j() - 1);
	n[1] = Node(i() + 1, j() - 1);
	n[2] = Node(i() + 1, j());
	n[3] = Node(i() + 1, j() + 1);
	n[4] = Node(i(), j() + 1);
	n[5] = Node(i() - 1, j() + 1);
	n[6] = Node(i() - 1, j());
	n[7] = Node(i() - 1, j() - 1);

	return n;
}

signed Node::getAdjacentCode(const Node & other) const
{
	if (other.i() == this->i())
	{
		if (other.j() == this->j() - 1)
			return 0;
		if (other.j() == this->j() + 1)
			return 4;
	}
	else if (other.i() == this->i() - 1)
	{
		if (other.j() == this->j() - 1)
			return 7;
		if (other.j() == this->j())
			return 6;
		if (other.j() == this->j() + 1)
			return 5;
	}
	else if (other.i() == this->i() + 1)
	{
		if (other.j() == this->j() - 1)
			return 1;
		if (other.j() == this->j())
			return 2;
		if (other.j() == this->j() + 1)
			return 3;
	}
	return -1;
}

ostream & operator<<(ostream & os, const Node & n)
{
	os << "(" << n.first << "," << n.second << ")";
	return os;
}

ostream & operator<<(ostream & os, const Edge & e)
{
	os << "(" << e.first << "," << e.second << ")";
	return os;
}

/********************************************************************************
*
* Nodes
*
*********************************************************************************/
	
Nodes Nodes::getNeighbours(const Node & node) const
{
	// Initializes all 8 neighbours as deleted from the graph
	Nodes result(8);

	Nodes adjacent(8);
	for (unsigned i = 0; i < 8; ++i)
		adjacent[i] = node.adjacent(i);

	const_iterator itLower = std::lower_bound(begin(), end(), Node(node.i() - 1, node.j() - 1));
	const_iterator itUpper = std::upper_bound(begin(), end(), Node(node.i() + 1, node.j() + 2));

	// Fils result with neighbours which still exist on the graph 
	for (const_iterator it = itLower; it != itUpper; ++it)
	{
		const Node & curr = *it;
		signed code = node.getAdjacentCode(curr);
		if (code != -1)
			result[code] = curr;
	}

	return result;
}

Nodes & Nodes::remove(const vector<Node> & ns)
{
	// cerr << "Removing nodes: " << ns << endl;
	Nodes diff(this->size());
	Nodes::iterator itEnd = set_difference(this->begin(), this->end(), 
		ns.begin(), ns.end(), diff.begin());
	diff.erase(itEnd, diff.end());
	swap(diff);
	return *this;
}

bool Nodes::isRemoved(const Node & n)
{
	if (n.deleted())
		return true;
	return find(begin(), end(), n) == end();
}

Edges Edges::getClosedCircleEdges(const Node & node) const
{
	Edges result(13);

	vector<Node> n = node.getNeighbours();

	Edge lowerBoundEdge(Node(node.i() - 1, node.j() - 1), Node(node.i() - 1, node.j()));
	Edge upperBoundEdge(Node(node.i() + 1, node.j() + 1), Node());

	const_iterator itLower = std::lower_bound(begin(), end(), lowerBoundEdge);
	const_iterator itUpper = std::upper_bound(begin(), end(), upperBoundEdge);

	for (const_iterator it = itLower; it != itUpper; ++it)
	{
		const Edge & e = *it;
			
		// Add all edges originating at (i,j)
		if (e.src() == node)
		{
			if (e.dest() == n[1])
				result[1] = e;
			else if (e.dest() == n[2])
				result[2] = e;
			else if (e.dest() == n[3])
				result[3] = e;
			else if (e.dest() == n[4])
				result[4] = e;
		}
		// Add all edges ending at (i,j)
		else if (e.dest() == node)
		{
			if (e.src() == n[0])
				result[0] = result[8] = e;
			else if (e.src() == n[7])
				result[7] = e;
			else if (e.src() == n[6])
				result[6] = e;
			else if (e.src() == n[5])
				result[5] = e;
		}
		// Add diagonals surrounding (i,j)	
		else if (e.src() == n[6]) 
		{
			if (e.dest() == n[0])
				result[12] = e;
			else if (e.dest() == n[4])
				result[11] = e;
		}
		else if (e.dest() == n[2]) 
		{
			if (e.src() == n[0])
				result[9] = e;
			else if (e.src() == n[4])
				result[10] = e;
		}
	}

	return result;
}

Edges Edges::get16CircleEdges(const Node & node) const
{
	Edges result(16);

	Nodes adjacent(8);
	for (unsigned i = 0; i < 8; ++i)
		adjacent[i] = node.adjacent(i);

	for (const_iterator it = begin(), itEnd = end(); it != itEnd; ++it)
	{
		const Edge & e = *it;

		if (e.src() == node)
			result[node.getAdjacentCode(e.dest())] = e;
		else if (e.dest() == node)
			result[node.getAdjacentCode(e.src())] = e;
		else if (e.src() == adjacent[0] && e.dest() == adjacent[1])
			result[8] = e;
		else if (e.src() == adjacent[1] && e.dest() == adjacent[2])
			result[9] = e;
		else if (e.src() == adjacent[2] && e.dest() == adjacent[3])
			result[10] = e;
		else if (e.src() == adjacent[4] && e.dest() == adjacent[3])
			result[11] = e;
		else if (e.src() == adjacent[5] && e.dest() == adjacent[4])
			result[12] = e;
		else if (e.src() == adjacent[6] && e.dest() == adjacent[5])
			result[13] = e;
		else if (e.src() == adjacent[7] && e.dest() == adjacent[6])
			result[14] = e;
		else if (e.src() == adjacent[7] && e.dest() == adjacent[0])
			result[15] = e;
	}

	return result;
}

Edges Edges::get20CircleEdges(const Node & node) const
{
	Edges result(20);

	Nodes adjacent(8);
	for (unsigned i = 0; i < 8; ++i)
		adjacent[i] = node.adjacent(i);

	for (const_iterator it = begin(), itEnd = end(); it != itEnd; ++it)
	{
		const Edge & e = *it;

		if (e.src() == node)
			result[node.getAdjacentCode(e.dest())] = e;
		else if (e.dest() == node)
			result[node.getAdjacentCode(e.src())] = e;
		else if (e.src() == adjacent[0] && e.dest() == adjacent[1])
			result[8] = e;
		else if (e.src() == adjacent[1] && e.dest() == adjacent[2])
			result[9] = e;
		else if (e.src() == adjacent[2] && e.dest() == adjacent[3])
			result[10] = e;
		else if (e.src() == adjacent[4] && e.dest() == adjacent[3])
			result[11] = e;
		else if (e.src() == adjacent[5] && e.dest() == adjacent[4])
			result[12] = e;
		else if (e.src() == adjacent[6] && e.dest() == adjacent[5])
			result[13] = e;
		else if (e.src() == adjacent[7] && e.dest() == adjacent[6])
			result[14] = e;
		else if (e.src() == adjacent[7] && e.dest() == adjacent[0])
			result[15] = e;
		else if (e.src() == adjacent[0] && e.dest() == adjacent[2])
			result[16] = e;
		else if (e.src() == adjacent[4] && e.dest() == adjacent[2])
			result[17] = e;
		else if (e.src() == adjacent[6] && e.dest() == adjacent[6])
			result[18] = e;
		else if (e.src() == adjacent[6] && e.dest() == adjacent[0])
			result[19] = e;
	}

	return result;
}

Edges & Edges::remove(const vector<Edge> & es)
{
	// cerr << "Removing edges: " << es << endl;
	Edges diff(this->size());
	Edges::iterator itEnd = set_difference(this->begin(), this->end(), 
		es.begin(), es.end(), diff.begin());
	diff.erase(itEnd, diff.end());
	swap(diff);
	return *this;
}

bool Edges::isRemoved(const Edge & e)
{
	if (e.deleted())
		return true;
	return find(begin(), end(), e) == end();
}

void Edges::add(const Edge & e)
{
	push_back(e);
	make_vector_set(*this);
}

void Edges::add(const vector<Edge> & es)
{
	insert(end(), es.begin(), es.end());
	make_vector_set(*this);
}